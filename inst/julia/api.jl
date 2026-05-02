using SDDP, HiGHS, JuMP

# ─── Config ──────────────────────────────────────────────────────────────────

# Build a HyperParams from a Dict passed in from R via JuliaCall.
# R integer vectors arrive as Vector{Int32}; all casts to Int64/Float64 are
# done here so the rest of the Julia code uses concrete types throughout.
function build_config(p::AbstractDict)
    # JuliaCall passes OrderedDict{Symbol, Any}; normalise to Dict{String, Any}
    p = Dict{String, Any}(string(k) => v for (k, v) in p)
    to_i64(x::AbstractArray)  = Vector{Int64}(Int64.(x))
    to_i64(x::Number)         = [Int64(x)]
    to_f64(x::AbstractArray)  = Vector{Float64}(Float64.(x))
    to_f64(x::Number)         = [Float64(x)]
    to_vi64(v) = Vector{Vector{Int64}}([
        isa(x, AbstractArray) ? Vector{Int64}(Int64.(collect(x))) : [Int64(x)]
        for x in v])

    lambda  = isa(p["lambda"], AbstractVector) ? to_f64(p["lambda"]) : [Float64(p["lambda"])]
    corrmat = Matrix{Float64}(p["corrmat"])

    HyperParams(
        tau              = Int64(p["tau"]),
        nOrigins         = Int64(p["nOrigins"]),
        nDestinations    = Int64(p["nDestinations"]),
        nCarriers        = Int64(p["nCarriers"]),
        nSpotCarriers    = Int64(p["nSpotCarriers"]),
        Bids             = to_vi64(p["Bids"]),
        Winners          = to_vi64(p["Winners"]),
        entry_stock_0    = to_i64(p["entry_stock_0"]),
        exit_stock_0     = to_i64(p["exit_stock_0"]),
        exit_short_0     = to_i64(p["exit_short_0"]),
        entry_capacity   = to_i64(p["entry_capacity"]),
        exit_capacity    = to_i64(p["exit_capacity"]),
        entry_store_coef = to_f64(p["entry_store_coef"]),
        exit_store_coef  = to_f64(p["exit_store_coef"]),
        exit_short_coef  = to_f64(p["exit_short_coef"]),
        transport_coef   = to_f64(p["transport_coef"]),
        spot_coef        = to_f64(p["spot_coef"]),
        carrier_capacity = to_i64(p["carrier_capacity"]),
        lambda           = lambda,
        corrmat          = corrmat,
        n_scenarios      = haskey(p, "n_scenarios") ? Int64(p["n_scenarios"]) : 10,
    )
end

# ─── Train ────────────────────────────────────────────────────────────────────

# Shared HiGHS optimizer settings.  Primal simplex (strategy=4) avoids the
# OPTIMAL+INFEASIBLE_POINT status that dual simplex can emit when cut
# coefficients span many orders of magnitude.  Scale strategy 4 (max-value
# scaling) normalises the LP matrix before each solve.
# dual_feasibility_tolerance tightened to 1e-9 so that LP duals with
# |value| < 1e-9 are rounded to zero before being used as Benders cut
# coefficients — this prevents the 1e-12 "ghost" slopes that cause the
# lower bound to overshoot the upper bound on long horizons (τ=26, 52).
function _make_optimizer()
    optimizer_with_attributes(
        HiGHS.Optimizer,
        "output_flag"                  => false,
        "simplex_scale_strategy"       => 4,
        "simplex_strategy"             => 4,
        "primal_feasibility_tolerance" => 1e-7,
        "dual_feasibility_tolerance"   => 1e-9,
    )
end

# Rebuild a fresh PolicyGraph with the given config and optimizer.
function _make_model(config::HyperParams)
    SDDP.PolicyGraph(
        (sp, stage) -> transportation_t(sp, stage; config = config),
        SDDP.LinearGraph(config.tau);
        sense       = :Min,
        lower_bound = 0.0,
        optimizer   = _make_optimizer(),
    )
end

# Train the SDDP policy graph using a crash-resilient batched loop.
#
# Strategy
# --------
# Each batch of `batch_size` iterations is wrapped in a try-catch.  After
# every successful batch the Benders cuts are checkpointed to a temp file.
# If HiGHS returns OPTIMAL+INFEASIBLE_POINT and SDDP throws, we:
#   1. Rebuild a fresh PolicyGraph (clears all JuMP state).
#   2. Reload cuts from the last good checkpoint (preserving convergence).
#   3. Retry the *same* batch up to MAX_RETRIES times — new RNG state usually
#      yields a different (better-conditioned) scenario draw that avoids crash.
#   4. If MAX_RETRIES exhausted, discard the checkpoint (cut coefficients may
#      be numerically bad), rebuild completely fresh, and advance past the
#      failed iteration to break the retry loop.
#
# cut_deletion_minimum = 10 aggressively prunes dominated cuts, keeping the
# LP matrix small and reducing the chance of near-zero coefficient buildup.
# This is important for capacity-optimisation runs where capacity changes
# between gradient steps can make older cuts near-degenerate.
#
# Returns the trained PolicyGraph — held as a Julia-side proxy in R.
function train_model(config::HyperParams, iterations::Int64;
                     batch_size::Int64 = 1,
                     max_retries::Int64 = 3)
    model       = _make_model(config)
    ckpt        = tempname() * ".cuts.json"
    completed   = 0
    has_ckpt    = false
    retries     = 0

    while completed < iterations
        to_do = min(batch_size, iterations - completed)
        try
            SDDP.train(model;
                iteration_limit      = to_do,
                cut_type             = SDDP.SINGLE_CUT,
                parallel_scheme      = SDDP.Serial(),
                cut_deletion_minimum = 10,
                add_to_existing_cuts = completed > 0,
            )
            completed += to_do
            SDDP.write_cuts_to_file(model, ckpt)
            has_ckpt = true
            retries  = 0   # reset on success
        catch e
            retries += 1
            if retries <= max_retries
                @warn "SDDP numerical crash at completed=$completed (retry $retries/$max_retries): $(sprint(showerror, e)). Rebuilding from checkpoint."
                model = _make_model(config)
                if has_ckpt
                    SDDP.read_cuts_from_file(model, ckpt)
                end
                # Do NOT advance `completed` — retry same batch with new RNG state.
            else
                @warn "SDDP numerical crash at completed=$completed: max retries ($max_retries) exhausted. Discarding checkpoint and advancing."
                model    = _make_model(config)   # completely fresh, no cuts
                has_ckpt = false
                completed += to_do               # skip this batch to break loop
                retries  = 0
            end
        end
    end

    has_ckpt && isfile(ckpt) && rm(ckpt; force = true)
    model
end

# ─── Warm-start helpers ───────────────────────────────────────────────────────

# Persist Benders cuts to a JSON file so they can be reloaded into a new model.
function write_cuts(model::SDDP.PolicyGraph, path::String)
    SDDP.write_cuts_to_file(model, path)
    nothing
end

# Build a new model with updated config, seed it with cuts from a previous run,
# then continue training.  Cuts are valid across capacity changes because they
# are functions of the inventory state variables, not of carrier_capacity.
function train_model_warm(config::HyperParams, iterations::Int64, cuts_path::String)
    model = _make_model(config)
    SDDP.read_cuts_from_file(model, cuts_path)
    SDDP.train(model;
        iteration_limit      = iterations,
        cut_type             = SDDP.SINGLE_CUT,
        parallel_scheme      = SDDP.Serial(),
        cut_deletion_minimum = 100,
        add_to_existing_cuts = true,
    )
    model
end

# ─── Bound ────────────────────────────────────────────────────────────────────

# Return the SDDP lower bound (dual / Benders) after training.
function get_bound(model::SDDP.PolicyGraph)
    SDDP.calculate_bound(model)
end

# ─── Simulate ─────────────────────────────────────────────────────────────────

# Run out-of-bag simulations using the trained model.
# Generates fresh OOB scenarios from the uncertainty model and returns
# a Dict with plain arrays that JuliaCall converts to R-native types.
function simulate_model(model::SDDP.PolicyGraph, config::HyperParams, trials::Int64)
    # Generate OOB scenarios: one trajectory per trial
    oob = [
        [(t, sample_scenarios(1, config.lambda, config.corrmat)...) for t in 1:config.tau]
        for _ in 1:trials
    ]

    sims = SDDP.simulate(
        model, trials,
        [:move, :inflow, :outflow, :entry, :exitp, :exitm];
        sampling_scheme = SDDP.Historical(oob),
        parallel_scheme = SDDP.Serial(),
    )

    # Total cost per simulation (length = trials)
    obj = Float64[sum(node[:stage_objective] for node in sim) for sim in sims]

    # Realized noise: (trials × tau) × (nOrigins + nDestinations)
    noise = vcat([
        vcat([collect(ξ)' for ξ in map(node -> node[:noise_term], sim)]...)
        for sim in sims
    ]...)

    # Entry inventory .in at each stage: (trials × tau) × nOrigins
    entry = vcat([
        vcat([[node[:entry][i].in for i in 1:config.nOrigins]' for node in sim]...)
        for sim in sims
    ]...)

    # Exit stock .in at each stage: (trials × tau) × nDestinations
    exitp = vcat([
        vcat([[node[:exitp][j].in for j in 1:config.nDestinations]' for node in sim]...)
        for sim in sims
    ]...)

    # Exit shortage .in at each stage: (trials × tau) × nDestinations
    exitm = vcat([
        vcat([[node[:exitm][j].in for j in 1:config.nDestinations]' for node in sim]...)
        for sim in sims
    ]...)

    # Move allocations: (trials × tau) × (nLanes + nSpotLanes)
    moves = vcat([
        vcat([[node[:move][k] for k in 1:(config.nLanes + config.nSpotLanes)]' for node in sim]...)
        for sim in sims
    ]...)

    Dict(
        "obj"           => obj,
        "noise"         => Float64.(noise),
        "entry"         => Float64.(entry),
        "exitp"         => Float64.(exitp),
        "exitm"         => Float64.(exitm),
        "moves"         => Float64.(moves),
        "trials"        => trials,
        "tau"           => config.tau,
        "nOrigins"      => config.nOrigins,
        "nDestinations" => config.nDestinations,
        "nLanes"        => config.nLanes,
        "nSpotLanes"    => config.nSpotLanes,
    )
end

# ─── Capacity duals ───────────────────────────────────────────────────────────

# Simulate n_samples trajectories and collect the dual of each carrier capacity
# constraint at every stage.  Returns a flat Float64 vector of length
# (nCarriers + nSpotCarriers) * tau with the same (k-1)*tau + stage indexing as
# carrier_capacity, so it can be used directly as ∂V/∂carrier_capacity.
#
# Uses the default InSampleMonteCarlo sampling scheme (same noise distribution as
# training) and runs sequentially to avoid thread-safety issues with dual().
function simulate_cap_duals(model::SDDP.PolicyGraph, config::HyperParams, n_samples::Int64)
    nK    = config.nCarriers + config.nSpotCarriers
    duals = zeros(Float64, nK * config.tau)

    sims = SDDP.simulate(
        model, n_samples,
        custom_recorders = Dict{Symbol, Function}(
            :cap_duals => sp -> Float64[JuMP.dual(sp[:cap][k]) for k in 1:nK]
        ),
    )

    for sim in sims
        for (t, node) in enumerate(sim)
            for k in 1:nK
                duals[(k - 1) * config.tau + t] += node[:cap_duals][k]
            end
        end
    end
    duals ./= n_samples

    Dict("duals" => duals, "nK" => Int64(nK), "tau" => Int64(config.tau))
end

# ─── Batch train + simulate ───────────────────────────────────────────────────

# Convenience: train and simulate a list of configs in one call.
# Returns a Vector of result Dicts (one per instance).
function batch_run(configs::Vector{HyperParams}, iterations::Int64, trials::Int64)
    map(configs) do config
        model = train_model(config, iterations)
        simulate_model(model, config, trials)
    end
end
