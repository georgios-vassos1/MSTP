using SDDP, HiGHS

# ─── Config ──────────────────────────────────────────────────────────────────

# Build a HyperParams from a Dict passed in from R via JuliaCall.
# R integer vectors arrive as Vector{Int32}; all casts to Int64/Float64 are
# done here so the rest of the Julia code uses concrete types throughout.
function build_config(p::Dict)
    to_i64(x)  = Int64.(x)
    to_f64(x)  = Float64.(x)
    to_vi64(v) = [Int64.(x) for x in v]   # Vector{Any} of int vecs → Vector{Vector{Int64}}

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

# Build and train the SDDP policy graph.
# Returns the trained PolicyGraph — held as a Julia-side proxy in R.
function train_model(config::HyperParams, iterations::Int64)
    optimizer = HiGHS.Optimizer

    model = SDDP.PolicyGraph(
        (sp, stage) -> transportation_t(sp, stage; config = config),
        SDDP.LinearGraph(config.tau);
        sense       = :Min,
        lower_bound = 0.0,
        optimizer   = optimizer,
    )

    SDDP.train(model;
        iteration_limit = iterations,
        cut_type        = SDDP.SINGLE_CUT,
        parallel_scheme = SDDP.Threaded(),
    )

    model
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
        parallel_scheme = SDDP.Threaded(),
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

# ─── Batch train + simulate ───────────────────────────────────────────────────

# Convenience: train and simulate a list of configs in one call.
# Returns a Vector of result Dicts (one per instance).
function batch_run(configs::Vector{HyperParams}, iterations::Int64, trials::Int64)
    map(configs) do config
        model = train_model(config, iterations)
        simulate_model(model, config, trials)
    end
end
