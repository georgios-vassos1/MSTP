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

# ─── Noise reshaping (R supplies flat arrays) ──────────────────────────────────
#
# R owns all scenario sampling and passes the draws as flat Float64 vectors plus
# their dimensions; these helpers rebuild the nested structures SDDP needs. Index
# conventions (must match the R flatteners in R/engine.R):
#   support:   stage-major, then scenario, then dim
#     idx(t,s,d) = ((t-1)*n_scen + (s-1))*dim + d
#   trajectory: path-major, then stage, then dim
#     idx(p,t,d) = ((p-1)*tau + (t-1))*dim + d

# Per-stage in-sample support: Vector over stages of Vector of outcome vectors.
function build_support(flat::Vector{Float64}, tau::Int64, n_scen::Int64, dim::Int64)
    @assert length(flat) == tau * n_scen * dim "support flat length mismatch"
    [[[flat[((t - 1) * n_scen + (s - 1)) * dim + d] for d in 1:dim]
      for s in 1:n_scen] for t in 1:tau]
end

# Full-horizon trajectories as an SDDP Historical scenario list: Vector over
# paths of Vector of (stage, noise-vector) tuples.
function build_historical(flat::Vector{Float64}, n_traj::Int64, tau::Int64, dim::Int64)
    @assert length(flat) == n_traj * tau * dim "trajectory flat length mismatch"
    [[(t, [flat[((p - 1) * tau + (t - 1)) * dim + d] for d in 1:dim])
      for t in 1:tau] for p in 1:n_traj]
end

# ─── Train ────────────────────────────────────────────────────────────────────

# Shared HiGHS optimizer settings.  Primal simplex (strategy=4) avoids the
# OPTIMAL+INFEASIBLE_POINT status that dual simplex can emit when cut
# coefficients span many orders of magnitude.  Scale strategy 4 (max-value
# scaling) normalises the LP matrix before each solve.
#
# NOTE on dual_feasibility_tolerance: this is the simplex dual-feasibility
# (optimality) tolerance; it is NOT a post-solve rounding of the dual values
# returned by JuMP.dual.  It was previously assumed to zero out |dual| < 1e-9
# before those duals became Benders cut slopes, but the numerical-stability
# reports in the SDDP logs still show cut coefficients down to ~2e-15 and an
# explicit "numerical stability issues detected" warning (e.g. 27_run3.log),
# and the lower bound still overshoots the simulated upper bound by a
# statistically significant margin on long horizons / large instances
# (τ=52: LB−UB ≈ 7 standard errors).  The ghost-slope overshoot is therefore
# NOT fully mitigated by this setting.  `validate_bound` (below) detects the
# residual overshoot so an invalid bound is never reported silently.
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

# Build a fresh PolicyGraph over the caller-supplied per-stage support.
# `support[t]` is the finite noise support for stage t (see build_support).
function _make_model(config::HyperParams, support::Vector{Vector{Vector{Float64}}})
    length(support) == config.tau ||
        throw(ArgumentError("support has $(length(support)) stages, expected tau=$(config.tau)"))
    SDDP.PolicyGraph(
        (sp, stage) -> transportation_t(sp, stage; config = config,
                                        stage_support = support[stage]),
        SDDP.LinearGraph(config.tau);
        sense       = :Min,
        lower_bound = 0.0,
        optimizer   = _make_optimizer(),
    )
end

# Train the SDDP policy graph deterministically.
#
# The engine performs no RNG: the forward pass replays the caller-supplied
# `forward` trajectories via the sequential Historical scheme, and the backward
# pass (default CompleteSampler) enumerates the caller-supplied per-stage
# support. Given the same support and trajectories the trained cuts — and hence
# the lower bound — are reproducible.
#
# The earlier RNG-driven crash-retry loop is gone: with a fixed support a retry
# reproduces the same subproblems, so it could not escape a crash. The numerical
# crashes it worked around were caused by random initial inventory and many
# near-identical carriers per lane, both eliminated by the corrected generator
# (see CHANGES-sddp-bound-validity.md). A residual crash is surfaced, not hidden.
#
# Returns the trained PolicyGraph — held as a Julia-side proxy in R.
function train_model(config::HyperParams,
                     support::Vector{Vector{Vector{Float64}}},
                     forward::Vector{<:Vector{<:Tuple}},
                     iterations::Int64;
                     cut_deletion_minimum::Int64 = 10,
                     duality_handler = SDDP.ContinuousConicDuality())
    model = _make_model(config, support)
    try
        SDDP.train(model;
            iteration_limit      = iterations,
            sampling_scheme      = SDDP.Historical(forward),
            cut_type             = SDDP.SINGLE_CUT,
            parallel_scheme      = SDDP.Serial(),
            cut_deletion_minimum = cut_deletion_minimum,
            duality_handler      = duality_handler,
        )
    catch e
        error("SDDP training failed: $(sprint(showerror, e)). With a fixed " *
              "support this is deterministic; check instance conditioning " *
              "(see CHANGES-sddp-bound-validity.md).")
    end
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
function train_model_warm(config::HyperParams,
                          support::Vector{Vector{Vector{Float64}}},
                          forward::Vector{<:Vector{<:Tuple}},
                          iterations::Int64, cuts_path::String)
    model = _make_model(config, support)
    SDDP.read_cuts_from_file(model, cuts_path)
    SDDP.train(model;
        iteration_limit      = iterations,
        sampling_scheme      = SDDP.Historical(forward),
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

# ─── Bound validity check ─────────────────────────────────────────────────────
#
# A valid SDDP lower bound for a minimisation satisfies  bound ≤ E[policy cost].
# We estimate E[policy cost] by an in-sample Monte-Carlo simulation (the same
# noise distribution as training, via the default InSampleMonteCarlo scheme) and
# flag the bound as INVALID when it exceeds the simulated mean by more than `z`
# standard errors — i.e. by more than Monte-Carlo noise can explain.
#
# This is a detector, not a repair: it surfaces the ghost-cut overshoot
# (LB > UB) documented above for long horizons / large instances, where
# degenerate, ill-conditioned LP duals produce numerically invalid Benders cuts.
# `margin_se` is how many standard errors the bound sits ABOVE the mean
# (positive ⇒ overshoot); `valid` is false once that exceeds `z`.

# Pure decision rule (no solver, no SDDP): given the lower bound and a sample of
# realised policy costs, decide whether bound ≤ E[cost] within Monte-Carlo noise.
# Factored out so the invariant can be unit-tested deterministically.
function bound_validity(bound::Real, costs::AbstractVector{<:Real}; z::Real = 3.0)
    n  = length(costs)
    n == 0 && throw(ArgumentError("costs must be non-empty"))
    m  = sum(costs) / n
    sd = n > 1 ? sqrt(sum((c - m)^2 for c in costs) / (n - 1)) : 0.0
    se = sd / sqrt(n)
    margin_se = se > 0 ? (bound - m) / se : (bound > m ? Inf : (bound < m ? -Inf : 0.0))
    (mean = m, se = se, margin_se = margin_se, n = n, valid = bound <= m + z * se)
end

function validate_bound(model::SDDP.PolicyGraph,
                        oob::Vector{<:Vector{<:Tuple}}; z::Float64 = 3.0)
    bound = SDDP.calculate_bound(model)
    sims  = SDDP.simulate(model, length(oob);
                          sampling_scheme = SDDP.Historical(oob),
                          parallel_scheme = SDDP.Serial())
    costs = Float64[sum(node[:stage_objective] for node in sim) for sim in sims]
    r = bound_validity(bound, costs; z = z)
    if !r.valid
        @warn "SDDP lower bound exceeds in-sample simulated mean by $(round(r.margin_se, digits=2)) SE " *
              "(bound=$(round(bound, digits=1)), sim_mean=$(round(r.mean, digits=1)) ± $(round(r.se, digits=1)) SE): " *
              "bound is numerically invalid (ghost-cut overshoot)."
    end
    Dict(
        "bound"     => bound,
        "sim_mean"  => r.mean,
        "sim_se"    => r.se,
        "margin_se" => r.margin_se,
        "trials"    => r.n,
        "valid"     => r.valid,
    )
end

# ─── Simulate ─────────────────────────────────────────────────────────────────

# Run out-of-bag simulations using the trained model.
# Generates fresh OOB scenarios from the uncertainty model and returns
# a Dict with plain arrays that JuliaCall converts to R-native types.
function simulate_model(model::SDDP.PolicyGraph, config::HyperParams,
                        oob::Vector{<:Vector{<:Tuple}})
    # Out-of-sample trajectories are supplied by the caller (sampled in R).
    trials = length(oob)

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

# Replay the caller-supplied out-of-sample trajectories and collect the dual of
# each carrier capacity constraint at every stage.  Returns a flat Float64 vector
# of length (nCarriers + nSpotCarriers) * tau with the same (k-1)*tau + stage
# indexing as carrier_capacity, so it can be used directly as ∂V/∂carrier_capacity.
#
# Runs sequentially to avoid thread-safety issues with dual().
function simulate_cap_duals(model::SDDP.PolicyGraph, config::HyperParams,
                            oob::Vector{<:Vector{<:Tuple}})
    nK        = config.nCarriers + config.nSpotCarriers
    n_samples = length(oob)
    duals     = zeros(Float64, nK * config.tau)

    sims = SDDP.simulate(
        model, n_samples,
        custom_recorders = Dict{Symbol, Function}(
            :cap_duals => sp -> Float64[JuMP.dual(sp[:cap][k]) for k in 1:nK]
        );
        sampling_scheme = SDDP.Historical(oob),
        parallel_scheme = SDDP.Serial(),
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

# ─── Flat-array bridge (R boundary) ─────────────────────────────────────────────
#
# R passes scenario draws as flat Float64 vectors plus their dimensions; these
# shims reshape (build_support / build_historical) and delegate to the typed core
# functions above. They exist only to keep the R↔Julia marshaling in one place.

function train_flat(config::HyperParams,
                    support_flat, tau, n_scen, dim,
                    forward_flat, n_fwd, iterations;
                    cut_deletion_minimum::Int64 = 10)
    tau = Int64(tau); n_scen = Int64(n_scen); dim = Int64(dim); n_fwd = Int64(n_fwd)
    support = build_support(Vector{Float64}(support_flat), tau, n_scen, dim)
    forward = build_historical(Vector{Float64}(forward_flat), n_fwd, tau, dim)
    train_model(config, support, forward, Int64(iterations);
                cut_deletion_minimum = cut_deletion_minimum)
end

function train_warm_flat(config::HyperParams,
                         support_flat, tau, n_scen, dim,
                         forward_flat, n_fwd, iterations, cuts_path::String)
    tau = Int64(tau); n_scen = Int64(n_scen); dim = Int64(dim); n_fwd = Int64(n_fwd)
    support = build_support(Vector{Float64}(support_flat), tau, n_scen, dim)
    forward = build_historical(Vector{Float64}(forward_flat), n_fwd, tau, dim)
    train_model_warm(config, support, forward, Int64(iterations), cuts_path)
end

function simulate_flat(model::SDDP.PolicyGraph, config::HyperParams,
                       oob_flat, n_oob, tau, dim)
    oob = build_historical(Vector{Float64}(oob_flat), Int64(n_oob), Int64(tau), Int64(dim))
    simulate_model(model, config, oob)
end

function validate_flat(model::SDDP.PolicyGraph, oob_flat, n_oob, tau, dim; z::Float64 = 3.0)
    oob = build_historical(Vector{Float64}(oob_flat), Int64(n_oob), Int64(tau), Int64(dim))
    validate_bound(model, oob; z = z)
end

function cap_duals_flat(model::SDDP.PolicyGraph, config::HyperParams,
                        oob_flat, n_oob, tau, dim)
    oob = build_historical(Vector{Float64}(oob_flat), Int64(n_oob), Int64(tau), Int64(dim))
    simulate_cap_duals(model, config, oob)
end
