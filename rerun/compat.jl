# Compatibility shim for the frozen Julia reproduction oracle (rerun/*.jl).
#
# The engine went R-first and data-driven (R supplies scenarios; deterministic).
# These scripts predate that and call the OLD engine API. This shim loads the
# current engine and re-exposes the old signatures as thin methods that sample
# scenarios IN JULIA (via the retained `sample_scenarios`) and delegate to the
# new data-driven functions -- so the oracle scripts run unchanged against the
# current engine. It is NOT part of the package (only rerun/*.jl include it).
#
# Note: results are the current engine's (fixed-support, dual-simplex, bounded
# states) -- a new baseline, not digit-identical to the old unseeded runs.

include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))

using Statistics

# Per-stage in-sample support (tau stages x n_scen outcomes), sampled in Julia.
_mk_support(c::HyperParams, n) =
    [[Float64.(v) for v in sample_scenarios(Int64(n), c.lambda, c.corrmat)] for _ in 1:c.tau]

# n full-horizon trajectories as an SDDP Historical scenario list.
_mk_traj(c::HyperParams, n) =
    [[(t, Float64.(sample_scenarios(Int64(1), c.lambda, c.corrmat)[1])) for t in 1:c.tau]
     for _ in 1:Int64(n)]

# --- old API restored -------------------------------------------------------

# 1-arg model builder (used by the convergence scripts' custom train loops):
# register a Julia-sampled support.
_make_model(config::HyperParams) = _make_model(config, _mk_support(config, config.n_scenarios))

# 2-arg train: sample support + one forward trajectory per iteration.
train_model(config::HyperParams, iters::Int64; kwargs...) =
    train_model(config, _mk_support(config, config.n_scenarios),
                _mk_traj(config, iters), iters; kwargs...)

# OOB simulate by trial count.
simulate_model(model::SDDP.PolicyGraph, config::HyperParams, trials::Int64) =
    simulate_model(model, config, _mk_traj(config, trials))

# Capacity duals by sample count.
simulate_cap_duals(model::SDDP.PolicyGraph, config::HyperParams, n::Int64) =
    simulate_cap_duals(model, config, _mk_traj(config, n))

# Bound validity (the rule moved to R; re-inlined here for the oracle). Returns
# the same Dict the scripts expect.
function validate_bound(model::SDDP.PolicyGraph, config::HyperParams;
                        trials::Int64 = 200, z::Float64 = 3.0)
    oob   = _mk_traj(config, trials)
    bound = SDDP.calculate_bound(model)
    sims  = SDDP.simulate(model, length(oob);
                          sampling_scheme = SDDP.Historical(oob),
                          parallel_scheme = SDDP.Serial())
    costs = Float64[sum(node[:stage_objective] for node in sim) for sim in sims]
    n  = length(costs); m = sum(costs) / n
    sd = n > 1 ? sqrt(sum((c - m)^2 for c in costs) / (n - 1)) : 0.0
    se = sd / sqrt(n)
    Dict("bound" => bound, "sim_mean" => m, "sim_se" => se,
         "margin_se" => (se > 0 ? (bound - m) / se : 0.0),
         "trials" => n, "valid" => bound <= m + z * se)
end
