# Pure Julia (no JuliaCall ⇒ no libomp clash): load the seed-42 instance JSON,
# train the SDDP policy, and assert the lower bound is numerically valid.
# Exits non-zero if invalid. This is the RED acceptance gate for reviewer #3.
#
# Usage:  julia tests/acceptance/check_bound_acceptance.jl [in.json]
using JSON3, LinearAlgebra, Statistics
include(joinpath(@__DIR__, "..", "..", "inst", "julia", "mstp.jl"))

path = length(ARGS) >= 1 ? ARGS[1] : "/tmp/inst52_seed42.json"
ic   = JSON3.read(read(path, String))
toi(x) = x isa Number ? Int64[Int64(x)]     : Vector{Int64}(Int64.(collect(x)))
tof(x) = x isa Number ? Float64[Float64(x)] : Vector{Float64}(Float64.(collect(x)))

cfg = HyperParams(
    tau = Int64(ic.tau), nOrigins = Int64(ic.nOrigins), nDestinations = Int64(ic.nDestinations),
    nCarriers = Int64(ic.nCarriers), nSpotCarriers = Int64(ic.nSpotCarriers),
    Bids = [toi(b) for b in ic.Bids], Winners = [toi(w) for w in ic.Winners],
    entry_stock_0 = toi(ic.entry_stock_0), exit_stock_0 = toi(ic.exit_stock_0), exit_short_0 = toi(ic.exit_short_0),
    entry_capacity = toi(ic.entry_capacity), exit_capacity = toi(ic.exit_capacity),
    entry_store_coef = tof(ic.entry_store_coef), exit_store_coef = tof(ic.exit_store_coef), exit_short_coef = tof(ic.exit_short_coef),
    transport_coef = tof(ic.transport_coef), spot_coef = tof(ic.spot_coef), carrier_capacity = toi(ic.carrier_capacity),
    lambda = fill(700.0, Int64(ic.nOrigins) + Int64(ic.nDestinations)),
    corrmat = gen_cov_mat(2, Int64(ic.nOrigins), 0.0), n_scenarios = 10)

model = train_model(cfg, Int64(300))
r     = validate_bound(model, cfg; trials = Int64(200), z = 3.0)
println("bound=$(round(r["bound"])) sim_mean=$(round(r["sim_mean"])) se=$(round(r["sim_se"])) " *
        "margin_se=$(round(r["margin_se"], digits=2)) valid=$(r["valid"])")
exit(r["valid"] ? 0 : 1)
