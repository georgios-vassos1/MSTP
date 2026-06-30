# Regression tests for the SDDP lower-bound validity check.
#
# Run:  julia --project=. inst/../tests/julia/test_bound_validity.jl
#   or: julia tests/julia/test_bound_validity.jl   (from the package root)
#
# Two layers:
#   1. `bound_validity` — pure decision rule, deterministic, no solver. These
#      are the fast unit tests that guard against regressing the invariant
#      "a valid lower bound does not exceed E[policy cost] beyond MC noise".
#   2. `validate_bound` — end-to-end on a tiny, well-conditioned SDDP instance
#      that should NOT overshoot; a smoke/non-regression that the detector
#      integrates with real SDDP output.

using Test, LinearAlgebra

include(joinpath(@__DIR__, "..", "..", "inst", "julia", "mstp.jl"))

@testset "bound_validity (pure rule)" begin
    # Bound well below the cost sample → valid, and sits below the mean.
    r = bound_validity(90.0, fill(100.0, 50); z = 3.0)
    @test r.valid
    @test r.margin_se < 0
    @test r.mean ≈ 100.0

    # Bound far above the mean with negligible noise → INVALID (the bug).
    costs = [100.0 + 0.5*sin(i) for i in 1:200]   # mean ≈ 100, tiny spread
    r = bound_validity(130.0, costs; z = 3.0)
    @test !r.valid
    @test r.margin_se > 3.0

    # Bound above the mean but within Monte-Carlo noise → still valid.
    # 20*sin spread over 100 samples ⇒ SE ≈ 1.4, so a bound ~1 SE above the
    # mean must be accepted (this is the "negligible overshoot" regime).
    noisy = [100.0 + 20.0*sin(i) for i in 1:100]
    r = bound_validity(101.0, noisy; z = 3.0)
    @test r.valid
    @test 0 < r.margin_se < 3

    # Zero-variance sample: bound above mean is invalid, at/below is valid.
    @test !bound_validity(101.0, fill(100.0, 10)).valid
    @test  bound_validity(100.0, fill(100.0, 10)).valid
    @test  bound_validity(99.0,  fill(100.0, 10)).valid

    # Degenerate input is rejected, not silently mis-judged.
    @test_throws ArgumentError bound_validity(1.0, Float64[])
end

# ── Tiny instance builder (well-conditioned; should not overshoot) ─────────────
function tiny_config(; tau = 3, lambda = 5.0)
    nO, nD, nC = 2, 2, 2
    nOD = nO + nD
    HyperParams(
        tau              = tau,
        nOrigins         = nO,
        nDestinations    = nD,
        nCarriers        = nC,
        nSpotCarriers    = nC,
        Bids             = [[1, 2], [3, 4]],
        Winners          = [[1], [2]],
        entry_stock_0    = fill(0, nO),
        exit_stock_0     = fill(0, nD),
        exit_short_0     = fill(0, nD),
        entry_capacity   = fill(Int64(ceil(lambda * tau)), nO),
        exit_capacity    = fill(Int64(ceil(lambda * tau)), nD),
        entry_store_coef = fill(20.0, nO),
        exit_store_coef  = fill(10.0, nD),
        exit_short_coef  = fill(30.0, nD),
        transport_coef   = fill(7.0, nO * nD),
        spot_coef        = fill(10.0, nO * nD * nC * tau),
        carrier_capacity = fill(Int64(10), (nC + nC) * tau),
        lambda           = fill(Float64(lambda), nOD),
        corrmat          = Matrix{Float64}(I, nOD, nOD),
        n_scenarios      = 5,
    )
end

@testset "validate_bound (tiny SDDP, should be valid)" begin
    config = tiny_config()
    model  = train_model(config, Int64(40))
    res    = validate_bound(model, config; trials = Int64(80), z = 3.0)

    @test haskey(res, "valid") && haskey(res, "bound") && haskey(res, "margin_se")
    @test res["trials"] == 80
    @test isfinite(res["bound"]) && isfinite(res["sim_mean"])
    # Well-conditioned tiny instance: the bound must be a valid lower bound.
    @test res["valid"]
end
