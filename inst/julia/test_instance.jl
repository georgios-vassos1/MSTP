using JSON3, HDF5, Statistics, DataStructures

# Load the MSTP engine
include(joinpath(@__DIR__, "mstp.jl"))

println("=== MSTP engine loaded ===\n")

# ── Load instance ──────────────────────────────────────────────────────────────
instance_path = expanduser("~/tsproj/instances/instance_12x6x6x20_001.json")
ic = JSON3.read(read(instance_path, String))

println("Instance: tau=$(ic[:tau][1]), nOrigins=$(ic[:nOrigins][1]), nDestinations=$(ic[:nDestinations][1]), nCarriers=$(ic[:nCarriers][1])")

# ── Build config ───────────────────────────────────────────────────────────────
config = HyperParams(
    tau              = ic[:tau][1],
    nOrigins         = ic[:nOrigins][1],
    nDestinations    = ic[:nDestinations][1],
    nCarriers        = ic[:nCarriers][1],
    nSpotCarriers    = ic[:nCarriers][1],
    Bids             = [Vector{Int64}(b) for b in ic[:Bids]],
    Winners          = [Vector{Int64}(w) for w in ic[:winners]],
    entry_stock_0    = Vector{Int64}(ic[:entry_stock_0]),
    exit_stock_0     = Vector{Int64}(ic[:exit_stock_0]),
    exit_short_0     = Vector{Int64}(ic[:exit_short_0]),
    entry_capacity   = Vector{Int64}(ic[:entry_capacity]),
    exit_capacity    = Vector{Int64}(ic[:exit_capacity]),
    entry_store_coef = Vector{Float64}(ic[:entry_store_coef]),
    exit_store_coef  = Vector{Float64}(ic[:exit_store_coef]),
    exit_short_coef  = Vector{Float64}(ic[:exit_short_coef]),
    transport_coef   = Vector{Float64}(ic[:transport_coef]),
    spot_coef        = Vector{Float64}(ic[:spot_coef]),
    carrier_capacity = Vector{Int64}(ic[:carrier_capacity]),
    lambda           = [2000.0],
    corrmat          = gen_cov_mat(2, 6, 0.4),
    n_scenarios      = 10,
)

println("Config built: nLanes=$(config.nLanes), nSpotLanes=$(config.nSpotLanes)")
println("  CarrierIdx blocks: $(length(config.CarrierIdx))")
println("  from_SG lengths:   $(length.(config.from_SG))")
println("  corrmat size:      $(size(config.corrmat))")

# ── Train (small iteration count for smoke test) ───────────────────────────────
println("\nTraining (100 iterations) with HiGHS...")
model = train_model(config, 100, "highs")
println("Training complete.")

# ── Simulate ───────────────────────────────────────────────────────────────────
println("\nSimulating (100 trials)...")
results = simulate_model(model, config, 100)

obj   = results["obj"]
noise = results["noise"]

println("obj:   length=$(length(obj))  mean=$(round(mean(obj), digits=2))  std=$(round(std(obj), digits=2))")
println("noise: size=$(size(noise))")

# ── Compare against reference (same instance, full 1500-iter run) ──────────────
println("\n── Reference (1500 iters, 1000 trials) ──")
ref_path = expanduser("~/tsproj/output/12x6x6x20_10x1500_sims.h5")
h5open(ref_path, "r") do f
    ref_obj = read(f, "obj_oob10e3_001")
    ref_ksi = read(f, "ksi_oob10e3_001")
    println("ref obj: length=$(length(ref_obj))  mean=$(round(mean(ref_obj), digits=2))  std=$(round(std(ref_obj), digits=2))")
    println("ref ksi: size=$(size(ref_ksi))")
    println("\nNote: means won't match exactly (different iterations + stochastic OOB),")
    println("but should be in the same ballpark (~1.3-1.5M range).")
end

println("\n=== Smoke test passed ===")
