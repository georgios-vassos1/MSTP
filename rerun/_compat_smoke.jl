# Smoke test for rerun/compat.jl: every restored old-API method resolves and
# runs against the current engine on the tiny t2x2 instance.
using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics
include(expanduser("~/drayage/MSTP/rerun/compat.jl"))
const INST = expanduser("~/drayage/MSTP/rerun/inst")
toi(x)= x isa Number ? Int64[Int64(x)] : Vector{Int64}(Int64.(collect(x)))
tof(x)= x isa Number ? Float64[Float64(x)] : Vector{Float64}(Float64.(collect(x)))
function loadcfg(name; lambda=5.5, rho=0.4, n_scen=10)
  ic=JSON3.read(read(joinpath(INST,"$name.json"),String)); nO=Int64(ic.nOrigins); nD=Int64(ic.nDestinations)
  HyperParams(tau=Int64(ic.tau),nOrigins=nO,nDestinations=nD,nCarriers=Int64(ic.nCarriers),nSpotCarriers=Int64(ic.nSpotCarriers),
    Bids=[toi(b) for b in ic.Bids],Winners=[toi(w) for w in ic.Winners],
    entry_stock_0=toi(ic.entry_stock_0),exit_stock_0=toi(ic.exit_stock_0),exit_short_0=toi(ic.exit_short_0),
    entry_capacity=toi(ic.entry_capacity),exit_capacity=toi(ic.exit_capacity),
    entry_store_coef=tof(ic.entry_store_coef),exit_store_coef=tof(ic.exit_store_coef),exit_short_coef=tof(ic.exit_short_coef),
    transport_coef=tof(ic.transport_coef),spot_coef=tof(ic.spot_coef),carrier_capacity=toi(ic.carrier_capacity),
    lambda=fill(Float64(lambda),nO+nD),corrmat=gen_cov_mat(2,max(nO,nD),rho)[1:(nO+nD),1:(nO+nD)],n_scenarios=n_scen)
end

cfg = loadcfg("t2x2")
m1  = _make_model(cfg);                 println("PASS _make_model(cfg) 1-arg")
SDDP.train(m1; iteration_limit=3, print_level=0); println("PASS SDDP.train on _make_model")
m2  = train_model(cfg, Int64(10));      println("PASS train_model(cfg, iters) LB=", round(SDDP.calculate_bound(m2), digits=1))
vb  = validate_bound(m2, cfg; trials=Int64(20)); println("PASS validate_bound valid=", vb["valid"], " margin=", round(vb["margin_se"],digits=2))
s   = simulate_model(m2, cfg, Int64(20)); println("PASS simulate_model(...,trials) meanUB=", round(mean(s["obj"]),digits=1))
d   = simulate_cap_duals(m2, cfg, Int64(20)); println("PASS simulate_cap_duals(...,n) len=", length(d["duals"]))
println(">>> COMPAT SMOKE OK")
