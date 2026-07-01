using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics
include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))
const INST=expanduser("~/drayage/MSTP/rerun/inst")
toi(x)= x isa Number ? Int64[Int64(x)] : Vector{Int64}(Int64.(collect(x)))
tof(x)= x isa Number ? Float64[Float64(x)] : Vector{Float64}(Float64.(collect(x)))
function loadcfg(name)
  ic=JSON3.read(read(joinpath(INST,"$name.json"),String)); nO=Int64(ic.nOrigins); nD=Int64(ic.nDestinations)
  HyperParams(tau=Int64(ic.tau),nOrigins=nO,nDestinations=nD,nCarriers=Int64(ic.nCarriers),nSpotCarriers=Int64(ic.nSpotCarriers),
    Bids=[toi(b) for b in ic.Bids],Winners=[toi(w) for w in ic.Winners],
    entry_stock_0=toi(ic.entry_stock_0),exit_stock_0=toi(ic.exit_stock_0),exit_short_0=toi(ic.exit_short_0),
    entry_capacity=toi(ic.entry_capacity),exit_capacity=toi(ic.exit_capacity),
    entry_store_coef=tof(ic.entry_store_coef),exit_store_coef=tof(ic.exit_store_coef),exit_short_coef=tof(ic.exit_short_coef),
    transport_coef=tof(ic.transport_coef),spot_coef=tof(ic.spot_coef),carrier_capacity=toi(ic.carrier_capacity),
    lambda=fill(700.0,nO+nD),corrmat=gen_cov_mat(2,nO,0.0),n_scenarios=20)
end
for (nm,it) in (("scal_20_42",50),("scal_40_42",20))
  c=loadcfg(nm); tt=@elapsed (m=train_model(c,Int64(it))); ts=@elapsed (r=validate_bound(m,c;trials=Int64(50),z=3.0))
  sd=r["sim_se"]*sqrt(r["trials"])
  println(">>> $nm iters=$it t_train=$(round(tt,digits=1))s t_iter=$(round(tt/it,digits=2))s t_sim=$(round(ts,digits=1))s LB=$(round(r["bound"])) UB=$(round(r["sim_mean"]))±$(round(sd)) gap%=$(round(100*(r["sim_mean"]-r["bound"])/abs(r["bound"]),digits=2)) valid=$(r["valid"])"); flush(stdout)
end
println("DONE")
