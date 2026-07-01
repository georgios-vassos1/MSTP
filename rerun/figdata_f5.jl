using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics
include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))
const INST=expanduser("~/drayage/MSTP/rerun/inst"); const FD=expanduser("~/drayage/MSTP/rerun/figdata")
toi(x)= x isa Number ? Int64[Int64(x)] : Vector{Int64}(Int64.(collect(x)))
tof(x)= x isa Number ? Float64[Float64(x)] : Vector{Float64}(Float64.(collect(x)))
function loadcfg(name; lambda, rho, n_scen=20)
  ic=JSON3.read(read(joinpath(INST,"$name.json"),String)); nO=Int64(ic.nOrigins); nD=Int64(ic.nDestinations)
  HyperParams(tau=Int64(ic.tau),nOrigins=nO,nDestinations=nD,nCarriers=Int64(ic.nCarriers),nSpotCarriers=Int64(ic.nSpotCarriers),
    Bids=[toi(b) for b in ic.Bids],Winners=[toi(w) for w in ic.Winners],
    entry_stock_0=toi(ic.entry_stock_0),exit_stock_0=toi(ic.exit_stock_0),exit_short_0=toi(ic.exit_short_0),
    entry_capacity=toi(ic.entry_capacity),exit_capacity=toi(ic.exit_capacity),
    entry_store_coef=tof(ic.entry_store_coef),exit_store_coef=tof(ic.exit_store_coef),exit_short_coef=tof(ic.exit_short_coef),
    transport_coef=tof(ic.transport_coef),spot_coef=tof(ic.spot_coef),carrier_capacity=toi(ic.carrier_capacity),
    lambda=fill(Float64(lambda),nO+nD),corrmat=gen_cov_mat(2,nO,rho),n_scenarios=n_scen)
end
mkcap(cfg,x)=HyperParams(tau=cfg.tau,nOrigins=cfg.nOrigins,nDestinations=cfg.nDestinations,nCarriers=cfg.nCarriers,nSpotCarriers=cfg.nSpotCarriers,
  Bids=cfg.Bids,Winners=cfg.Winners,entry_stock_0=cfg.entry_stock_0,exit_stock_0=cfg.exit_stock_0,exit_short_0=cfg.exit_short_0,
  entry_capacity=cfg.entry_capacity,exit_capacity=cfg.exit_capacity,entry_store_coef=cfg.entry_store_coef,exit_store_coef=cfg.exit_store_coef,
  exit_short_coef=cfg.exit_short_coef,transport_coef=cfg.transport_coef,spot_coef=cfg.spot_coef,
  carrier_capacity=Int64.(round.(x)),lambda=cfg.lambda,corrmat=cfg.corrmat,n_scenarios=cfg.n_scenarios)
function traj(cfg; outer=15, a0=20.0)
  v=0.2*mean(cfg.transport_coef); vv=fill(v,length(cfg.carrier_capacity))
  x=Float64.(cfg.carrier_capacity); xmax=copy(x)
  its=Int[]; objs=Float64[]; gns=Float64[]; steps=Float64[]
  for k in 1:outer
    cf=mkcap(cfg,x); m=train_model(cf,Int64(100)); lb=SDDP.calculate_bound(m)
    d=simulate_cap_duals(m,cf,Int64(200))["duals"]; grad=d.+v; gn=norm(grad); g=grad./(gn+1e-12)
    xn=clamp.(x .- (a0/sqrt(k)).*g, 1.0, xmax); st=norm(xn.-x)
    push!(its,k); push!(objs,lb+sum(vv.*x)); push!(gns,gn); push!(steps,st); x=xn
    println("F5 iter $k gn=$(round(gn,digits=2))"); flush(stdout)
  end
  Dict("iter"=>its,"obj"=>objs,"gradnorm"=>gns,"step"=>steps)
end
d=traj(loadcfg("capreg_l20u95"; lambda=20.0, rho=0.4))
open(joinpath(FD,"capopt_traj.json"),"w") do io; JSON3.write(io,d); end
println("DONE F5")
