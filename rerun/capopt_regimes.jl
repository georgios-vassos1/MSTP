# Does the SDDP capacity-opt advantage over the LP proxy reappear under RICHER
# uncertainty? Sweep (lambda/CV, utilisation-tightness, cross-correlation).
# Advantage_pp = (LP_total - SDDP_total)/x0_total * 100  (positive ⇒ SDDP better).
using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics
include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))
const INST = expanduser("~/drayage/MSTP/rerun/inst")
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
function det_lp(cfg,x)
  T=cfg.tau; nO=cfg.nOrigins; nD=cfg.nDestinations; nK=cfg.nCarriers+cfg.nSpotCarriers; lam=cfg.lambda
  m=Model(HiGHS.Optimizer); JuMP.set_silent(m)
  @variable(m, mv[1:(cfg.nLanes+cfg.nSpotLanes),1:T]>=0); @variable(m, ent[1:nO,1:T+1]>=0)
  @variable(m, exp_[1:nD,1:T+1]>=0); @variable(m, exm[1:nD,1:T+1]>=0)
  for i in 1:nO; JuMP.fix(ent[i,1],cfg.entry_stock_0[i];force=true); end
  for j in 1:nD; JuMP.fix(exp_[j,1],cfg.exit_stock_0[j];force=true); JuMP.fix(exm[j,1],cfg.exit_short_0[j];force=true); end
  @constraint(m,cap[k=1:nK,t=1:T], sum(mv[cfg.CarrierIdx[k],t]) <= x[(k-1)*T+t])
  for t in 1:T
    for i in 1:nO
      @constraint(m, sum(mv[cfg.from_SG[i],t])+sum(mv[cfg.from_SP[i],t]) <= ent[i,t]+lam[i])
      @constraint(m, ent[i,t+1]==ent[i,t]+lam[i]-sum(mv[cfg.from_SG[i],t])-sum(mv[cfg.from_SP[i],t]))
    end
    for j in 1:nD
      @constraint(m, exp_[j,t]+sum(mv[cfg.to_SG[j],t])+sum(mv[cfg.to_SP[j],t]) <= cfg.exit_capacity[j])
      @constraint(m, exp_[j,t+1]-exm[j,t+1]==exp_[j,t]-exm[j,t]+sum(mv[cfg.to_SG[j],t])+sum(mv[cfg.to_SP[j],t])-lam[nO+j])
    end
  end
  obj=AffExpr(0.0)
  for t in 1:T
    for i in 1:nO; add_to_expression!(obj,cfg.entry_store_coef[i],ent[i,t]); end
    for j in 1:nD; add_to_expression!(obj,cfg.exit_store_coef[j],exp_[j,t]); add_to_expression!(obj,cfg.exit_short_coef[j],exm[j,t]); end
    add_to_expression!(obj, sum(cfg.transport_coef[1:cfg.nLanes].*mv[1:cfg.nLanes,t]))
    kdx=vec(((1:cfg.nSpotCarriers)' .-1).*nO.*nD.*T .+ (cfg.SpotIdx[1].-1).*T .+ t)
    add_to_expression!(obj, sum(cfg.spot_coef[kdx].*mv[(cfg.nLanes+1):(cfg.nLanes+cfg.nSpotLanes),t]))
  end
  @objective(m,Min,obj); JuMP.optimize!(m)
  d=zeros(nK*T); for k in 1:nK,t in 1:T; d[(k-1)*T+t]=JuMP.dual(cap[k,t]); end; d
end
function pgd(cfg,v,vvec;source,outer=15,a0=20.0)
  x=Float64.(cfg.carrier_capacity); xmax=copy(x)
  for k in 1:outer
    grad = source==:lp ? det_lp(cfg,x).+v : (cf=mkcap(cfg,x); simulate_cap_duals(train_model(cf,Int64(100)),cf,Int64(200))["duals"].+v)
    x=clamp.(x .- (a0/sqrt(k)).*(grad./(norm(grad)+1e-12)), 1.0, xmax)
  end
  x
end
totcost(cfg,vvec,x)=(c=mkcap(cfg,x); r=validate_bound(train_model(c,Int64(100)),c;trials=Int64(300),z=3.0); r["sim_mean"]+sum(vvec.*x))

for (name,lam) in (("capreg_l50u80",50.0),("capreg_l50u95",50.0),("capreg_l20u80",20.0),("capreg_l20u95",20.0)), rho in (0.0,0.4)
  cfg=loadcfg(name; lambda=lam, rho=rho)
  v=0.2*mean(cfg.transport_coef); vvec=fill(v,length(cfg.carrier_capacity)); x0=Float64.(cfg.carrier_capacity)
  xlp=pgd(cfg,v,vvec;source=:lp); xsd=pgd(cfg,v,vvec;source=:sddp)
  t0=totcost(cfg,vvec,x0); tl=totcost(cfg,vvec,xlp); ts=totcost(cfg,vvec,xsd)
  println(">>> $name rho=$rho | LP=$(round(100*(tl-t0)/t0,digits=1))% SDDP=$(round(100*(ts-t0)/t0,digits=1))% adv=$(round(100*(tl-ts)/t0,digits=1))pp | xLP=$(round(mean(xlp),digits=1)) xSDDP=$(round(mean(xsd),digits=1)) x0=$(round(mean(x0),digits=1))")
  flush(stdout)
end
println("DONE")
