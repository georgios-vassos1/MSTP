# Reliable LP-proxy vs SDDP-dual capacity-opt comparison using COMMON RANDOM
# NUMBERS: both capacity vectors are evaluated on the SAME fixed OOB scenario
# set (variance of the difference is largely cancelled), with a higher-fidelity
# policy (300 iters) and 500 shared scenarios. Sweeps benign→rich uncertainty.
#
# Advantage_pp = (LP_total - SDDP_total)/x0_total*100  (positive ⇒ SDDP better).
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
function det_lp_duals(cfg,x)
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
pgd(cfg,v;source,outer=12,cold=100,a0=20.0)=begin
  x=Float64.(cfg.carrier_capacity); xmax=copy(x)
  for k in 1:outer
    grad = source==:lp ? det_lp_duals(cfg,x).+v : (cf=mkcap(cfg,x); simulate_cap_duals(train_model(cf,Int64(cold)),cf,Int64(200))["duals"].+v)
    x=clamp.(x .- (a0/sqrt(k)).*(grad./(norm(grad)+1e-12)), 1.0, xmax)
  end
  x
end
# Fixed OOB scenario set (common random numbers), matching simulate_model's construction.
make_oob(cfg,n)=[[(t, sample_scenarios(1,cfg.lambda,cfg.corrmat)...) for t in 1:cfg.tau] for _ in 1:n]
# Operational cost of the (300-iter) policy at capacity x, on the SHARED oob.
function op_on(cfg,x,oob)
  cf=mkcap(cfg,x); model=train_model(cf,Int64(300))
  sims=SDDP.simulate(model,length(oob),[:move]; sampling_scheme=SDDP.Historical(oob), parallel_scheme=SDDP.Serial())
  mean(sum(node[:stage_objective] for node in sim) for sim in sims)
end
total(cfg,vvec,x,oob)=op_on(cfg,x,oob)+sum(vvec.*x)

function run_regime(name,lam,rho; check=false)
  cfg=loadcfg(name; lambda=lam, rho=rho); v=0.2*mean(cfg.transport_coef); vvec=fill(v,length(cfg.carrier_capacity))
  x0=Float64.(cfg.carrier_capacity)
  xlp=pgd(cfg,v;source=:lp); xsd=pgd(cfg,v;source=:sddp)
  oob=make_oob(cfg,500)
  t0=total(cfg,vvec,x0,oob); tl=total(cfg,vvec,xlp,oob); ts=total(cfg,vvec,xsd,oob)
  adv=100*(tl-ts)/t0
  println(">>> $name rho=$rho | LP=$(round(100*(tl-t0)/t0,digits=1))% SDDP=$(round(100*(ts-t0)/t0,digits=1))% adv=$(round(adv,digits=1))pp (CRN) | xLP=$(round(mean(xlp),digits=1)) xSDDP=$(round(mean(xsd),digits=1))")
  flush(stdout)
  if check  # CRN stability: re-evaluate SAME x on a 2nd shared oob; adv should be close
    oob2=make_oob(cfg,500)
    t0b=total(cfg,vvec,x0,oob2); tlb=total(cfg,vvec,xlp,oob2); tsb=total(cfg,vvec,xsd,oob2)
    println("    [CRN check] adv(rep2)=$(round(100*(tlb-tsb)/t0b,digits=1))pp  (compare to $(round(adv,digits=1))pp)"); flush(stdout)
  end
end

println("=== CRN capacity-opt: cross-correlation impact at TIGHT capacity (l20u95) ===")
# Isolate ρ_cross at a scarce-capacity regime where variance bites (λ=20, util 0.95).
# (ρ=0.4 already measured: adv 32.5pp. Complete the sweep with ρ=0.0, 0.2.)
run_regime("capreg_l20u95",20.0,0.0; check=true)
run_regime("capreg_l20u95",20.0,0.2)
println("DONE")
