# Regime-dependence grid: how the VALUE of stochastic optimization varies with
# capacity scarcity (utilisation) and supply-demand correlation, on a FIXED
# 6x6, tau=12 topology.  Paired CRN: myopic, SDDP, and perfect-foresight are all
# evaluated on the SAME realised demand per trial (from simulate_model's noise).
#   gain_of_recourse = (myopic - SDDP)/myopic     (value vs greedy)
#   regret           = (SDDP - clairvoyant)/clair (loss vs perfect info)
using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics
include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))
const INST=expanduser("~/drayage/MSTP/rerun/inst")
toi(x)= x isa Number ? Int64[Int64(x)] : Vector{Int64}(Int64.(collect(x)))
tof(x)= x isa Number ? Float64[Float64(x)] : Vector{Float64}(Float64.(collect(x)))
function loadcfg(name; lambda, rho, capscale=1.0, n_scen=20)
  ic=JSON3.read(read(joinpath(INST,"$name.json"),String)); nO=Int64(ic.nOrigins); nD=Int64(ic.nDestinations)
  cap=Int64.(max.(1, round.(toi(ic.carrier_capacity).*capscale)))
  HyperParams(tau=Int64(ic.tau),nOrigins=nO,nDestinations=nD,nCarriers=Int64(ic.nCarriers),nSpotCarriers=Int64(ic.nSpotCarriers),
    Bids=[toi(b) for b in ic.Bids],Winners=[toi(w) for w in ic.Winners],
    entry_stock_0=toi(ic.entry_stock_0),exit_stock_0=toi(ic.exit_stock_0),exit_short_0=toi(ic.exit_short_0),
    entry_capacity=toi(ic.entry_capacity),exit_capacity=toi(ic.exit_capacity),
    entry_store_coef=tof(ic.entry_store_coef),exit_store_coef=tof(ic.exit_store_coef),exit_short_coef=tof(ic.exit_short_coef),
    transport_coef=tof(ic.transport_coef),spot_coef=tof(ic.spot_coef),carrier_capacity=cap,
    lambda=fill(Float64(lambda),nO+nD),corrmat=gen_cov_mat(2,nO,rho),n_scenarios=n_scen)
end
spot_kdx(cfg,t)=vec(((1:cfg.nSpotCarriers)' .-1).*cfg.nOrigins.*cfg.nDestinations.*cfg.tau .+ (cfg.SpotIdx[1].-1).*cfg.tau .+ t)
function clair(cfg,Q,D)
  T=cfg.tau;nO=cfg.nOrigins;nD=cfg.nDestinations;nK=cfg.nCarriers+cfg.nSpotCarriers
  m=Model(HiGHS.Optimizer);JuMP.set_silent(m)
  @variable(m,mv[1:(cfg.nLanes+cfg.nSpotLanes),1:T]>=0);@variable(m,ent[1:nO,1:T+1]>=0);@variable(m,ep[1:nD,1:T+1]>=0);@variable(m,em[1:nD,1:T+1]>=0)
  for i in 1:nO;JuMP.fix(ent[i,1],cfg.entry_stock_0[i];force=true);end
  for j in 1:nD;JuMP.fix(ep[j,1],cfg.exit_stock_0[j];force=true);JuMP.fix(em[j,1],cfg.exit_short_0[j];force=true);end
  @constraint(m,[k=1:nK,t=1:T],sum(mv[cfg.CarrierIdx[k],t])<=cfg.carrier_capacity[(k-1)*T+t])
  for t in 1:T
    for i in 1:nO
      @constraint(m,sum(mv[cfg.from_SG[i],t])+sum(mv[cfg.from_SP[i],t])<=ent[i,t]+Q[i,t])
      @constraint(m,ent[i,t+1]==ent[i,t]+Q[i,t]-sum(mv[cfg.from_SG[i],t])-sum(mv[cfg.from_SP[i],t]))
    end
    for j in 1:nD
      @constraint(m,ep[j,t]+sum(mv[cfg.to_SG[j],t])+sum(mv[cfg.to_SP[j],t])<=cfg.exit_capacity[j])
      @constraint(m,ep[j,t+1]-em[j,t+1]==ep[j,t]-em[j,t]+sum(mv[cfg.to_SG[j],t])+sum(mv[cfg.to_SP[j],t])-D[j,t])
    end
  end
  o=AffExpr(0.0)
  for t in 1:T
    for i in 1:nO;add_to_expression!(o,cfg.entry_store_coef[i],ent[i,t]);end
    for j in 1:nD;add_to_expression!(o,cfg.exit_store_coef[j],ep[j,t]);add_to_expression!(o,cfg.exit_short_coef[j],em[j,t]);end
    add_to_expression!(o,sum(cfg.transport_coef[1:cfg.nLanes].*mv[1:cfg.nLanes,t]))
    add_to_expression!(o,sum(cfg.spot_coef[spot_kdx(cfg,t)].*mv[(cfg.nLanes+1):(cfg.nLanes+cfg.nSpotLanes),t]))
  end
  @objective(m,Min,o);JuMP.optimize!(m);JuMP.objective_value(m)
end
# Greedy single-stage LP with one-step holding/shortage look-ahead (matches 20_vss.R).
function myopic_stage(cfg,ent,ep,em,Q,D,t)
  nO=cfg.nOrigins;nD=cfg.nDestinations;nK=cfg.nCarriers+cfg.nSpotCarriers
  m=Model(HiGHS.Optimizer);JuMP.set_silent(m)
  @variable(m,mv[1:(cfg.nLanes+cfg.nSpotLanes)]>=0);@variable(m,eo[1:nO]>=0);@variable(m,xp[1:nD]>=0);@variable(m,xm[1:nD]>=0)
  @constraint(m,[k=1:nK],sum(mv[cfg.CarrierIdx[k]])<=cfg.carrier_capacity[(k-1)*cfg.tau+t])
  for i in 1:nO;@constraint(m,eo[i]+sum(mv[cfg.from_SG[i]])+sum(mv[cfg.from_SP[i]])==ent[i]+Q[i]);end
  for j in 1:nD
    @constraint(m,sum(mv[cfg.to_SG[j]])+sum(mv[cfg.to_SP[j]])<=cfg.exit_capacity[j]-ep[j])
    @constraint(m,xp[j]-xm[j]-sum(mv[cfg.to_SG[j]])-sum(mv[cfg.to_SP[j]])==ep[j]-em[j]-D[j])
  end
  o=AffExpr(0.0)
  add_to_expression!(o,sum(cfg.transport_coef[1:cfg.nLanes].*mv[1:cfg.nLanes]))
  add_to_expression!(o,sum(cfg.spot_coef[spot_kdx(cfg,t)].*mv[(cfg.nLanes+1):(cfg.nLanes+cfg.nSpotLanes)]))
  for i in 1:nO;add_to_expression!(o,cfg.entry_store_coef[i],eo[i]);end
  for j in 1:nD;add_to_expression!(o,cfg.exit_store_coef[j],xp[j]);add_to_expression!(o,cfg.exit_short_coef[j],xm[j]);end
  @objective(m,Min,o);JuMP.optimize!(m)
  sc=sum(cfg.entry_store_coef.*ent)+sum(cfg.exit_store_coef.*ep)+sum(cfg.exit_short_coef.*em)
  (sc+JuMP.objective_value(m), JuMP.value.(eo), JuMP.value.(xp), JuMP.value.(xm))
end
function myopic_rollout(cfg,Q,D)
  ent=zeros(cfg.nOrigins);ep=zeros(cfg.nDestinations);em=zeros(cfg.nDestinations);tot=0.0
  for t in 1:cfg.tau; c,ent,ep,em=myopic_stage(cfg,ent,ep,em,Q[:,t],D[:,t],t); tot+=c; end
  tot
end
QD(n,tr,T,nO,nD)=begin Q=zeros(nO,T);D=zeros(nD,T);for t in 1:T,i in 1:nO;Q[i,t]=n[(tr-1)*T+t,i];end;for t in 1:T,j in 1:nD;D[j,t]=n[(tr-1)*T+t,nO+j];end;(Q,D) end
const BASE="capreg_l20u80"; const LAM=20.0; const NTR=300; const ITER=500
println("REGIME GRID on $BASE (6x6, tau=12), SDDP $ITER iters, $NTR CRN trials")
for capscale in (1.33,1.0,0.84)          # nominal util ~ 0.80/scale = 0.60, 0.80, 0.95
  util=round(0.80/capscale,digits=2)
  for rho in (0.0,0.4)
    cfg=loadcfg(BASE; lambda=LAM, rho=rho, capscale=capscale)
    mdl=train_model(cfg,Int64(ITER)); s=simulate_model(mdl,cfg,Int64(NTR)); obj=s["obj"]; noise=s["noise"]
    reg=Float64[]; gof=Float64[]
    for tr in 1:length(obj)
      Q,D=QD(noise,tr,cfg.tau,cfg.nOrigins,cfg.nDestinations)
      lp=clair(cfg,Q,D); my=myopic_rollout(cfg,Q,D)
      push!(reg, abs(lp)<1e-9 ? 0.0 : 100*(obj[tr]-lp)/abs(lp))
      push!(gof, abs(my)<1e-9 ? 0.0 : 100*(my-obj[tr])/abs(my))
    end
    println(">>> util=$util rho=$rho | gain_of_recourse=$(round(mean(gof),digits=1))% regret=$(round(mean(reg),digits=2))% | meanSDDP=$(round(mean(obj),digits=0)) meanMyopic=$(round(mean(myopic_rollout(cfg,QD(noise,1,cfg.tau,cfg.nOrigins,cfg.nDestinations)...)),digits=0))")
    flush(stdout)
  end
end
println("DONE")
