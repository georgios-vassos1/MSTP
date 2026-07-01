# Tables 2 & 6 (clairvoyant-LP regret, regret-based stopping) and VSS
# (gain-of-recourse vs myopic) on the FIXED generator, pure Julia.
# Reimplements MSTP/R/recourse.R (compute_regret + .stage_lp/compute_vss).
using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics
include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))
const INST = expanduser("~/drayage/MSTP/rerun/inst")
const SMOKE = get(ENV,"MODE","full")=="smoke"
const NIT = SMOKE ? 30 : 1000
const NTR = SMOKE ? 20 : 500
toi(x)= x isa Number ? Int64[Int64(x)] : Vector{Int64}(Int64.(collect(x)))
tof(x)= x isa Number ? Float64[Float64(x)] : Vector{Float64}(Float64.(collect(x)))
function loadcfg(name; lambda=5.5, rho=0.4, n_scen=20)
  ic=JSON3.read(read(joinpath(INST,"$name.json"),String)); nO=Int64(ic.nOrigins); nD=Int64(ic.nDestinations)
  HyperParams(tau=Int64(ic.tau),nOrigins=nO,nDestinations=nD,nCarriers=Int64(ic.nCarriers),nSpotCarriers=Int64(ic.nSpotCarriers),
    Bids=[toi(b) for b in ic.Bids],Winners=[toi(w) for w in ic.Winners],
    entry_stock_0=toi(ic.entry_stock_0),exit_stock_0=toi(ic.exit_stock_0),exit_short_0=toi(ic.exit_short_0),
    entry_capacity=toi(ic.entry_capacity),exit_capacity=toi(ic.exit_capacity),
    entry_store_coef=tof(ic.entry_store_coef),exit_store_coef=tof(ic.exit_store_coef),exit_short_coef=tof(ic.exit_short_coef),
    transport_coef=tof(ic.transport_coef),spot_coef=tof(ic.spot_coef),carrier_capacity=toi(ic.carrier_capacity),
    lambda=fill(Float64(lambda),nO+nD),corrmat=gen_cov_mat(2,max(nO,nD),rho)[1:(nO+nD),1:(nO+nD)],n_scenarios=n_scen)
end
# spot cost for lane l (contracted-lane index), spot carrier c, stage t
spot_at(cfg,l,c,t)= cfg.spot_coef[(c-1)*cfg.nOrigins*cfg.nDestinations*cfg.tau + (l-1)*cfg.tau + t]

# Clairvoyant perfect-foresight multi-period LP with realised demand Q(nO×T), D(nD×T).
function clairvoyant(cfg, Q, D)
  T=cfg.tau; nO=cfg.nOrigins; nD=cfg.nDestinations; nK=cfg.nCarriers+cfg.nSpotCarriers
  m=Model(HiGHS.Optimizer); JuMP.set_silent(m)
  @variable(m, mv[1:(cfg.nLanes+cfg.nSpotLanes),1:T]>=0); @variable(m, ent[1:nO,1:T+1]>=0)
  @variable(m, exp_[1:nD,1:T+1]>=0); @variable(m, exm[1:nD,1:T+1]>=0)
  for i in 1:nO; JuMP.fix(ent[i,1],cfg.entry_stock_0[i];force=true); end
  for j in 1:nD; JuMP.fix(exp_[j,1],cfg.exit_stock_0[j];force=true); JuMP.fix(exm[j,1],cfg.exit_short_0[j];force=true); end
  @constraint(m,[k=1:nK,t=1:T], sum(mv[cfg.CarrierIdx[k],t]) <= cfg.carrier_capacity[(k-1)*T+t])
  for t in 1:T
    for i in 1:nO
      @constraint(m, sum(mv[cfg.from_SG[i],t])+sum(mv[cfg.from_SP[i],t]) <= ent[i,t]+Q[i,t])
      @constraint(m, ent[i,t+1]==ent[i,t]+Q[i,t]-sum(mv[cfg.from_SG[i],t])-sum(mv[cfg.from_SP[i],t]))
    end
    for j in 1:nD
      @constraint(m, exp_[j,t]+sum(mv[cfg.to_SG[j],t])+sum(mv[cfg.to_SP[j],t]) <= cfg.exit_capacity[j])
      @constraint(m, exp_[j,t+1]-exm[j,t+1]==exp_[j,t]-exm[j,t]+sum(mv[cfg.to_SG[j],t])+sum(mv[cfg.to_SP[j],t])-D[j,t])
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
  @objective(m,Min,obj); JuMP.optimize!(m); JuMP.objective_value(m)
end

# Myopic greedy: single-stage LP per period (objective includes next-state holding/shortage), state carried forward.
function myopic(cfg, Q, D)
  T=cfg.tau; nO=cfg.nOrigins; nD=cfg.nDestinations; nK=cfg.nCarriers+cfg.nSpotCarriers
  ein=Float64.(cfg.entry_stock_0); xpin=Float64.(cfg.exit_stock_0); xmin=Float64.(cfg.exit_short_0); tot=0.0
  for t in 1:T
    m=Model(HiGHS.Optimizer); JuMP.set_silent(m)
    @variable(m, mv[1:(cfg.nLanes+cfg.nSpotLanes)]>=0)
    @variable(m, eo[1:nO]>=0); @variable(m, xpo[1:nD]>=0); @variable(m, xmo[1:nD]>=0)
    @constraint(m,[k=1:nK], sum(mv[cfg.CarrierIdx[k]]) <= cfg.carrier_capacity[(k-1)*T+t])
    for i in 1:nO
      @constraint(m, eo[i]==ein[i]+Q[i,t]-sum(mv[cfg.from_SG[i]])-sum(mv[cfg.from_SP[i]]))
      @constraint(m, sum(mv[cfg.from_SG[i]])+sum(mv[cfg.from_SP[i]]) <= ein[i]+Q[i,t])
    end
    for j in 1:nD
      @constraint(m, xpin[j]+sum(mv[cfg.to_SG[j]])+sum(mv[cfg.to_SP[j]]) <= cfg.exit_capacity[j])
      @constraint(m, xpo[j]-xmo[j]==xpin[j]-xmin[j]+sum(mv[cfg.to_SG[j]])+sum(mv[cfg.to_SP[j]])-D[j,t])
    end
    obj=AffExpr(0.0)
    add_to_expression!(obj, sum(cfg.transport_coef[1:cfg.nLanes].*mv[1:cfg.nLanes]))
    kdx=vec(((1:cfg.nSpotCarriers)' .-1).*nO.*nD.*T .+ (cfg.SpotIdx[1].-1).*T .+ t)
    add_to_expression!(obj, sum(cfg.spot_coef[kdx].*mv[(cfg.nLanes+1):(cfg.nLanes+cfg.nSpotLanes)]))
    for i in 1:nO; add_to_expression!(obj,cfg.entry_store_coef[i],eo[i]); end
    for j in 1:nD; add_to_expression!(obj,cfg.exit_store_coef[j],xpo[j]); add_to_expression!(obj,cfg.exit_short_coef[j],xmo[j]); end
    @objective(m,Min,obj); JuMP.optimize!(m)
    # reported stage cost = current-state holding + optimised (move+next-state) cost
    tot += sum(cfg.entry_store_coef.*ein)+sum(cfg.exit_store_coef.*xpin)+sum(cfg.exit_short_coef.*xmin)+JuMP.objective_value(m)
    ein=max.(0.0,JuMP.value.(eo)); xpin=max.(0.0,JuMP.value.(xpo)); xmin=max.(0.0,JuMP.value.(xmo))
  end
  tot += sum(cfg.entry_store_coef.*ein)+sum(cfg.exit_store_coef.*xpin)+sum(cfg.exit_short_coef.*xmin)
  tot
end

# Extract per-trial realised (Q,D) matrices from simulate_model noise: (trials*T)×(nO+nD)
function trial_QD(noise, tr, T, nO, nD)
  Q=zeros(nO,T); Dm=zeros(nD,T)
  for t in 1:T
    row=(tr-1)*T+t
    for i in 1:nO; Q[i,t]=noise[row,i]; end
    for j in 1:nD; Dm[j,t]=noise[row,nO+j]; end
  end
  Q,Dm
end
function regret_of(cfg, iters, trials)
  m=train_model(cfg,Int64(iters)); s=simulate_model(m,cfg,Int64(trials))
  obj=s["obj"]; noise=s["noise"]; T=cfg.tau; nO=cfg.nOrigins; nD=cfg.nDestinations
  reg=Float64[]
  for tr in 1:length(obj)
    Q,Dm=trial_QD(noise,tr,T,nO,nD); lp=clairvoyant(cfg,Q,Dm)
    push!(reg, abs(lp)<1e-9 ? 0.0 : (obj[tr]-lp)/abs(lp))
  end
  reg, obj, noise
end

q95(v)=sort(v)[max(1,ceil(Int,0.95*length(v)))]
println("=== TABLE 2 regret across topologies (tau=4, lambda=5.5, rho=0.4, 1000 iters, 500 OOB) ===")
for topo in (SMOKE ? ("t2x2",) : ("t2x1","t1x2","t2x2"))
  reg,_,_ = regret_of(loadcfg(topo), NIT, NTR)
  rp = 100 .* reg
  println(">>> $topo mean=$(round(mean(rp),digits=3))% median=$(round(median(rp),digits=3))% q95=$(round(q95(rp),digits=3))%"); flush(stdout)
end

println("\n=== VSS gain-of-recourse (2x2, 500 OOB) ===")
let cfg=loadcfg("t2x2")
  m=train_model(cfg,Int64(NIT)); s=simulate_model(m,cfg,Int64(NTR)); obj=s["obj"]; noise=s["noise"]
  gains=Float64[]
  for tr in 1:length(obj)
    Q,Dm=trial_QD(noise,tr,cfg.tau,cfg.nOrigins,cfg.nDestinations); my=myopic(cfg,Q,Dm)
    push!(gains, (my-obj[tr])/my)
  end
  gp=100 .* gains
  println(">>> VSS 2x2 mean=$(round(mean(gp),digits=1))% median=$(round(median(gp),digits=1))% posfrac=$(round(mean(gains.>0),digits=3))"); flush(stdout)
end

println("\n=== TABLE 6 regret-based stopping (2x2, 500 OOB) ===")
for it in (SMOKE ? (100,) : (100,250,500,1000,1500))
  reg,_,_ = regret_of(loadcfg("t2x2"), it, NTR); rp=100 .* reg
  println(">>> iters=$it mean=$(round(mean(rp),digits=3))% sd=$(round(std(rp),digits=3))% q95=$(round(q95(rp),digits=3))%"); flush(stdout)
end
println("DONE")
