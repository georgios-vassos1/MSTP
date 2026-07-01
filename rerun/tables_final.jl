# Final table numbers with ALL columns (adds t_train and sd, plus LB/UB for the
# regret and stopping tables) on the FIXED generator. Pure Julia.
#   julia rerun/tables_final.jl         # fast tables (2,4,5,6,VSS + 6x6 of T3)
#   BIG=1 julia rerun/tables_final.jl   # also 20x20x100 and 40x40x100 (slow)
using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics
include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))
const INST = expanduser("~/drayage/MSTP/rerun/inst")
const BIG  = get(ENV,"BIG","0") == "1"
toi(x)= x isa Number ? Int64[Int64(x)] : Vector{Int64}(Int64.(collect(x)))
tof(x)= x isa Number ? Float64[Float64(x)] : Vector{Float64}(Float64.(collect(x)))
function loadcfg(name; lambda, rho, m_mult=1.0, n_scen=20)
  ic=JSON3.read(read(joinpath(INST,"$name.json"),String)); nO=Int64(ic.nOrigins); nD=Int64(ic.nDestinations)
  HyperParams(tau=Int64(ic.tau),nOrigins=nO,nDestinations=nD,nCarriers=Int64(ic.nCarriers),nSpotCarriers=Int64(ic.nSpotCarriers),
    Bids=[toi(b) for b in ic.Bids],Winners=[toi(w) for w in ic.Winners],
    entry_stock_0=toi(ic.entry_stock_0),exit_stock_0=toi(ic.exit_stock_0),exit_short_0=toi(ic.exit_short_0),
    entry_capacity=toi(ic.entry_capacity),exit_capacity=toi(ic.exit_capacity),
    entry_store_coef=tof(ic.entry_store_coef),exit_store_coef=tof(ic.exit_store_coef),exit_short_coef=tof(ic.exit_short_coef),
    transport_coef=tof(ic.transport_coef),spot_coef=tof(ic.spot_coef).*m_mult,carrier_capacity=toi(ic.carrier_capacity),
    lambda=fill(Float64(lambda),nO+nD),corrmat=gen_cov_mat(2,max(nO,nD),rho)[1:(nO+nD),1:(nO+nD)],n_scenarios=n_scen)
end
# Train + certify with timing and dispersion.
function cert(cfg, iters; trials=200)
  t_train=@elapsed (m=train_model(cfg,Int64(iters)))
  t_sim  =@elapsed (r=validate_bound(m,cfg;trials=Int64(trials),z=3.0))
  n=r["trials"]; sd=r["sim_se"]*sqrt(n)
  (m=m, lb=r["bound"], ub=r["sim_mean"], sd=sd, gap=100*(r["sim_mean"]-r["bound"])/abs(r["bound"]),
   valid=r["valid"], t_train=t_train, t_sim=t_sim, iters=iters)
end

# ---- Clairvoyant LP + myopic (for regret / VSS), demand from realised noise ----
function clairvoyant(cfg,Q,D)
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
function myopic(cfg,Q,D)
  T=cfg.tau; nO=cfg.nOrigins; nD=cfg.nDestinations; nK=cfg.nCarriers+cfg.nSpotCarriers
  ein=Float64.(cfg.entry_stock_0); xpin=Float64.(cfg.exit_stock_0); xmin=Float64.(cfg.exit_short_0); tot=0.0
  for t in 1:T
    m=Model(HiGHS.Optimizer); JuMP.set_silent(m)
    @variable(m, mv[1:(cfg.nLanes+cfg.nSpotLanes)]>=0); @variable(m, eo[1:nO]>=0); @variable(m, xpo[1:nD]>=0); @variable(m, xmo[1:nD]>=0)
    @constraint(m,[k=1:nK], sum(mv[cfg.CarrierIdx[k]]) <= cfg.carrier_capacity[(k-1)*T+t])
    for i in 1:nO
      @constraint(m, eo[i]==ein[i]+Q[i,t]-sum(mv[cfg.from_SG[i]])-sum(mv[cfg.from_SP[i]]))
      @constraint(m, sum(mv[cfg.from_SG[i]])+sum(mv[cfg.from_SP[i]]) <= ein[i]+Q[i,t])
    end
    for j in 1:nD
      @constraint(m, xpin[j]+sum(mv[cfg.to_SG[j]])+sum(mv[cfg.to_SP[j]]) <= cfg.exit_capacity[j])
      @constraint(m, xpo[j]-xmo[j]==xpin[j]-xmin[j]+sum(mv[cfg.to_SG[j]])+sum(mv[cfg.to_SP[j]])-D[j,t])
    end
    obj=AffExpr(0.0); add_to_expression!(obj, sum(cfg.transport_coef[1:cfg.nLanes].*mv[1:cfg.nLanes]))
    kdx=vec(((1:cfg.nSpotCarriers)' .-1).*nO.*nD.*T .+ (cfg.SpotIdx[1].-1).*T .+ t)
    add_to_expression!(obj, sum(cfg.spot_coef[kdx].*mv[(cfg.nLanes+1):(cfg.nLanes+cfg.nSpotLanes)]))
    for i in 1:nO; add_to_expression!(obj,cfg.entry_store_coef[i],eo[i]); end
    for j in 1:nD; add_to_expression!(obj,cfg.exit_store_coef[j],xpo[j]); add_to_expression!(obj,cfg.exit_short_coef[j],xmo[j]); end
    @objective(m,Min,obj); JuMP.optimize!(m)
    tot += sum(cfg.entry_store_coef.*ein)+sum(cfg.exit_store_coef.*xpin)+sum(cfg.exit_short_coef.*xmin)+JuMP.objective_value(m)
    ein=max.(0.0,JuMP.value.(eo)); xpin=max.(0.0,JuMP.value.(xpo)); xmin=max.(0.0,JuMP.value.(xmo))
  end
  tot + sum(cfg.entry_store_coef.*ein)+sum(cfg.exit_store_coef.*xpin)+sum(cfg.exit_short_coef.*xmin)
end
function QD(noise,tr,T,nO,nD)
  Q=zeros(nO,T); Dm=zeros(nD,T)
  for t in 1:T, i in 1:nO; Q[i,t]=noise[(tr-1)*T+t,i]; end
  for t in 1:T, j in 1:nD; Dm[j,t]=noise[(tr-1)*T+t,nO+j]; end
  Q,Dm
end
q95(v)=sort(v)[max(1,ceil(Int,0.95*length(v)))]

# Regret table row: trains, then per-trial regret vs clairvoyant. Returns LB, UB, sd, t_train, regret stats.
function regret_row(cfg, iters, trials)
  c=cert(cfg,iters;trials=trials); s=simulate_model(c.m,cfg,Int64(trials)); obj=s["obj"]; noise=s["noise"]
  reg=Float64[]; for tr in 1:length(obj); Q,Dm=QD(noise,tr,cfg.tau,cfg.nOrigins,cfg.nDestinations); lp=clairvoyant(cfg,Q,Dm); push!(reg, abs(lp)<1e-9 ? 0.0 : 100*(obj[tr]-lp)/abs(lp)); end
  (lb=c.lb, ub=mean(obj), sd=std(obj), t_train=c.t_train, rmean=mean(reg), rmed=median(reg), rsd=std(reg), rq95=q95(reg))
end

println("=== TABLE 2 regret (tau=4, lambda=5.5, rho=0.4, 1000 iter, 500 OOB) ===")
for topo in ("t2x1","t1x2","t2x2")
  r=regret_row(loadcfg(topo; lambda=5.5, rho=0.4), 1000, 500)
  println(">>> T2 $topo t_train=$(round(r.t_train,digits=1))s LB=$(round(r.lb)) UB=$(round(r.ub))±$(round(r.sd)) regret: mean=$(round(r.rmean,digits=3))% med=$(round(r.rmed,digits=3))% q95=$(round(r.rq95,digits=3))%"); flush(stdout)
end

println("\n=== VSS (2x2, 500 OOB) ===")
let cfg=loadcfg("t2x2"; lambda=5.5, rho=0.4)
  c=cert(cfg,1000;trials=500); s=simulate_model(c.m,cfg,Int64(500)); obj=s["obj"]; noise=s["noise"]
  g=Float64[]; for tr in 1:length(obj); Q,Dm=QD(noise,tr,cfg.tau,cfg.nOrigins,cfg.nDestinations); my=myopic(cfg,Q,Dm); push!(g,100*(my-obj[tr])/my); end
  println(">>> VSS 2x2 mean=$(round(mean(g),digits=1))% med=$(round(median(g),digits=1))% posfrac=$(round(mean(g.>0),digits=3))"); flush(stdout)
end

println("\n=== TABLE 4 horizons (6x6x20, seed 42, 300 iter, 200 OOB) ===")
for tau in (12,26,52)
  r=cert(loadcfg("hor_$(tau)_42"; lambda=700, rho=0.0), 300)
  println(">>> T4 tau=$tau t_train=$(round(r.t_train,digits=1))s LB=$(round(r.lb)) UB=$(round(r.ub))±$(round(r.sd)) UBper=$(round(r.ub/tau)) gap%=$(round(r.gap,digits=2)) valid=$(r.valid)"); flush(stdout)
end

println("\n=== TABLE 5 sensitivity (6x6x20, tau=12, seed-avg 42-44, 300 iter, 200 OOB) ===")
avgUB(f)=(rs=[f(s) for s in (42,43,44)]; (ub=mean(getfield.(rs,:ub)), sd=mean(getfield.(rs,:sd)), valid=all(getfield.(rs,:valid))))
for lam in (200,700,2000); r=avgUB(s->cert(loadcfg("sens_$(lam)_$(s)"; lambda=lam, rho=0.0),300)); println(">>> T5 lambda=$lam UB=$(round(r.ub))±$(round(r.sd)) valid=$(r.valid)"); flush(stdout); end
for rho in (0.0,0.2,0.4); r=avgUB(s->cert(loadcfg("sens_700_$(s)"; lambda=700, rho=rho),300)); println(">>> T5 rho=$rho UB=$(round(r.ub))±$(round(r.sd)) valid=$(r.valid)"); flush(stdout); end
for mm in (1.0,2.0,4.0); r=avgUB(s->cert(loadcfg("sens_700_$(s)"; lambda=700, rho=0.0, m_mult=mm),300)); println(">>> T5 m=$mm UB=$(round(r.ub))±$(round(r.sd)) valid=$(r.valid)"); flush(stdout); end

println("\n=== TABLE 6 stopping (2x2, 500 OOB) ===")
for it in (100,250,500,1000,1500)
  r=regret_row(loadcfg("t2x2"; lambda=5.5, rho=0.4), it, 500)
  println(">>> T6 iters=$it LB=$(round(r.lb)) UB=$(round(r.ub))±$(round(r.sd)) gap%=$(round(100*(r.ub-r.lb)/abs(r.lb),digits=2)) regret_mean=$(round(r.rmean,digits=3))% sd=$(round(r.rsd,digits=3))% q95=$(round(r.rq95,digits=3))%"); flush(stdout)
end

println("\n=== TABLE 3 scalability (seed 42) ===")
let r=cert(loadcfg("scal_6_42"; lambda=700, rho=0.0), 500)
  println(">>> T3 6x6x20 iters=500 t_train=$(round(r.t_train,digits=1))s t_sim=$(round(r.t_sim,digits=1))s LB=$(round(r.lb)) UB=$(round(r.ub))±$(round(r.sd)) gap%=$(round(r.gap,digits=2))"); flush(stdout)
end
if BIG
  for (nm,it) in (("scal_20_42",50),("scal_40_42",20))
    r=cert(loadcfg(nm; lambda=700, rho=0.0), it)
    println(">>> T3 $nm iters=$it t_train=$(round(r.t_train,digits=1))s t_sim=$(round(r.t_sim,digits=1))s LB=$(round(r.lb)) UB=$(round(r.ub))±$(round(r.sd)) gap%=$(round(r.gap,digits=2))"); flush(stdout)
  end
end
println("DONE")
