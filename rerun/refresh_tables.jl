# Paper-table refresh on the FIXED generator, pure Julia (no JuliaCall).
# Re-runs Tables 3 (scalability), 4 (horizons), 5 (sensitivity), 7 (capacity opt),
# seed-averaged, certifying every bound with validate_bound.
#
#   MODE=smoke julia rerun/refresh_tables.jl   # tiny, end-to-end sanity (~mins)
#   MODE=full  julia rerun/refresh_tables.jl   # overnight
#
# Reads rerun/inst/*.json (from gen_instances.R). Writes results to stdout as
# >>> lines and a final JSON-ish summary.
using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics

include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))
const SMOKE = get(ENV, "MODE", "full") == "smoke"
const INST  = expanduser("~/drayage/MSTP/rerun/inst")
toi(x)= x isa Number ? Int64[Int64(x)] : Vector{Int64}(Int64.(collect(x)))
tof(x)= x isa Number ? Float64[Float64(x)] : Vector{Float64}(Float64.(collect(x)))

# Load instance JSON → HyperParams, with optional sensitivity overrides.
function loadcfg(name; lambda=nothing, rho=0.0, m_mult=1.0, n_scen=20)
  ic=JSON3.read(read(joinpath(INST,"$name.json"),String)); nO=Int64(ic.nOrigins); nD=Int64(ic.nDestinations)
  lam = lambda===nothing ? 700.0 : Float64(lambda)
  HyperParams(tau=Int64(ic.tau),nOrigins=nO,nDestinations=nD,nCarriers=Int64(ic.nCarriers),nSpotCarriers=Int64(ic.nSpotCarriers),
    Bids=[toi(b) for b in ic.Bids],Winners=[toi(w) for w in ic.Winners],
    entry_stock_0=toi(ic.entry_stock_0),exit_stock_0=toi(ic.exit_stock_0),exit_short_0=toi(ic.exit_short_0),
    entry_capacity=toi(ic.entry_capacity),exit_capacity=toi(ic.exit_capacity),
    entry_store_coef=tof(ic.entry_store_coef),exit_store_coef=tof(ic.exit_store_coef),exit_short_coef=tof(ic.exit_short_coef),
    transport_coef=tof(ic.transport_coef),spot_coef=tof(ic.spot_coef).*m_mult,carrier_capacity=toi(ic.carrier_capacity),
    lambda=fill(lam,nO+nD),corrmat=gen_cov_mat(2,nO,rho),n_scenarios=n_scen)
end

# Deterministic multi-period LP at mean demand; returns (obj, capacity duals aligned to carrier_capacity).
function det_lp(cfg, x)
  T=cfg.tau; nO=cfg.nOrigins; nD=cfg.nDestinations; nK=cfg.nCarriers+cfg.nSpotCarriers; lam=cfg.lambda
  m=Model(HiGHS.Optimizer); JuMP.set_silent(m)
  @variable(m, mv[1:(cfg.nLanes+cfg.nSpotLanes),1:T]>=0)
  @variable(m, ent[1:nO,1:T+1]>=0); @variable(m, exp_[1:nD,1:T+1]>=0); @variable(m, exm[1:nD,1:T+1]>=0)
  for i in 1:nO; JuMP.fix(ent[i,1],cfg.entry_stock_0[i];force=true); end
  for j in 1:nD; JuMP.fix(exp_[j,1],cfg.exit_stock_0[j];force=true); JuMP.fix(exm[j,1],cfg.exit_short_0[j];force=true); end
  @constraint(m, cap[k=1:nK,t=1:T], sum(mv[cfg.CarrierIdx[k],t]) <= x[(k-1)*T+t])
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
  duals=zeros(nK*T); for k in 1:nK,t in 1:T; duals[(k-1)*T+t]=JuMP.dual(cap[k,t]); end
  (obj=JuMP.objective_value(m), duals=duals)
end

# Projected-gradient capacity optimisation. grad_source ∈ (:sddp,:lp).
function capopt(cfg, v; source, outer, cold, nsamp, a0=20.0)
  x=Float64.(cfg.carrier_capacity); xmax=copy(x); flo=1.0
  for k in 1:outer
    if source==:sddp
      cf=mkcap(cfg,x); m=train_model(cf,Int64(cold)); d=simulate_cap_duals(m,cf,Int64(nsamp))["duals"]
      grad=d.+v
    else
      _,d=det_lp(cfg,x); grad=d.+v
    end
    g=grad./(norm(grad)+1e-12); x=clamp.(x .- (a0/sqrt(k)).*g, flo, xmax)
  end
  x
end
mkcap(cfg,x)=HyperParams(tau=cfg.tau,nOrigins=cfg.nOrigins,nDestinations=cfg.nDestinations,nCarriers=cfg.nCarriers,nSpotCarriers=cfg.nSpotCarriers,
  Bids=cfg.Bids,Winners=cfg.Winners,entry_stock_0=cfg.entry_stock_0,exit_stock_0=cfg.exit_stock_0,exit_short_0=cfg.exit_short_0,
  entry_capacity=cfg.entry_capacity,exit_capacity=cfg.exit_capacity,entry_store_coef=cfg.entry_store_coef,exit_store_coef=cfg.exit_store_coef,
  exit_short_coef=cfg.exit_short_coef,transport_coef=cfg.transport_coef,spot_coef=cfg.spot_coef,
  carrier_capacity=Int64.(round.(x)),lambda=cfg.lambda,corrmat=cfg.corrmat,n_scenarios=cfg.n_scenarios)

# Train+certify one config; returns (lb,ub,se,gap,valid).
function cert(cfg, iters; trials=200)
  m=train_model(cfg,Int64(iters)); r=validate_bound(m,cfg;trials=Int64(trials),z=3.0)
  (lb=r["bound"],ub=r["sim_mean"],se=r["sim_se"],gap=100*(r["sim_mean"]-r["bound"])/abs(r["bound"]),valid=r["valid"])
end
avg(f,seeds)=(rs=[f(s) for s in seeds]; (lb=mean(getfield.(rs,:lb)),ub=mean(getfield.(rs,:ub)),
  gap=mean(getfield.(rs,:gap)),valid=all(getfield.(rs,:valid))))

# Paper protocol: Tables 3 & 4 use single seed 42; only Table 5 seed-averages.
const SEEDS_AVG = SMOKE ? [42] : [42,43,44]
println("MODE=", SMOKE ? "smoke" : "full")

# ---- Table 4: horizons (6x6x20), single seed 42 ----
println("\n=== TABLE 4 horizons ===")
for tau in (SMOKE ? (12,) : (12,26,52))
  it = SMOKE ? 20 : 300
  r=cert(loadcfg("hor_$(tau)_42"), it)
  println(">>> T4 tau=$tau LB=$(round(r.lb)) UB=$(round(r.ub)) gap%=$(round(r.gap,digits=2)) valid=$(r.valid)"); flush(stdout)
end

# ---- Table 5: sensitivity (6x6x20 tau=12) ----
println("\n=== TABLE 5 sensitivity ===")
for lam in (SMOKE ? (700,) : (200,700,2000))
  it = SMOKE ? 20 : 300
  r=avg(s->cert(loadcfg("sens_$(lam)_$(s)"; lambda=lam), it), SEEDS_AVG)
  println(">>> T5 lambda=$lam UB=$(round(r.ub)) valid=$(r.valid)"); flush(stdout)
end
if !SMOKE
  for rho in (0.0,0.2,0.4)
    r=avg(s->cert(loadcfg("sens_700_$(s)"; lambda=700, rho=rho), 300), SEEDS_AVG)
    println(">>> T5 rho=$rho UB=$(round(r.ub)) valid=$(r.valid)"); flush(stdout)
  end
  for mm in (1.0,2.0,4.0)
    r=avg(s->cert(loadcfg("sens_700_$(s)"; lambda=700, m_mult=mm), 300), SEEDS_AVG)
    println(">>> T5 m=$mm UB=$(round(r.ub)) valid=$(r.valid)"); flush(stdout)
  end
end

# ---- Table 3: scalability (single seed 42, paper iteration budgets) ----
println("\n=== TABLE 3 scalability ===")
let
  r=cert(loadcfg("scal_6_42"), SMOKE ? 20 : 500); println(">>> T3 6x6x20 LB=$(round(r.lb)) UB=$(round(r.ub)) gap%=$(round(r.gap,digits=2)) valid=$(r.valid)"); flush(stdout)
  if !SMOKE
    r=cert(loadcfg("scal_20_42"), 50); println(">>> T3 20x20x100 LB=$(round(r.lb)) UB=$(round(r.ub)) gap%=$(round(r.gap,digits=2)) valid=$(r.valid)"); flush(stdout)
    r=cert(loadcfg("scal_40_42"), 20); println(">>> T3 40x40x100 LB=$(round(r.lb)) UB=$(round(r.ub)) gap%=$(round(r.gap,digits=2)) valid=$(r.valid)"); flush(stdout)
  end
end

# ---- Table 7: capacity optimisation (6x6x20 lambda=50) ----
println("\n=== TABLE 7 capacity opt ===")
let cfg0=loadcfg("capopt_42"; lambda=50.0)
  v = 0.2*mean(cfg0.transport_coef)
  vvec = fill(v, length(cfg0.carrier_capacity))
  totcost(x)=(c=mkcap(cfg0,x); r=cert(c, SMOKE ? 50 : 100; trials= SMOKE ? 100 : 300); r.ub + sum(vvec.*x))
  x0=Float64.(cfg0.carrier_capacity)
  oc = SMOKE ? 3 : 15
  xlp  = capopt(cfg0, vvec; source=:lp,   outer=oc, cold=(SMOKE ? 50 : 100), nsamp=(SMOKE ? 50 : 200))
  xsd  = capopt(cfg0, vvec; source=:sddp, outer=oc, cold=(SMOKE ? 50 : 100), nsamp=(SMOKE ? 50 : 200))
  t0=totcost(x0); tl=totcost(xlp); ts=totcost(xsd)
  println(">>> T7 total x0=$(round(t0)) xLP=$(round(tl)) [$(round(100*(tl-t0)/t0,digits=1))%] xSDDP=$(round(ts)) [$(round(100*(ts-t0)/t0,digits=1))%] adv=$(round(100*(tl-ts)/t0,digits=1))pp")
  flush(stdout)
end
println("\nDONE")
