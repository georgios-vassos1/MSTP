# Verify (1) the Table-7 capacity-opt reversal is real (trajectory-logged PGD,
# sensible x*), and (2) diagnose the m=2/4 sensitivity invalidity (per-seed margin).
using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics
include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))
const INST = expanduser("~/drayage/MSTP/rerun/inst")
toi(x)= x isa Number ? Int64[Int64(x)] : Vector{Int64}(Int64.(collect(x)))
tof(x)= x isa Number ? Float64[Float64(x)] : Vector{Float64}(Float64.(collect(x)))
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
mkcap(cfg,x)=HyperParams(tau=cfg.tau,nOrigins=cfg.nOrigins,nDestinations=cfg.nDestinations,nCarriers=cfg.nCarriers,nSpotCarriers=cfg.nSpotCarriers,
  Bids=cfg.Bids,Winners=cfg.Winners,entry_stock_0=cfg.entry_stock_0,exit_stock_0=cfg.exit_stock_0,exit_short_0=cfg.exit_short_0,
  entry_capacity=cfg.entry_capacity,exit_capacity=cfg.exit_capacity,entry_store_coef=cfg.entry_store_coef,exit_store_coef=cfg.exit_store_coef,
  exit_short_coef=cfg.exit_short_coef,transport_coef=cfg.transport_coef,spot_coef=cfg.spot_coef,
  carrier_capacity=Int64.(round.(x)),lambda=cfg.lambda,corrmat=cfg.corrmat,n_scenarios=cfg.n_scenarios)
function det_lp(cfg, x)
  T=cfg.tau; nO=cfg.nOrigins; nD=cfg.nDestinations; nK=cfg.nCarriers+cfg.nSpotCarriers; lam=cfg.lambda
  m=Model(HiGHS.Optimizer); JuMP.set_silent(m)
  @variable(m, mv[1:(cfg.nLanes+cfg.nSpotLanes),1:T]>=0); @variable(m, ent[1:nO,1:T+1]>=0)
  @variable(m, exp_[1:nD,1:T+1]>=0); @variable(m, exm[1:nD,1:T+1]>=0)
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
function cert(cfg, iters; trials=200)
  m=train_model(cfg,Int64(iters)); r=validate_bound(m,cfg;trials=Int64(trials),z=3.0)
  (lb=r["bound"],ub=r["sim_mean"],se=r["sim_se"],margin=r["margin_se"],valid=r["valid"])
end

# ---- (1) capacity-opt trajectory ----
println("=== T7 capacity-opt trajectory (lambda=50) ===")
cfg0=loadcfg("capopt_42"; lambda=50.0); v=0.2*mean(cfg0.transport_coef); vvec=fill(v,length(cfg0.carrier_capacity))
totcost(x)=(c=mkcap(cfg0,x); r=cert(c,100;trials=300); r.ub + sum(vvec.*x))
function capopt_traj(label; source, outer=20, a0=20.0)
  x=Float64.(cfg0.carrier_capacity); xmax=copy(x)
  println("  [$label] x0 mean=$(round(mean(x),digits=1))")
  for k in 1:outer
    if source==:sddp
      cf=mkcap(cfg0,x); m=train_model(cf,Int64(100)); d=simulate_cap_duals(m,cf,Int64(200))["duals"]; grad=d.+v
    else
      _,d=det_lp(cfg0,x); grad=d.+v
    end
    gn=norm(grad); g=grad./(gn+1e-12); x=clamp.(x .- (a0/sqrt(k)).*g, 1.0, xmax)
    (k%5==0||k==1) && println("  [$label] it=$k |grad|=$(round(gn,digits=2)) x_mean=$(round(mean(x),digits=1))")
    flush(stdout)
  end
  x
end
x0=Float64.(cfg0.carrier_capacity)
xlp=capopt_traj("LP"; source=:lp); xsd=capopt_traj("SDDP"; source=:sddp)
t0=totcost(x0); tl=totcost(xlp); ts=totcost(xsd)
println(">>> x0 mean=$(round(mean(x0),digits=1)) total=$(round(t0)) | LP mean=$(round(mean(xlp),digits=1)) total=$(round(tl)) [$(round(100*(tl-t0)/t0,digits=1))%] | SDDP mean=$(round(mean(xsd),digits=1)) total=$(round(ts)) [$(round(100*(ts-t0)/t0,digits=1))%] | adv=$(round(100*(tl-ts)/t0,digits=1))pp")
flush(stdout)

# ---- (2) m=2/4 per-seed margin ----
println("\n=== T5 m sensitivity per-seed margins ===")
for mm in (1.0,2.0,4.0), s in (42,43,44)
  r=cert(loadcfg("sens_700_$(s)"; lambda=700, m_mult=mm), 300)
  println(">>> m=$mm seed=$s UB=$(round(r.ub)) margin_se=$(round(r.margin,digits=2)) valid=$(r.valid)"); flush(stdout)
end
println("DONE")
