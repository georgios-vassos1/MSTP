using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics
include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))
const INST=expanduser("~/drayage/MSTP/rerun/inst")
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
    kdx=vec(((1:cfg.nSpotCarriers)' .-1).*nO.*nD.*T .+ (cfg.SpotIdx[1].-1).*T .+ t)
    add_to_expression!(o,sum(cfg.spot_coef[kdx].*mv[(cfg.nLanes+1):(cfg.nLanes+cfg.nSpotLanes),t]))
  end
  @objective(m,Min,o);JuMP.optimize!(m);JuMP.objective_value(m)
end
QD(n,tr,T,nO,nD)=begin Q=zeros(nO,T);D=zeros(nD,T);for t in 1:T,i in 1:nO;Q[i,t]=n[(tr-1)*T+t,i];end;for t in 1:T,j in 1:nD;D[j,t]=n[(tr-1)*T+t,nO+j];end;(Q,D) end
q95(v)=sort(v)[max(1,ceil(Int,0.95*length(v)))]
cfg=loadcfg("hor_12_42"; lambda=700, rho=0.0)
for it in (25,50,100,250,500,1000)
  m=train_model(cfg,Int64(it)); s=simulate_model(m,cfg,Int64(200)); obj=s["obj"]; noise=s["noise"]; r=Float64[]
  for tr in 1:length(obj); Q,D=QD(noise,tr,cfg.tau,cfg.nOrigins,cfg.nDestinations); lp=clair(cfg,Q,D); push!(r, abs(lp)<1e-9 ? 0.0 : 100*(obj[tr]-lp)/abs(lp)); end
  println(">>> iters=$it regret mean=$(round(mean(r),digits=3))% med=$(round(median(r),digits=3))% q95=$(round(q95(r),digits=3))%"); flush(stdout)
end
println("DONE")
