using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics, Random
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
    lambda=fill(Float64(lambda),nO+nD),corrmat=gen_cov_mat(2,max(nO,nD),rho)[1:(nO+nD),1:(nO+nD)],n_scenarios=n_scen)
end
mkcap(cfg,x)=HyperParams(tau=cfg.tau,nOrigins=cfg.nOrigins,nDestinations=cfg.nDestinations,nCarriers=cfg.nCarriers,nSpotCarriers=cfg.nSpotCarriers,
  Bids=cfg.Bids,Winners=cfg.Winners,entry_stock_0=cfg.entry_stock_0,exit_stock_0=cfg.exit_stock_0,exit_short_0=cfg.exit_short_0,
  entry_capacity=cfg.entry_capacity,exit_capacity=cfg.exit_capacity,entry_store_coef=cfg.entry_store_coef,exit_store_coef=cfg.exit_store_coef,
  exit_short_coef=cfg.exit_short_coef,transport_coef=cfg.transport_coef,spot_coef=cfg.spot_coef,
  carrier_capacity=Int64.(round.(x)),lambda=cfg.lambda,corrmat=cfg.corrmat,n_scenarios=cfg.n_scenarios)
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
QD(noise,tr,T,nO,nD)=begin Q=zeros(nO,T);D=zeros(nD,T);for t in 1:T,i in 1:nO;Q[i,t]=noise[(tr-1)*T+t,i];end;for t in 1:T,j in 1:nD;D[j,t]=noise[(tr-1)*T+t,nO+j];end;(Q,D) end
function regret_arr(cfg,iters,trials)
  m=train_model(cfg,Int64(iters));s=simulate_model(m,cfg,Int64(trials));obj=s["obj"];noise=s["noise"];r=Float64[]
  for tr in 1:length(obj);Q,D=QD(noise,tr,cfg.tau,cfg.nOrigins,cfg.nDestinations);lp=clair(cfg,Q,D);push!(r,abs(lp)<1e-9 ? 0.0 : 100*(obj[tr]-lp)/abs(lp));end
  r
end
# F3: regret distributions per topology
d3=Dict{String,Vector{Float64}}(); for t in ("t2x1","t1x2","t2x2"); d3[t]=regret_arr(loadcfg(t;lambda=5.5,rho=0.4),1000,500); println("F3 $t done"); flush(stdout); end
open(joinpath(FD,"regret_topo.json"),"w") do io; JSON3.write(io,d3); end
# F4: regret vs iterations (2x2)
d4=Dict{String,Vector{Float64}}(); for it in (100,250,500,1000,1500); d4[string(it)]=regret_arr(loadcfg("t2x2";lambda=5.5,rho=0.4),it,500); println("F4 $it done"); flush(stdout); end
open(joinpath(FD,"regret_iters.json"),"w") do io; JSON3.write(io,d4); end
# F5: capacity-opt trajectory (scarce regime), SDDP-dual path
cfg=loadcfg("capreg_l20u95";lambda=20.0,rho=0.4); v=0.2*mean(cfg.transport_coef); vv=fill(v,length(cfg.carrier_capacity))
x=Float64.(cfg.carrier_capacity);xmax=copy(x);its=Int[];objs=Float64[];gns=Float64[];steps=Float64[]
for k in 1:15
  cf=mkcap(cfg,x);m=train_model(cf,Int64(100));lb=SDDP.calculate_bound(m)
  d=simulate_cap_duals(m,cf,Int64(200))["duals"];grad=d.+v;gn=norm(grad);g=grad./(gn+1e-12)
  xn=clamp.(x .- (20.0/sqrt(k)).*g,1.0,xmax);step=norm(xn.-x)
  push!(its,k);push!(objs,lb+sum(vv.*x));push!(gns,gn);push!(steps,step);x=xn
  println("F5 iter $k gn=$(round(gn,digits=2))");flush(stdout)
end
open(joinpath(FD,"capopt_traj.json"),"w") do io; JSON3.write(io,Dict("iter"=>its,"obj"=>objs,"gradnorm"=>gns,"step"=>steps)); end
println("DONE figdata")
