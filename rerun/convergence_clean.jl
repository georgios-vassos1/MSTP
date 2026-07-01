# Clean curvature sweep: which levers make SDDP's LOWER BOUND take many iters to
# converge, WITHOUT ill-conditioning the LP? Report %-of-final LB reached @10/@50.
using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics
include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))
const INST=expanduser("~/drayage/MSTP/rerun/inst")
toi(x)= x isa Number ? Int64[Int64(x)] : Vector{Int64}(Int64.(collect(x)))
tof(x)= x isa Number ? Float64[Float64(x)] : Vector{Float64}(Float64.(collect(x)))
function mkcfg(name; lambda, rho, spotcost=1.0, stratcap=1.0, spotcap=1.0, n_scen=20)
  ic=JSON3.read(read(joinpath(INST,"$name.json"),String)); nO=Int64(ic.nOrigins); nD=Int64(ic.nDestinations)
  nC=Int64(ic.nCarriers); nSC=Int64(ic.nSpotCarriers); T=Int64(ic.tau); cap=toi(ic.carrier_capacity); newcap=copy(cap)
  for k in 1:nC, t in 1:T; newcap[(k-1)*T+t]=max(1,round(Int,cap[(k-1)*T+t]*stratcap)); end
  for k in (nC+1):(nC+nSC), t in 1:T; newcap[(k-1)*T+t]=max(1,round(Int,cap[(k-1)*T+t]*spotcap)); end
  HyperParams(tau=T,nOrigins=nO,nDestinations=nD,nCarriers=nC,nSpotCarriers=nSC,
    Bids=[toi(b) for b in ic.Bids],Winners=[toi(w) for w in ic.Winners],
    entry_stock_0=toi(ic.entry_stock_0),exit_stock_0=toi(ic.exit_stock_0),exit_short_0=toi(ic.exit_short_0),
    entry_capacity=toi(ic.entry_capacity),exit_capacity=toi(ic.exit_capacity),
    entry_store_coef=tof(ic.entry_store_coef),exit_store_coef=tof(ic.exit_store_coef),exit_short_coef=tof(ic.exit_short_coef),
    transport_coef=tof(ic.transport_coef),spot_coef=tof(ic.spot_coef).*spotcost,carrier_capacity=newcap,
    lambda=fill(Float64(lambda),nO+nD),corrmat=gen_cov_mat(2,nO,rho),n_scenarios=n_scen)
end
function lb_traj(cfg, milestones)
  model=_make_model(cfg); ckpt=tempname()*".cuts.json"; completed=0; has_ckpt=false; out=Tuple{Int,Float64}[]
  for M in milestones
    while completed < M
      todo=min(25, M-completed); ok=false; tries=0
      while !ok && tries<6
        try
          SDDP.train(model; iteration_limit=todo, cut_type=SDDP.SINGLE_CUT, parallel_scheme=SDDP.Serial(),
                     cut_deletion_minimum=10, add_to_existing_cuts=(completed>0), print_level=0)
          completed+=todo; SDDP.write_cuts_to_file(model,ckpt); has_ckpt=true; ok=true
        catch e; tries+=1; has_ckpt && (model=_make_model(cfg); SDDP.read_cuts_from_file(model,ckpt)); end
      end
      ok || return out
    end
    lb=SDDP.calculate_bound(model); push!(out,(M,lb)); println("    .. $(cfg.tau)-stage @$M LB=$(round(lb,digits=0))"); flush(stdout)
  end
  out
end
const MS=[10,50,100,250]
function report(tag,cfg)
  local tr; try; tr=lb_traj(cfg,MS); catch e; println(">>> $tag CRASH:$(sprint(showerror,e))"); return; end
  isempty(tr) && (println(">>> $tag no-bound"); return)
  b=Dict(tr...); f=tr[end][2]; at(i)=get(b,i,f); pct(i)=abs(f)<1e-9 ? 100.0 : round(100*at(i)/f,digits=1)
  println(">>> $tag | %final @10=$(pct(10))% @50=$(pct(50))% @100=$(pct(100))% (final@$(tr[end][1])=$(round(f,digits=0)))"); flush(stdout)
end
println("CLEAN curvature sweep (LB %-of-final; lower %@10 = slower convergence)")
report("A base tau12 loose-spot ", mkcfg("capreg_l20u95";lambda=20.0,rho=0.0))
# (B,C limited-spot removed: HiGHS numerical fragility under stress)
report("D tau26 (hor_26_42)     ", mkcfg("hor_26_42";lambda=700.0,rho=0.0))
report("E tau52 (hor_52_42)     ", mkcfg("hor_52_42";lambda=700.0,rho=0.0))
println("DONE")
