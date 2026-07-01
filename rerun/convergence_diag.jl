# WHY does SDDP converge in ~10 iters, and what makes it need more?
# Diagnostic: track the SDDP LOWER-BOUND trajectory (SDDP.train's own log) as we
# dial up value-function CURVATURE. Hypothesis: a cheap/unbounded spot backstop
# linearises the cost-to-go (few cuts needed); limited+expensive spot, steep
# shortage penalties, and tight capacity make it convex (many cuts needed).
using JSON3, SDDP, JuMP, HiGHS, LinearAlgebra, Statistics
include(expanduser("~/drayage/MSTP/inst/julia/mstp.jl"))
const INST=expanduser("~/drayage/MSTP/rerun/inst")
toi(x)= x isa Number ? Int64[Int64(x)] : Vector{Int64}(Int64.(collect(x)))
tof(x)= x isa Number ? Float64[Float64(x)] : Vector{Float64}(Float64.(collect(x)))
# Build cfg from instance with curvature multipliers.
function mkcfg(name; lambda, rho, spotcost=1.0, shortmult=1.0, stratcap=1.0, spotcap=1.0, n_scen=20)
  ic=JSON3.read(read(joinpath(INST,"$name.json"),String)); nO=Int64(ic.nOrigins); nD=Int64(ic.nDestinations)
  nC=Int64(ic.nCarriers); nSC=Int64(ic.nSpotCarriers); T=Int64(ic.tau)
  cap=toi(ic.carrier_capacity)                       # length (nC+nSC)*T
  newcap=copy(cap)
  for k in 1:nC, t in 1:T; newcap[(k-1)*T+t]=max(1,round(Int, cap[(k-1)*T+t]*stratcap)); end
  for k in (nC+1):(nC+nSC), t in 1:T; newcap[(k-1)*T+t]=max(1,round(Int, cap[(k-1)*T+t]*spotcap)); end
  HyperParams(tau=T,nOrigins=nO,nDestinations=nD,nCarriers=nC,nSpotCarriers=nSC,
    Bids=[toi(b) for b in ic.Bids],Winners=[toi(w) for w in ic.Winners],
    entry_stock_0=toi(ic.entry_stock_0),exit_stock_0=toi(ic.exit_stock_0),exit_short_0=toi(ic.exit_short_0),
    entry_capacity=toi(ic.entry_capacity),exit_capacity=toi(ic.exit_capacity),
    entry_store_coef=tof(ic.entry_store_coef),exit_store_coef=tof(ic.exit_store_coef),
    exit_short_coef=tof(ic.exit_short_coef).*shortmult,
    transport_coef=tof(ic.transport_coef),spot_coef=tof(ic.spot_coef).*spotcost,carrier_capacity=newcap,
    lambda=fill(Float64(lambda),nO+nD),corrmat=gen_cov_mat(2,nO,rho),n_scenarios=n_scen)
end
# Robust incremental LB trajectory: train in tiny batches with checkpoint+retry
# (mirrors train_model's crash safety), record calculate_bound at each milestone.
function lb_traj(cfg, milestones)
  model=_make_model(cfg); ckpt=tempname()*".cuts.json"; completed=0; has_ckpt=false
  out=Tuple{Int,Float64}[]
  for M in milestones
    while completed < M
      todo=min(25, M-completed); ok=false; tries=0
      while !ok && tries<6
        try
          SDDP.train(model; iteration_limit=todo, cut_type=SDDP.SINGLE_CUT, parallel_scheme=SDDP.Serial(),
                     cut_deletion_minimum=10, add_to_existing_cuts=(completed>0), print_level=0)
          completed+=todo; SDDP.write_cuts_to_file(model,ckpt); has_ckpt=true; ok=true
        catch e
          tries+=1
          if has_ckpt; model=_make_model(cfg); SDDP.read_cuts_from_file(model,ckpt); end
        end
      end
      ok || return out   # gave up after retries
    end
    lb=SDDP.calculate_bound(model); push!(out,(M,lb)); println("    .. @$M LB=$(round(lb,digits=0))"); flush(stdout)
  end
  out
end
const MS=[10,50,100,250]
function report(tag,cfg)
  local tr
  try
    tr=lb_traj(cfg,MS)
  catch e
    println(">>> $tag | CRASH: $(sprint(showerror,e))"); flush(stdout); return
  end
  isempty(tr) && (println(">>> $tag | no bound recorded"); return)
  b=Dict(tr...); final=tr[end][2]
  at(i)= get(b,i, tr[end][2]); pct(i)= abs(final)<1e-9 ? 100.0 : round(100*at(i)/final,digits=1)
  println(">>> $tag | LB@10=$(round(at(10),digits=0)) @50=$(round(at(50),digits=0)) @100=$(round(at(100),digits=0)) @$(tr[end][1])=$(round(final,digits=0)) | %final @10=$(pct(10))% @50=$(pct(50))% @100=$(pct(100))% @250=$(pct(250))%")
  flush(stdout)
end
println("SDDP LB-trajectory diagnostic on capreg_l20u95 (6x6 tau=12), milestones=$MS, rho=0")
println("(convergence = %of-final LB reached early; low %@10 => needs many iterations)")
report("V0 baseline           ", mkcfg("capreg_l20u95";lambda=20.0,rho=0.0))
report("V2 spot x5 short x20  ", mkcfg("capreg_l20u95";lambda=20.0,rho=0.0,spotcost=5.0,shortmult=20.0))
report("V4 +spotcap x0.3      ", mkcfg("capreg_l20u95";lambda=20.0,rho=0.0,spotcost=5.0,shortmult=20.0,stratcap=0.6,spotcap=0.3))
println("DONE")
