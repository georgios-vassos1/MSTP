# FINAL clean convergence run: LB (SDDP.calculate_bound) and out-of-sample UB
# (mean policy cost on a FIXED CRN scenario set) vs iteration count, on the
# scarce 6x6 tau=12 instance.  batch=1 + checkpoint-per-iter absorbs the
# occasional infeasible-node crash so successful iterations always accumulate.
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
make_oob(cfg,n)=[[(t, sample_scenarios(1,cfg.lambda,cfg.corrmat)...) for t in 1:cfg.tau] for _ in 1:n]
function ub_on(model,oob)
  sims=SDDP.simulate(model,length(oob),[:move]; sampling_scheme=SDDP.Historical(oob), parallel_scheme=SDDP.Serial())
  c=[sum(node[:stage_objective] for node in sim) for sim in sims]; (mean(c), std(c)/sqrt(length(c)))
end
cfg=loadcfg("capreg_l20u95"; lambda=20.0, rho=0.0)
oob=make_oob(cfg,300)                          # fixed CRN evaluation set
MS=[2,5,10,15,20,30,40,50,65,80,100,125,150,200,250,300]
# ONE training trajectory: cuts accumulate monotonically; on any numerical crash
# rebuild the model and reload accumulated cuts (mirrors train_model's recovery,
# ALWAYS rebuilding on crash -- the fix vs the earlier inline loop). Record LB
# (calculate_bound) and CRN out-of-sample UB whenever we cross a milestone.
function run_trajectory(cfg, oob, MS)
  iters=Int[]; lbs=Float64[]; ubs=Float64[]; ses=Float64[]
  model=_make_model(cfg); ckpt=tempname()*".cuts.json"; completed=0; has_ckpt=false; retries=0; msidx=1
  maxM=maximum(MS)
  while completed < maxM && msidx <= length(MS)
    to_do=min(1, MS[msidx]-completed)
    try
      SDDP.train(model; iteration_limit=to_do, cut_type=SDDP.SINGLE_CUT, parallel_scheme=SDDP.Serial(),
                 cut_deletion_minimum=10, duality_handler=SDDP.ContinuousConicDuality(),
                 add_to_existing_cuts=(completed>0), print_level=0)
      completed+=to_do; SDDP.write_cuts_to_file(model,ckpt); has_ckpt=true; retries=0
    catch e
      retries+=1
      model=_make_model(cfg); has_ckpt && SDDP.read_cuts_from_file(model,ckpt)   # ALWAYS rebuild
      retries>8 && (completed+=to_do; retries=0)                                 # skip stuck iter
      continue
    end
    if completed >= MS[msidx]
      lb=SDDP.calculate_bound(model); ub,se=ub_on(model,oob)
      push!(iters,completed); push!(lbs,lb); push!(ubs,ub); push!(ses,se)
      println(">>> iter=$completed LB=$(round(lb,digits=0)) UB=$(round(ub,digits=0)) gap=$(round(100*(ub-lb)/lb,digits=1))%"); flush(stdout)
      msidx+=1
    end
  end
  iters,lbs,ubs,ses
end
iters,lbs,ubs,ses = run_trajectory(cfg,oob,MS)
open(expanduser("~/drayage/MSTP/rerun/figdata/convergence.json"),"w") do io
  JSON3.write(io, Dict("iters"=>iters,"lb"=>lbs,"ub"=>ubs,"se"=>ses,
    "instance"=>"capreg_l20u95","topology"=>"6x6","tau"=>cfg.tau,"rho_cross"=>0.0,"noob"=>length(oob)))
end
println("DONE wrote convergence.json")
