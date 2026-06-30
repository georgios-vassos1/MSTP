# Bound-overshoot acceptance gate (reviewer concern #3)

This is the slow, real-compute regression gate for the SDDP lower-bound
overshoot (LB > UB). It reproduces the bug on the **exact seed-42** 6×6×20 /
τ=52 / λ=700 instance with the 1e6 storage caps used in the production τ=52 run.

It is split into two steps to avoid the R+Julia `libomp` double-initialisation
crash (`OMP: Error #15`) seen when SDDP/HiGHS is driven through JuliaCall in the
same process as R:

1. `gen_seed42_instance.R` — pure R, writes the instance to JSON.
2. `check_bound_acceptance.jl` — pure Julia, loads the JSON, trains, and runs
   `validate_bound`. Exits non-zero if the bound is invalid.

## Run

```sh
# from the package root
Rscript tests/acceptance/gen_seed42_instance.R /tmp/inst52_seed42.json
julia  tests/acceptance/check_bound_acceptance.jl /tmp/inst52_seed42.json
```

## Status

**Currently RED.** Verified baseline (default `ContinuousConicDuality`, 300
iters, 200 trials): bound 4,066,030 vs sim mean 3,886,312 ± 28,847 →
margin **6.23 SE**, `valid=FALSE`. Matches the production result
(`TLPR/demo/results/26_results.rds`, τ=52, gap −7.1%).

A genuine fix must turn this GREEN (`valid=TRUE`, margin ≤ 3 SE). Note: tightening
the storage caps does NOT fix it (physical caps made it worse, 10.27 SE) — the
cause is the Benders cut slopes / LP-dual conditioning, not the rhs range.
