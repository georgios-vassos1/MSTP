# Reproducing the paper results

The R-first MSTP package regenerates the paper's SDDP results as a **new, seeded,
deterministic baseline**. This document says what is reproduced, how, and how it
relates to the frozen `rerun/` oracle.

## Determinism and the oracle

R owns all randomness (`R/scenarios.R`); the Julia engine has none. A spec list
plus a seed therefore reproduces the same numbers bit-for-bit — verified by
`tests/integration/reproduce_baseline.R` (same seed → identical table) and
`tests/integration/test_engine_determinism.R` (same seed → identical LB and
out-of-sample cost vector).

These numbers are **not** expected to match the `rerun/*.jl` oracle digit-for-digit:

- The oracle was produced by the **pre-refactor** engine with **unseeded** Julia
  RNG (see `rerun/RESULTS.md`; "seed 42" there is an instance-file label, not an
  RNG seed) and the default `InSampleMonteCarlo` training scheme.
- The R-first engine uses R-sampled support and `SDDP.Historical` training, so it
  is a different (deterministic) sampling regime.
- The oracle scripts are also engine-incompatible since the M2 refactor (they call
  the old `train_model` / `validate_bound` signatures). They are retained as
  **frozen provenance** — their captured numbers (`rerun/RESULTS.md`, and the
  gitignored `rerun/*.out`) are the regression reference — and are not re-run.

Parity with the oracle is therefore **statistical / qualitative**: valid bounds
(`LB <= UB` within sampling error), small certified gaps, mean regret below the
paper's threshold, positive gain of recourse — not digit equality.

## Bound tables (scalability, horizons, sensitivity, convergence)

One harness, `mstp_reproduce_bounds()` (`R/reproduce.R`), produces all four bound
tables; they differ only in which spec fields vary:

| Paper table | Vary across specs |
|---|---|
| Scalability (Table 3) | `nOrigins`/`nDestinations`/`nCarriers`, `iters` |
| Horizons (Table 4)    | `tau` |
| Sensitivity (Table 5) | `lambda`, `rho_cross`, `spot_mult` |
| Convergence (Table 6) | `iters` (instance fixed) |

Each row reports `LB`, `UB` (out-of-sample mean), `gap_pct = 100*(UB-LB)/|LB|`,
`sim_se`, `margin_se`, and `valid` (the R rule `mstp_bound_validity`, z = 3 SE).

Paper-scale example (seed 42; slow — minutes to tens of minutes; not run in CI):

```r
specs <- list(
  list(label="6x6x20",    tau=12L, nOrigins=6L,  nDestinations=6L,  nCarriers=20L,  lambda=700, iters=500L),
  list(label="20x20x100", tau=12L, nOrigins=20L, nDestinations=20L, nCarriers=100L, lambda=700, iters=50L),
  list(label="40x40x100", tau=12L, nOrigins=40L, nDestinations=40L, nCarriers=100L, lambda=700, iters=20L)
)
mstp_reproduce_bounds(specs, seed = 42L)
```

Verified here at tiny (fast) scale, `tests/integration/reproduce_baseline.R`
(seed 7): both a τ=4 and a τ=6 `2x2x3` instance reproduce identically across runs
and are certified valid (margins ≈ 2.5–2.7 SE). Bound validity depends on support
size: coarse support (`n_scenarios=5`) on a tiny instance can overshoot and is
correctly flagged invalid; the paper's `n_scenarios=10`+ on larger instances is
cleaner.

## Regret (Table 2) and gain of recourse / VSS

`mstp_reproduce_recourse()` (`R/reproduce.R`) generates a seed-deterministic
instance, trains and simulates the SDDP policy, and reports regret against the
perfect-foresight LP and gain of recourse against the myopic policy. Both
analyses (`compute_regret`, `compute_vss`, `R/recourse.R`) are **pure R** — the
former TLPR-based clairvoyant/myopic LPs were reimplemented in `R/lp_transport.R`,
so the package no longer depends on TLPR.

Correctness is gated, not assumed (`tests/integration/`):

- `reconstruct_cost.R` — the R cost model reproduces the SDDP stage objective to
  ~1e-12 (validated on a non-square 2×3 instance so the spot-cost index `nL`
  is distinguished from `nOD`).
- `recourse_check.R` — on every out-of-sample trajectory the clairvoyant optimum
  is ≤ the SDDP realised cost (regret ≥ 0) and ≤ the myopic cost; `compute_regret`
  matches the direct computation; the gain-of-recourse baseline is locked.

Paper-scale example (seed 42; slow):

```r
mstp_reproduce_recourse(
  list(label = "2x2", tau = 4L, nOrigins = 2L, nDestinations = 2L,
       nCarriers = 4L, lambda = 5.5, iters = 1000L, trials = 500L,
       rho_cross = 0.4), seed = 42L)
```
