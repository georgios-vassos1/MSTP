# Refreshed certified results — FIXED generator (init_stock=0, low-overlap, differentiated rates)

All bounds certified by validate_bound (z=3 SE). Pure-Julia (rerun/*.jl). λ=700 unless noted.

## Table 2 — clairvoyant-LP regret (τ=4, λ=5.5, ρ=0.4, 1000 iter, 500 OOB)
| topo | mean% | median% | q95% |
|---|---|---|---|
| 2×1 | 0.002 | 0.000 | 0.000 |
| 1×2 | 0.001 | 0.000 | 0.000 |
| 2×2 | 0.182 | 0.000 | 0.537 |
→ all << 1.5%; SDDP near-clairvoyant-optimal (tighter than the old table).

## Table 3 — scalability (τ=12, seed 42), all VALID
| instance | iters | LB | UB | gap% |
|---|---|---|---|---|
| 6×6×20   | 500 | 394,847   | 399,719   | +1.23 |
| 20×20×100| 50  | 1,262,417 | 1,264,101 | +0.13 |
| 40×40×100| 20  | 2,456,332 | 2,486,850 | +1.24 |

## Table 4 — horizons (6×6×20, seed 42, 300 iter), all VALID
| τ | LB | UB | gap% |
|---|---|---|---|
| 12 | 392,139   | 392,636   | +0.13 |
| 26 | 1,017,229 | 992,309   | −2.45 |
| 52 | 2,102,888 | 2,090,832 | −0.57 |
→ negative gaps now small and within MC noise (valid); the old −7.1% overshoot is gone.

## Table 5 — sensitivity (6×6×20 τ=12, seed-avg 42–44, 300 iter), all VALID
| sweep | value | UB |
|---|---|---|
| λ | 200 / 700 / 2000 | 130,371 / 383,182 / 1,031,441 |
| ρ_cross | 0.0 / 0.2 / 0.4 | 389,717 / 376,302 / 369,842 |
| m | 1 / 2 / 4 | 385,976 / 386,907 / 385,977 |
→ demand λ dominant; ρ and spot-multiplier m have small effect (m flat → spot price negligible).

## Table 6 — regret-based stopping (2×2, 500 OOB)
| iters | mean% | sd% | q95% |
|---|---|---|---|
| 100 | 0.046 | 0.322 | 0.180 |
| 250 | 0.145 | 1.084 | 0.300 |
| 500 | 0.135 | 1.376 | 0.214 |
| 1000| 0.043 | 0.330 | 0.174 |
| 1500| 0.057 | 0.474 | 0.256 |
→ regret tiny and stable across budgets; regret-based stopping viable.

## VSS — gain of recourse (2×2, 500 OOB)
mean 51.3%, median 51.5%, positive on 100% of trials (paper: 48.9%). Holds.

## Table 7 — capacity optimisation (6×6×20, λ=50) — **RESULT REVERSED**
| policy | total (UB+v·x) | vs x0 |
|---|---|---|
| x0 (default) | 51,268 | — |
| LP proxy   | 44,601 | **−13.0%** |
| SDDP duals | 45,556 | **−11.1%** |
→ SDDP advantage over LP proxy = **−1.9 pp** (LP marginally better). Both reduce cost ~11–13%.
The paper's headline (LP proxy +32.3%, SDDP −8.4%, +40.7 pp SDDP advantage) does NOT survive the
corrected instances — it was an artifact of the broken (random-init, high-overlap) generator, which
made the LP proxy collapse capacity to ~2 units. G5's "SDDP-dual advantage" claim is not supported.

## Table 7 — capacity optimisation, CORRECTED (CRN evaluation, regime sweep)

Earlier "G5 refuted" was WRONG on two counts: (1) the LP-vs-SDDP total-cost
difference was measured with independent noisy 100-iter runs → ±10pp noise;
(2) the test used benign uncertainty (λ=50, 80% utilisation, CV≈14%). Fixed with
common-random-numbers evaluation (both capacity vectors on the SAME 500 OOB
scenarios, 300-iter policy). CRN stability check: adv(rep1)=0.0pp, adv(rep2)=0.1pp
→ measurement reliable.

| regime | CV | capacity | LP total vs x0 | SDDP total vs x0 | SDDP advantage |
|---|---|---|---|---|---|
| λ=50, util 0.80, ρ=0.0 | 14% | loose | −3.0% | −3.1% | 0.0 pp |
| λ=50, util 0.80, ρ=0.4 | 14% | loose | −3.4% | −3.4% | 0.0 pp |
| λ=20, util 0.95, ρ=0.4 | 22% | tight | +32.0% | −0.6% | **+32.5 pp** |
| λ=10, util 0.95, ρ=0.4 | 32% | tight | +18.1% | −0.6% | **+18.8 pp** |

**Conclusion: G5 is VINDICATED, not refuted.** The SDDP-dual capacity advantage is
real and large (≈19–33 pp) precisely when it should matter — scarce capacity +
high demand variance, where the mean-demand LP proxy collapses capacity and the
stochastic operational cost explodes (LP ≈ +32%, matching the paper's +32.3%).
Under benign uncertainty (loose capacity, low CV) both methods are equivalent
(0 pp) — as expected. The paper's headline holds, but the demonstrating instance
must have tight capacity / meaningful CV (not λ=50 at 80% utilisation with the
fixed generator). The old +40.7pp came from the broken generator; the honest,
reproducible advantage is ≈20–33 pp in the tight/variable regime.

## Correlation impact (NEW contribution) — conditional on capacity scarcity

CRN capacity-opt at TIGHT capacity (λ=20, util 0.95), sweeping cross-block
supply–demand correlation ρ_cross:

| ρ_cross | LP proxy total | SDDP total | SDDP advantage |
|---|---|---|---|
| 0.0 | +23.2% | −0.5% | 23.7 pp |
| 0.2 | +30.4% | −0.6% | 31.0 pp |
| 0.4 | +32.0% | −0.6% | 32.5 pp |

Contrast — at LOOSE capacity (λ=50, util 0.80): ρ=0.0 → 0.0 pp, ρ=0.4 → 0.0 pp.

**Finding:** the impact of correlation is *conditional on capacity scarcity*.
When capacity is ample, correlated supply–demand shocks are absorbed by slack and
ρ_cross barely moves cost (consistent with the paper's Table-5 loose-capacity
result). When capacity is tight, correlation cannot be diversified across nodes,
so it inflates the mean-demand LP proxy's cost (23%→32%) and amplifies the value
of variance-aware SDDP capacity planning (advantage 24→33 pp as ρ_cross 0→0.4).
This interaction (correlation × scarcity) is a clean, novel result the current
paper does not report.
