# MSTP

R package for solving the **Multistage Stochastic Transportation Problem** arising in drayage procurement, as described in:

> Vassos, Lusby & Pinson (2025). *Multistage Stochastic Transportation Problem via SDDP.* arXiv:2505.01813

The package provides an R interface to a [SDDP.jl](https://github.com/odow/SDDP.jl) Julia engine. Julia runs as a transparent backend via [JuliaCall](https://non-contradiction.github.io/JuliaCall/); all data passes in memory — no JSON or HDF5 intermediaries.

---

## Problem

A shipper procures drayage capacity across a bipartite origin–destination network over a rolling 12-period horizon. At each period, uncertain container inflows at origins and outflows at destinations must be routed through a mix of contracted carriers (won through a combinatorial auction) and spot-market carriers. The goal is to minimise total holding, shortage, and transport costs in expectation.

Uncertainty is modelled as correlated Poisson flows (Gaussian copula), and the multistage stochastic program is solved with Stochastic Dual Dynamic Programming (SDDP).

---

## Installation

### Prerequisites

1. **Julia ≥ 1.9** — [julialang.org](https://julialang.org/downloads/)
2. **R packages**

```r
install.packages(c("JuliaCall", "jsonlite", "Matrix", "highs",
                   "data.table", "ggplot2", "parallel"))
```

3. **TLPR** — local LP construction library (must be installed separately):

```r
devtools::install_local("~/drayage/TLPR")
```

4. **Install this package**

```r
devtools::install_local("~/drayage/MSTP")
```

### Julia packages (first time only)

```r
library(MSTP)
setup_engine(install = TRUE)   # installs SDDP.jl, HiGHS.jl, etc.
```

Subsequent sessions only need `setup_engine()`.

---

## Quick start

```r
library(MSTP)

# 1. Start the Julia engine
setup_engine()

# 2. Generate (or load) an instance
inst <- generate_instance(
  tau           = 12L,   # time periods
  nOrigins      = 6L,
  nDestinations = 6L,
  nCarriers     = 20L,
  seed          = 1L,
  lambda        = 50.0   # Poisson demand rate; calibrates carrier capacities
)

# 3. Build the SDDP configuration
config <- mstp_config(inst, lambda = 50.0)

# 4. Train the policy (returns a Julia proxy)
model <- mstp_train(config, iterations = 300L)

# 5. Simulate out-of-bag performance
sims <- mstp_simulate(model, config, trials = 1000L)

cat("Mean cost:", mean(sims$obj), "\n")
```

Loading saved instances from JSON files:

```r
instances <- load_instances("path/to/instances/", n_instances = 100L)
```

---

## Workflow

```
generate_instance()  ──►  mstp_config()  ──►  mstp_train()  ──►  mstp_simulate()
       │                        ▲                   │                    │
  (or load JSON)                │             mstp_write_cuts()          │
                                │                   │             ┌──────┴──────────┐
                    mstp_update_capacity()          ▼             │                 │
                                │           mstp_train_warm()  compute_regret()  compute_vss()
                                │                   │             plot_regret()
                                └── mstp_capacity_duals() ◄───────┘
                                    (projected gradient loop)
```

| Step | Function | Description |
|------|----------|-------------|
| Instance | `generate_instance()` | Random MSTP instance; capacities calibrated to `lambda` at `target_util` utilisation |
| Instance | `load_instances()` | Load batch from JSON files |
| Config | `mstp_config()` | Build Julia `HyperParams` from an instance |
| Config | `mstp_gen_corrmat()` | Generate block-diagonal correlation matrix |
| Train | `mstp_train()` | Train SDDP policy (returns opaque Julia proxy) |
| Train | `mstp_train_warm()` | Rebuild model with new config, load saved cuts, continue training |
| Train | `mstp_write_cuts()` | Persist Benders cuts to JSON for warm-starting |
| Simulate | `mstp_simulate()` | Out-of-bag simulation; returns costs, inventory, allocations |
| Cap. opt. | `mstp_capacity_duals()` | Extract ∂V/∂carrier_capacity from SDDP constraint duals |
| Cap. opt. | `mstp_update_capacity()` | Return instance copy with updated `carrier_capacity` vector |
| Analysis | `compute_regret()` | SDDP cost vs clairvoyant LP — how far from perfect information |
| Analysis | `compute_vss()` | SDDP cost vs myopic policy — gain of recourse |
| Analysis | `sensitivity_inflow()` | Optimal cost distribution under random inflow draws |
| Analysis | `sensitivity_spot()` | Optimal cost distribution under random spot rate draws |
| Analysis | `plot_sensitivity()` | Density plots for sensitivity analysis |
| Analysis | `plot_regret()` | Density plots for regret distributions |
| Logs | `parse_logs()` | Parse SDDP log files (bias, CI ratio, wall time) |
| Logs | `summarise_logs()` | Summary statistics from parsed logs |

---

## Key concepts

**Gain of recourse** (`compute_vss`): compares the adaptive SDDP policy against a myopic (greedy) policy that optimises each period in isolation. A gain of ~27% means the stochastic policy cuts costs by roughly a quarter relative to myopic dispatch.

**Regret** (`compute_regret`): compares SDDP cost against the clairvoyant LP that knows all future realisations. This bounds how far the learned policy is from the theoretical optimum. On 6×6×20 instances (τ=12, λ=50) mean regret is below 1.5% at 300 iterations (~90 s), and stable from 250 iterations onward — making regret a reliable stopping criterion.

**SDDP convergence** (`parse_logs` / `summarise_logs`): the simulation CI narrows with iterations. On the 6×6×20 benchmark: training completes in ~90 s at 300 iterations; scalability experiments reach 40×40×100 instances in under 15 minutes at 20 iterations using the open-source HiGHS solver.

**Capacity optimisation** (`mstp_capacity_duals` / `mstp_update_capacity` / `mstp_write_cuts` / `mstp_train_warm`): first-stage carrier capacity can be optimised via projected subgradient descent. The dual of each carrier capacity constraint (`∂V/∂carrier_capacity`) is averaged over simulated trajectories and used as the gradient. Benders cuts are functions of inventory state variables (not of `carrier_capacity`), so they transfer correctly across capacity changes — `mstp_write_cuts` / `mstp_train_warm` exploit this to warm-start each outer iteration cheaply. On a 6×6×20 instance (τ=12, λ=50), 15 gradient steps reduce total procurement cost by 8.4% relative to default capacities — versus a 32.3% cost *increase* from a deterministic LP proxy, a 40.7 pp advantage attributable to variance-aware gradient estimation.

---

## Architecture

```
R (MSTP package)
├── engine.R        Julia session management and R↔Julia API wrappers
├── geninst.R       Instance generation
├── instances.R     Batch JSON loading, R-side correlation matrix
├── sensitivity.R   Inflow and spot-rate sensitivity analysis (parallel LP)
├── recourse.R      Regret (clairvoyant LP) and VSS (myopic policy)
├── logs.R          SDDP log parsing and summarisation
└── utils.R         HiGHS LP adapter, TLPR model builders

Julia (inst/julia/)
├── mstp.jl         Entry point (sources all modules)
├── types.jl        HyperParams struct (plain types, R-compatible)
├── model.jl        SDDP stage subproblem (JuMP)
├── utils.jl        Correlated Poisson sampling
├── api.jl          build_config / train_model (batched, crash-resilient) / train_model_warm / write_cuts / simulate_model / simulate_cap_duals / batch_run
└── setup.jl        Julia package installation
```

Julia objects (`HyperParams`, `SDDP.PolicyGraph`) are held as opaque proxies in R and passed back when needed — no serialisation overhead.

---

## Dependencies

| Package | Role |
|---------|------|
| JuliaCall | R↔Julia bridge |
| SDDP.jl | Stochastic dual dynamic programming |
| HiGHS / HiGHS.jl | LP/MIP solver (open-source) |
| TLPR | Local LP construction utilities |
| Matrix | Sparse matrix support |
| parallel | Multi-core LP batch solving |
| jsonlite | JSON instance I/O |
