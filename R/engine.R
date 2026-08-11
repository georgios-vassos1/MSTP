# ─── Julia session management ─────────────────────────────────────────────────

#' Initialise the Julia engine (called automatically on first use)
#'
#' Starts the Julia session and loads the MSTP engine. Run `setup_engine()`
#' once per R session; subsequent calls are no-ops. To install Julia packages
#' for the first time, call `setup_engine(install = TRUE)`.
#'
#' @param install If TRUE, runs `inst/julia/setup.jl` to add Julia packages.
#' @param ... Passed to `JuliaCall::julia_setup()`.
#' @export
setup_engine <- function(install = FALSE, ...) {
  if (.mstp$loaded) return(invisible(NULL))

  JuliaCall::julia_setup(...)

  if (install) {
    setup_script <- system.file("julia/setup.jl", package = "MSTP")
    JuliaCall::julia_source(setup_script)
  }

  entry <- system.file("julia/mstp.jl", package = "MSTP")
  JuliaCall::julia_source(entry)

  .mstp$loaded <- TRUE
  invisible(NULL)
}

.ensure_engine <- function() {
  if (!.mstp$loaded) setup_engine()
}

# ─── Config ───────────────────────────────────────────────────────────────────

#' Build an MSTP configuration (HyperParams) from an instance list
#'
#' Passes instance parameters to Julia and returns an opaque Julia proxy
#' representing the `HyperParams` struct. The proxy can be passed to
#' `mstp_train()` and is otherwise invisible to the user.
#'
#' @param instance  Named list as returned by `generate_instance()`.
#' @param lambda    Poisson intensity for the uncertainty model. Either a
#'                  single numeric (broadcast to all flow dimensions) or a
#'                  numeric vector of length `nOrigins + nDestinations`.
#' @param corrmat   Correlation matrix of size `(nOrigins+nDestinations)^2`.
#'                  Use `mstp_corrmat()` to generate one.
#' @param n_scenarios Number of training scenarios (support size per stage).
#' @return An object of class `mstp_config`: a list with the opaque Julia proxy
#'   `jl` (the `HyperParams`) plus the R-side sampling parameters (`tau`, `dim`,
#'   `lambda`, `corrmat`, `n_scenarios`) that [mstp_train()] and [mstp_simulate()]
#'   use to draw scenarios in R.
#' @export
mstp_config <- function(instance, lambda = 2000.0, corrmat = NULL,
                         n_scenarios = 20L) {
  .ensure_engine()

  # Normalise field names from old JSON instances (lowercase winners/bids)
  if (is.null(instance$Winners) && !is.null(instance$winners))
    instance$Winners <- instance$winners
  if (is.null(instance$Bids) && !is.null(instance$bids))
    instance$Bids <- instance$bids
  if (is.null(instance$nSpotCarriers))
    instance$nSpotCarriers <- instance$nCarriers

  dim <- as.integer(instance$nOrigins + instance$nDestinations)
  if (is.null(corrmat)) corrmat <- diag(dim)
  lambda <- if (length(lambda) == 1L) rep(as.numeric(lambda), dim) else as.numeric(lambda)
  stopifnot(
    length(lambda) == dim,
    is.matrix(corrmat), nrow(corrmat) == dim, ncol(corrmat) == dim
  )

  params <- c(
    instance[c("tau", "nOrigins", "nDestinations", "nCarriers", "nSpotCarriers",
               "Bids", "Winners",
               "entry_stock_0", "exit_stock_0", "exit_short_0",
               "entry_capacity", "exit_capacity",
               "entry_store_coef", "exit_store_coef", "exit_short_coef",
               "transport_coef", "spot_coef", "carrier_capacity")],
    list(
      lambda      = lambda,
      corrmat     = corrmat,
      n_scenarios = as.integer(n_scenarios)
    )
  )

  structure(
    list(
      jl          = JuliaCall::julia_call("build_config", params),
      tau         = as.integer(instance$tau),
      dim         = dim,
      lambda      = lambda,
      corrmat     = corrmat,
      n_scenarios = as.integer(n_scenarios)
    ),
    class = "mstp_config"
  )
}

# The correlation matrix is built in R by mstp_corrmat() (R/scenarios.R), the
# single source of truth for the copula dependence structure. The former
# mstp_gen_corrmat() wrapper around Julia's RNG-based gen_cov_mat was removed:
# R owns all randomness, and mstp_corrmat() takes explicit, documented
# within-block correlations instead of drawing them.

# ─── Train ────────────────────────────────────────────────────────────────────

#' Train an SDDP policy
#'
#' Builds and trains the multistage stochastic transportation model using
#' SDDP. Returns an opaque Julia proxy for the trained `PolicyGraph`, which
#' can be passed to `mstp_simulate()`.
#'
#' @param config      `mstp_config` object returned by `mstp_config()`.
#' @param iterations  Number of SDDP training iterations (default 1500). One
#'   forward trajectory is drawn per iteration.
#' @param seed        Optional integer. When supplied, the per-stage support and
#'   the forward trajectories are drawn reproducibly (from derived sub-streams
#'   `seed` and `seed + 1`), so the trained policy — and its lower bound — are
#'   deterministic.
#' @return An opaque Julia proxy (SDDP.PolicyGraph).
#' @export
mstp_train <- function(config, iterations = 1500L, seed = NULL) {
  .ensure_engine()
  stopifnot(inherits(config, "mstp_config"))
  iterations <- as.integer(iterations)

  Om  <- mstp_scenario_support(config$tau, config$n_scenarios, config$lambda,
                               config$corrmat, seed = seed)
  fwd <- mstp_scenario_paths(iterations, config$tau, config$lambda,
                             config$corrmat,
                             seed = if (is.null(seed)) NULL else seed + 1L)

  JuliaCall::julia_call(
    "train_flat", config$jl,
    .flatten_support(Om), config$tau, config$n_scenarios, config$dim,
    .flatten_paths(fwd), iterations, iterations
  )
}

# ─── Bound ────────────────────────────────────────────────────────────────────

#' Return the SDDP lower bound after training
#'
#' Calls \code{SDDP.calculate_bound(model)} on the Julia side and returns the
#' dual / Benders lower bound on the optimal expected cost.
#'
#' @param model Julia proxy returned by \code{mstp_train()}.
#' @return Numeric scalar — SDDP lower bound.
#' @export
mstp_bound <- function(model) {
  .ensure_engine()
  as.numeric(JuliaCall::julia_call("get_bound", model))
}

#' Check whether the SDDP lower bound is numerically valid
#'
#' A valid SDDP lower bound for a minimisation satisfies
#' \code{bound <= E[policy cost]}. This retrieves the lower bound, runs an
#' out-of-sample simulation to estimate \code{E[policy cost]}, and applies the
#' R decision rule [mstp_bound_validity()] (valid iff the bound does not exceed
#' the simulated mean by more than \code{z} standard errors). Use it to detect
#' the ghost-cut overshoot (LB > UB) that arises from ill-conditioned LP duals
#' on long horizons / large instances, so an invalid bound is never reported
#' silently.
#'
#' @param model   Julia proxy returned by \code{mstp_train()}.
#' @param config  \code{mstp_config} object returned by \code{mstp_config()}.
#' @param trials  Number of out-of-sample simulations (default 200).
#' @param z       Standard-error margin allowed before flagging (default 3).
#' @param seed    Optional integer for reproducible evaluation scenarios.
#' @return A named list: \code{bound}, \code{sim_mean}, \code{sim_se},
#'   \code{margin_se} (SEs the bound sits above the mean; positive = overshoot),
#'   \code{trials}, and \code{valid} (logical).
#' @export
mstp_validate_bound <- function(model, config, trials = 200L, z = 3.0, seed = NULL) {
  .ensure_engine()
  stopifnot(inherits(config, "mstp_config"))
  bound <- mstp_bound(model)
  sim   <- mstp_simulate(model, config, trials = as.integer(trials), seed = seed)
  r     <- mstp_bound_validity(bound, sim$obj, z = as.numeric(z))
  if (!isTRUE(r$valid)) {
    warning(sprintf(
      "SDDP lower bound exceeds simulated mean by %.2f SE (bound=%.1f, sim_mean=%.1f +/- %.1f SE): numerically invalid (ghost-cut overshoot).",
      r$margin_se, bound, r$mean, r$se), call. = FALSE)
  }
  list(bound = bound, sim_mean = r$mean, sim_se = r$se,
       margin_se = r$margin_se, trials = r$n, valid = r$valid)
}

# ─── Cuts / warm start ──────────────────────────────────────────────────────

#' Write SDDP Benders cuts to a file
#'
#' Persists the cuts from a trained model so they can be reloaded into a new
#' model via `mstp_train_warm()`. Cuts are keyed by stage and state-variable
#' names, so they remain valid when only `carrier_capacity` changes between
#' calls (the inventory state variables are unchanged).
#'
#' @param model Julia proxy returned by `mstp_train()` or `mstp_train_warm()`.
#' @param path  File path for the cuts JSON (created or overwritten).
#' @return `path`, invisibly.
#' @export
mstp_write_cuts <- function(model, path) {
  .ensure_engine()
  JuliaCall::julia_call("write_cuts", model, as.character(path))
  invisible(path)
}

#' Train an SDDP policy with a warm start from existing cuts
#'
#' Builds a new `PolicyGraph` from `config`, loads Benders cuts from
#' `cuts_path`, then continues training. Because cuts are functions of the
#' inventory state variables (not of `carrier_capacity`), they transfer
#' correctly across capacity changes and dramatically reduce the iterations
#' needed for the lower bound to tighten.
#'
#' @param config      `mstp_config` object returned by `mstp_config()`.
#' @param iterations  Additional SDDP iterations to run after loading cuts.
#' @param cuts_path   Path to a cuts file written by `mstp_write_cuts()`.
#' @param seed        Optional integer seed for the scenarios drawn in R
#'   (support and forward trajectories).
#' @return An opaque Julia proxy (SDDP.PolicyGraph).
#' @export
mstp_train_warm <- function(config, iterations, cuts_path, seed = NULL) {
  .ensure_engine()
  stopifnot(inherits(config, "mstp_config"))
  iterations <- as.integer(iterations)
  Om  <- mstp_scenario_support(config$tau, config$n_scenarios, config$lambda,
                               config$corrmat, seed = seed)
  fwd <- mstp_scenario_paths(iterations, config$tau, config$lambda,
                             config$corrmat,
                             seed = if (is.null(seed)) NULL else seed + 1L)
  JuliaCall::julia_call(
    "train_warm_flat", config$jl,
    .flatten_support(Om), config$tau, config$n_scenarios, config$dim,
    .flatten_paths(fwd), iterations, iterations, as.character(cuts_path)
  )
}

#' Extract expected carrier capacity duals from a trained SDDP model
#'
#' Simulates `n_samples` out-of-sample trajectories (drawn in R) and averages
#' the dual variable of each carrier capacity constraint at every stage. The
#' returned vector has the same flat layout as `instance$carrier_capacity` —
#' index `(k-1)*tau + t` for carrier `k` at stage `t` — so it can be used
#' directly as the gradient `∂V/∂carrier_capacity` in a capacity optimisation
#' loop.
#'
#' Because duals of a minimisation `<=` constraint are non-positive, the
#' gradient of `V(x) + v·x` with respect to `x` is `cap_duals + v`.
#'
#' @param model     Julia proxy returned by `mstp_train()`.
#' @param config    `mstp_config` object returned by `mstp_config()`.
#' @param n_samples Number of simulation trajectories used to average the duals.
#' @param seed      Optional integer for reproducible trajectories (derived
#'   sub-stream `seed + 4`).
#' @return Numeric vector of length `(nCarriers + nSpotCarriers) * tau`.
#' @export
mstp_capacity_duals <- function(model, config, n_samples = 100L, seed = NULL) {
  .ensure_engine()
  stopifnot(inherits(config, "mstp_config"))
  oob <- mstp_scenario_paths(as.integer(n_samples), config$tau, config$lambda,
                             config$corrmat,
                             seed = if (is.null(seed)) NULL else seed + 4L)
  result <- JuliaCall::julia_call("cap_duals_flat", model, config$jl,
                                   .flatten_paths(oob), as.integer(n_samples),
                                   config$tau, config$dim)
  as.numeric(result$duals)
}

#' Optimise carrier capacities: LP proxy vs SDDP dual gradient (CRN)
#'
#' Runs the capacity-optimisation comparison in the Julia core: projected
#' gradient descent from the default capacities using (a) the duals of a
#' mean-demand LP proxy and (b) the SDDP stage capacity duals, then scores all
#' three capacity vectors (default, LP proxy, SDDP) on a shared out-of-sample
#' scenario set (common random numbers). All scenarios are sampled in R, so the
#' result is reproducible for a given `seed`.
#'
#' @param config `mstp_config` object.
#' @param v Per-unit reservation price (e.g. `0.2 * mean(transport rate)`).
#' @param outer Outer projected-gradient iterations (default 12).
#' @param cold SDDP iterations per gradient step for the SDDP-dual source
#'   (cold start; default 100).
#' @param eval_iters SDDP iterations for the out-of-sample scoring policy
#'   (default 300).
#' @param n_eval Shared out-of-sample scenarios for scoring (default 500).
#' @param n_capdual Trajectories used to average the SDDP capacity duals
#'   (default 200).
#' @param seed Optional integer base seed (derived sub-streams).
#' @return A named list: total procurement costs (`x0_total`, `lp_total`,
#'   `sddp_total`), percentage change vs default (`lp_pct`, `sddp_pct`), the SDDP
#'   advantage in percentage points (`adv_pp`), and mean capacities per plan.
#' @export
mstp_capacity_optimize <- function(config, v, outer = 12L, cold = 100L,
                                   eval_iters = 300L, n_eval = 500L,
                                   n_capdual = 200L, seed = NULL) {
  .ensure_engine()
  stopifnot(inherits(config, "mstp_config"), is.numeric(v), length(v) == 1L)
  outer <- as.integer(outer); cold <- as.integer(cold); eval_iters <- as.integer(eval_iters)
  n_capdual <- as.integer(n_capdual); n_eval <- as.integer(n_eval)
  n_fwd <- max(cold, eval_iters)
  s <- function(off) if (is.null(seed)) NULL else as.integer(seed) + off

  Om  <- mstp_scenario_support(config$tau, config$n_scenarios, config$lambda,
                               config$corrmat, seed = s(0L))
  fwd <- mstp_scenario_paths(n_fwd,     config$tau, config$lambda, config$corrmat, seed = s(1L))
  cap <- mstp_scenario_paths(n_capdual, config$tau, config$lambda, config$corrmat, seed = s(2L))
  ev  <- mstp_scenario_paths(n_eval,    config$tau, config$lambda, config$corrmat, seed = s(3L))

  res <- JuliaCall::julia_call(
    "capopt_regime_flat", config$jl,
    .flatten_support(Om), .flatten_paths(fwd), .flatten_paths(cap), .flatten_paths(ev),
    config$tau, config$n_scenarios, config$dim, n_fwd, n_capdual, n_eval,
    v = as.numeric(v), outer = outer, cold = cold, eval_iters = eval_iters)
  lapply(res, as.numeric)
}

#' Record the SDDP-dual projected-gradient capacity trajectory (Figure 5 data)
#'
#' Runs the SDDP-dual projected-gradient capacity optimisation and records, per
#' outer iteration, the total procurement cost, gradient norm, and capacity step
#' size. All scenarios are sampled in R (reproducible for a given `seed`).
#'
#' @param config `mstp_config` object.
#' @param v Per-unit reservation price.
#' @param outer,cold,eval_iters,n_eval,n_capdual As in [mstp_capacity_optimize()].
#' @param seed Optional integer base seed.
#' @return A named list of numeric vectors (length `outer`): `objective`,
#'   `grad_norm`, `step_norm`.
#' @export
mstp_capacity_trajectory <- function(config, v, outer = 15L, cold = 100L,
                                     eval_iters = 300L, n_eval = 500L,
                                     n_capdual = 200L, seed = NULL) {
  .ensure_engine()
  stopifnot(inherits(config, "mstp_config"), is.numeric(v), length(v) == 1L)
  outer <- as.integer(outer); cold <- as.integer(cold); eval_iters <- as.integer(eval_iters)
  n_capdual <- as.integer(n_capdual); n_eval <- as.integer(n_eval)
  n_fwd <- max(cold, eval_iters)
  s <- function(off) if (is.null(seed)) NULL else as.integer(seed) + off

  Om  <- mstp_scenario_support(config$tau, config$n_scenarios, config$lambda,
                               config$corrmat, seed = s(0L))
  fwd <- mstp_scenario_paths(n_fwd,     config$tau, config$lambda, config$corrmat, seed = s(1L))
  cap <- mstp_scenario_paths(n_capdual, config$tau, config$lambda, config$corrmat, seed = s(2L))
  ev  <- mstp_scenario_paths(n_eval,    config$tau, config$lambda, config$corrmat, seed = s(3L))

  res <- JuliaCall::julia_call(
    "capopt_traj_flat", config$jl,
    .flatten_support(Om), .flatten_paths(fwd), .flatten_paths(cap), .flatten_paths(ev),
    config$tau, config$n_scenarios, config$dim, n_fwd, n_capdual, n_eval,
    v = as.numeric(v), outer = outer, cold = cold, eval_iters = eval_iters)
  lapply(res, as.numeric)
}

#' Return a copy of an instance with updated carrier capacities
#'
#' Replaces `instance$carrier_capacity` with the supplied vector (rounded to
#' integers) and returns the modified instance list. Use this to build a new
#' `mstp_config()` at a candidate capacity `x` during optimisation.
#'
#' @param instance Named list as returned by `generate_instance()` or
#'   `load_instances()`.
#' @param carrier_capacity Numeric vector of length
#'   `(nCarriers + nSpotCarriers) * tau` — the new flat capacity vector.
#' @return A copy of `instance` with `carrier_capacity` replaced.
#' @export
mstp_update_capacity <- function(instance, carrier_capacity) {
  instance$carrier_capacity <- as.integer(round(carrier_capacity))
  instance
}

# ─── Simulate ─────────────────────────────────────────────────────────────────

#' Simulate the trained policy out-of-bag
#'
#' Draws `trials` out-of-sample scenario paths in R (reproducible when `seed` is
#' supplied) and runs the trained SDDP policy through them. Returns an R list
#' with objective values, noise realisations, inventory trajectories, and
#' carrier allocations — all as plain numeric matrices ready for analysis.
#'
#' @param model   Julia proxy returned by `mstp_train()`.
#' @param config  `mstp_config` object returned by `mstp_config()`.
#' @param trials  Number of OOB simulation replications (default 1000).
#' @param seed    Optional integer for reproducible OOB scenarios (derived
#'   sub-stream `seed + 2`).
#' @return A named list:
#'   \describe{
#'     \item{obj}{Numeric vector of length `trials` — total cost per run.}
#'     \item{noise}{Matrix `(trials*tau) × (nOrigins+nDestinations)` —
#'                  realised inflows/outflows.}
#'     \item{entry}{Matrix `(trials*tau) × nOrigins` — entry inventory.}
#'     \item{exitp}{Matrix `(trials*tau) × nDestinations` — exit stock.}
#'     \item{exitm}{Matrix `(trials*tau) × nDestinations` — exit shortage.}
#'     \item{moves}{Matrix `(trials*tau) × (nLanes+nSpotLanes)` — allocations.}
#'     \item{tau, trials, nOrigins, nDestinations, nLanes, nSpotLanes}{Dims.}
#'   }
#' @export
mstp_simulate <- function(model, config, trials = 1000L, seed = NULL) {
  .ensure_engine()
  stopifnot(inherits(config, "mstp_config"))
  oob <- mstp_scenario_paths(as.integer(trials), config$tau, config$lambda,
                             config$corrmat,
                             seed = if (is.null(seed)) NULL else seed + 2L)
  result <- JuliaCall::julia_call("simulate_flat", model, config$jl,
                                   .flatten_paths(oob), as.integer(trials),
                                   config$tau, config$dim)

  # Convert scalar dims to plain R integers for convenience
  result$tau           <- as.integer(result$tau)
  result$trials        <- as.integer(result$trials)
  result$nOrigins      <- as.integer(result$nOrigins)
  result$nDestinations <- as.integer(result$nDestinations)
  result$nLanes        <- as.integer(result$nLanes)
  result$nSpotLanes    <- as.integer(result$nSpotLanes)

  result
}
