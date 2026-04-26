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
#'                  Use `mstp_gen_corrmat()` to generate one.
#' @param n_scenarios Number of training scenarios sampled per SDDP iteration.
#' @return An opaque Julia proxy (HyperParams).
#' @export
mstp_config <- function(instance, lambda = 2000.0, corrmat = NULL,
                         n_scenarios = 10L) {
  .ensure_engine()

  # Normalise field names from old JSON instances (lowercase winners/bids)
  if (is.null(instance$Winners) && !is.null(instance$winners))
    instance$Winners <- instance$winners
  if (is.null(instance$Bids) && !is.null(instance$bids))
    instance$Bids <- instance$bids
  if (is.null(instance$nSpotCarriers))
    instance$nSpotCarriers <- instance$nCarriers

  if (is.null(corrmat)) {
    nOD     <- instance$nOrigins + instance$nDestinations
    corrmat <- diag(nOD)
  }

  params <- c(
    instance[c("tau", "nOrigins", "nDestinations", "nCarriers", "nSpotCarriers",
               "Bids", "Winners",
               "entry_stock_0", "exit_stock_0", "exit_short_0",
               "entry_capacity", "exit_capacity",
               "entry_store_coef", "exit_store_coef", "exit_short_coef",
               "transport_coef", "spot_coef", "carrier_capacity")],
    list(
      lambda      = as.numeric(lambda),
      corrmat     = corrmat,
      n_scenarios = as.integer(n_scenarios)
    )
  )

  JuliaCall::julia_call("build_config", params)
}

#' Generate a correlation matrix for the uncertainty model
#'
#' Wraps Julia's `gen_cov_mat`. Produces a block-diagonal covariance matrix
#' with intra-block correlations drawn uniformly in [0.6, 0.8] and
#' cross-block correlation set to `cross_corr`.
#'
#' @param n_blocks   Number of blocks (typically 2: one for origins, one for
#'                   destinations).
#' @param block_size Number of nodes per block (nOrigins or nDestinations).
#' @param cross_corr Off-diagonal cross-block correlation (default 0.4).
#' @return A numeric matrix of size `(n_blocks*block_size)^2`.
#' @export
mstp_gen_corrmat <- function(n_blocks, block_size, cross_corr = 0.4) {
  .ensure_engine()
  JuliaCall::julia_call("gen_cov_mat",
                         as.integer(n_blocks),
                         as.integer(block_size),
                         as.numeric(cross_corr))
}

# ─── Train ────────────────────────────────────────────────────────────────────

#' Train an SDDP policy
#'
#' Builds and trains the multistage stochastic transportation model using
#' SDDP. Returns an opaque Julia proxy for the trained `PolicyGraph`, which
#' can be passed to `mstp_simulate()`.
#'
#' @param config      Julia proxy returned by `mstp_config()`.
#' @param iterations  Number of SDDP training iterations (default 1500).
#' @return An opaque Julia proxy (SDDP.PolicyGraph).
#' @export
mstp_train <- function(config, iterations = 1500L) {
  .ensure_engine()
  JuliaCall::julia_call("train_model", config, as.integer(iterations))
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

# ─── Simulate ─────────────────────────────────────────────────────────────────

#' Simulate the trained policy out-of-bag
#'
#' Generates fresh OOB scenarios and runs the trained SDDP policy through
#' them. Returns an R list with objective values, noise realisations, inventory
#' trajectories, and carrier allocations — all as plain numeric matrices ready
#' for analysis.
#'
#' @param model   Julia proxy returned by `mstp_train()`.
#' @param config  Julia proxy returned by `mstp_config()`.
#' @param trials  Number of OOB simulation replications (default 1000).
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
#' Extract expected carrier capacity duals from a trained SDDP model
#'
#' Runs `n_samples` forward simulations and collects the dual variable of each
#' carrier capacity constraint at every stage. The returned vector has the same
#' flat layout as `instance$carrier_capacity` — index `(k-1)*tau + t` for
#' carrier `k` at stage `t` — so it can be used directly as the gradient
#' `∂V/∂carrier_capacity` in a capacity optimisation loop.
#'
#' Because duals of a minimisation ≤ constraint are non-positive, the gradient
#' of `V(x) + v·x` with respect to `x` is `cap_duals + v`.
#'
#' @param model     Julia proxy returned by `mstp_train()`.
#' @param config    Julia proxy returned by `mstp_config()`.
#' @param n_samples Number of simulation trajectories used to average the duals.
#' @return Numeric vector of length `(nCarriers + nSpotCarriers) * tau`.
#' @export
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
#' @param config      Julia proxy returned by `mstp_config()`.
#' @param iterations  Additional SDDP iterations to run after loading cuts.
#' @param cuts_path   Path to a cuts file written by `mstp_write_cuts()`.
#' @return An opaque Julia proxy (SDDP.PolicyGraph).
#' @export
mstp_train_warm <- function(config, iterations, cuts_path) {
  .ensure_engine()
  JuliaCall::julia_call("train_model_warm", config, as.integer(iterations),
                         as.character(cuts_path))
}

mstp_capacity_duals <- function(model, config, n_samples = 100L) {
  .ensure_engine()
  result <- JuliaCall::julia_call("simulate_cap_duals", model, config,
                                   as.integer(n_samples))
  as.numeric(result$duals)
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

mstp_simulate <- function(model, config, trials = 1000L) {
  .ensure_engine()
  result <- JuliaCall::julia_call("simulate_model", model, config,
                                   as.integer(trials))

  # Convert scalar dims to plain R integers for convenience
  result$tau           <- as.integer(result$tau)
  result$trials        <- as.integer(result$trials)
  result$nOrigins      <- as.integer(result$nOrigins)
  result$nDestinations <- as.integer(result$nDestinations)
  result$nLanes        <- as.integer(result$nLanes)
  result$nSpotLanes    <- as.integer(result$nSpotLanes)

  result
}
