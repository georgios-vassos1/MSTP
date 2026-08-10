# Deterministic scenario generation for the MSTP SDDP engine.
#
# This module is the single source of randomness for the package. R owns all
# scenario draws; the Julia layer consumes them as fixed data (no Julia-side
# RNG). Given a seed, every quantity downstream is reproducible.
#
# The uncertainty model is a Gaussian copula with Poisson marginals, mirroring
# the reference Julia pipeline N(0, Sigma) -> Phi -> Poisson^{-1}
# (inst/julia/utils.jl:24-26): draw latent normals with the supplied correlation
# matrix, map through the standard-normal CDF, then invert per-dimension Poisson
# CDFs. CDF values are clamped away from {0, 1} so Poisson quantiles stay finite.

# CDF clamp bound; matches `_PCLAMP` in inst/julia/utils.jl:20. Gaussian tails
# beyond ~+/-8 sigma round to 0/1 in floating point and would yield infinite
# Poisson quantiles (LP infeasibility).
.PCLAMP <- 1e-10

# Evaluate `code` under a temporary RNG seed without disturbing the caller's
# global RNG stream. If `seed` is NULL, `code` runs against the current stream.
# Restores `.Random.seed` (or removes it if it did not exist) on exit.
.with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  stopifnot(is.numeric(seed), length(seed) == 1L)
  has_old <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_old) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      },
      add = TRUE
    )
  }
  set.seed(as.integer(seed))
  force(code)
}

#' Block-structured Gaussian-copula correlation matrix
#'
#' Builds the `(nOrigins + nDestinations)`-dimensional correlation matrix used by
#' the copula sampler. The origin block carries a common off-diagonal
#' correlation `rho_oo`, the destination block `rho_dd`, and every origin-vs-
#' destination (cross) entry `rho_cross`. The diagonal is 1.
#'
#' This is the single source of truth for the correlation structure. Unlike the
#' legacy Julia `gen_cov_mat()` (inst/julia/utils.jl:10), the within-block values
#' are explicit arguments rather than random draws, matching the fixed values
#' documented in the paper (rho_OO = 0.6, rho_DD = 0.7).
#'
#' @param nOrigins Integer >= 1. Number of origin (entry) locations.
#' @param nDestinations Integer >= 1. Number of destination (exit) locations.
#' @param rho_oo Numeric in [0, 1). Origin-origin off-diagonal correlation.
#' @param rho_dd Numeric in [0, 1). Destination-destination off-diagonal correlation.
#' @param rho_cross Numeric in [0, 1). Origin-destination (cross-block) correlation.
#'
#' @return A symmetric `d x d` numeric matrix (`d = nOrigins + nDestinations`)
#'   with unit diagonal. Verified positive definite before return.
#' @export
mstp_corrmat <- function(nOrigins, nDestinations,
                         rho_oo = 0.6, rho_dd = 0.7, rho_cross = 0.0) {
  stopifnot(
    is.numeric(nOrigins), length(nOrigins) == 1L, nOrigins >= 1,
    is.numeric(nDestinations), length(nDestinations) == 1L, nDestinations >= 1,
    is.numeric(rho_oo), length(rho_oo) == 1L, rho_oo >= 0, rho_oo < 1,
    is.numeric(rho_dd), length(rho_dd) == 1L, rho_dd >= 0, rho_dd < 1,
    is.numeric(rho_cross), length(rho_cross) == 1L, rho_cross >= 0, rho_cross < 1
  )
  nOrigins <- as.integer(nOrigins)
  nDestinations <- as.integer(nDestinations)
  d <- nOrigins + nDestinations

  Sigma <- matrix(rho_cross, nrow = d, ncol = d)
  o <- seq_len(nOrigins)
  dd <- nOrigins + seq_len(nDestinations)
  Sigma[o, o] <- rho_oo
  Sigma[dd, dd] <- rho_dd
  diag(Sigma) <- 1.0

  # Positive definiteness is required by the latent MVN draw; fail loudly rather
  # than let MASS::mvrnorm warn and silently distort the copula.
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= .Machine$double.eps^0.5) {
    stop(sprintf(
      "correlation matrix is not positive definite (min eigenvalue %.3g); reduce rho_cross",
      min(ev)
    ), call. = FALSE)
  }
  Sigma
}

#' Sample correlated-Poisson scenario vectors (Gaussian copula)
#'
#' Draws `n` independent scenario vectors from a Gaussian copula with Poisson
#' marginals. Mirrors `sample_scenarios()` in inst/julia/utils.jl:24-26.
#'
#' @param n Integer >= 1. Number of scenario vectors to draw.
#' @param lambda Numeric, length 1 or `d`. Poisson rate(s); length-1 is recycled
#'   to all `d` dimensions. All entries must be > 0.
#' @param corrmat Symmetric `d x d` correlation matrix (see [mstp_corrmat()]).
#' @param seed Optional integer. If supplied, draws are seeded and the caller's
#'   global RNG stream is left untouched (see Details).
#'
#' @return An integer matrix with `n` rows and `d` columns. Column order is
#'   `[inflow_1..nOrigins, outflow_1..nDestinations]`, matching
#'   inst/julia/model.jl:50-51. Entries are non-negative.
#'
#' @details Randomness is scoped: with `seed` non-NULL the function saves,
#'   sets, and restores `.Random.seed`, so repeated calls with the same seed are
#'   identical and no global RNG side effect leaks to the caller.
#' @export
mstp_sample_copula <- function(n, lambda, corrmat, seed = NULL) {
  stopifnot(
    is.numeric(n), length(n) == 1L, n >= 1,
    is.numeric(lambda), all(lambda > 0),
    is.matrix(corrmat), nrow(corrmat) == ncol(corrmat)
  )
  n <- as.integer(n)
  d <- nrow(corrmat)
  if (length(lambda) == 1L) lambda <- rep(lambda, d)
  if (length(lambda) != d) {
    stop(sprintf("length(lambda) must be 1 or %d, got %d", d, length(lambda)),
         call. = FALSE)
  }

  .with_seed(seed, {
    z <- MASS::mvrnorm(n = n, mu = rep(0, d), Sigma = corrmat)
    z <- matrix(z, nrow = n, ncol = d)          # normalise shape when n == 1
    u <- stats::pnorm(z)
    u <- pmin(pmax(u, .PCLAMP), 1 - .PCLAMP)     # clamp like inst/julia/utils.jl:25
    lam <- matrix(lambda, nrow = n, ncol = d, byrow = TRUE)
    x <- stats::qpois(u, lam)
    matrix(as.integer(x), nrow = n, ncol = d)
  })
}

#' Per-stage in-sample support for SDDP training
#'
#' Builds the finite noise support Omega that the SDDP stage subproblems are
#' parameterised over. Under stagewise independence each of the `tau` stages
#' gets an independent draw of `n_scenarios` outcomes, reproducing the per-stage
#' sampling in inst/julia/model.jl:48 but with the draws owned by R.
#'
#' @param tau Integer >= 1. Number of stages (planning horizon).
#' @param n_scenarios Integer >= 1. Support size per stage.
#' @param lambda,corrmat As in [mstp_sample_copula()].
#' @param seed Optional integer seed (scoped; see [mstp_sample_copula()]).
#'
#' @return A list of length `tau`; element `t` is an `n_scenarios x d` integer
#'   matrix of noise outcomes for stage `t`. Equal probabilities are implied.
#' @export
mstp_scenario_support <- function(tau, n_scenarios, lambda, corrmat, seed = NULL) {
  stopifnot(
    is.numeric(tau), length(tau) == 1L, tau >= 1,
    is.numeric(n_scenarios), length(n_scenarios) == 1L, n_scenarios >= 1
  )
  tau <- as.integer(tau)
  n_scenarios <- as.integer(n_scenarios)
  .with_seed(seed, {
    lapply(seq_len(tau), function(t) {
      mstp_sample_copula(n_scenarios, lambda, corrmat, seed = NULL)
    })
  })
}

#' Sample full-horizon scenario paths for Historical sampling
#'
#' Draws `n_paths` independent length-`tau` trajectories from the copula, for use
#' as the deterministic forward-pass or out-of-sample scheme fed to
#' `SDDP.Historical` in the Julia layer.
#'
#' @param n_paths Integer >= 1. Number of trajectories.
#' @param tau Integer >= 1. Stages per trajectory.
#' @param lambda,corrmat As in [mstp_sample_copula()].
#' @param seed Optional integer seed (scoped; see [mstp_sample_copula()]).
#'
#' @return A list of length `n_paths`; each element is a `tau x d` integer matrix
#'   whose row `t` is the stage-`t` noise vector.
#' @export
mstp_scenario_paths <- function(n_paths, tau, lambda, corrmat, seed = NULL) {
  stopifnot(
    is.numeric(n_paths), length(n_paths) == 1L, n_paths >= 1,
    is.numeric(tau), length(tau) == 1L, tau >= 1
  )
  n_paths <- as.integer(n_paths)
  tau <- as.integer(tau)
  .with_seed(seed, {
    lapply(seq_len(n_paths), function(p) {
      mstp_sample_copula(tau, lambda, corrmat, seed = NULL)
    })
  })
}
