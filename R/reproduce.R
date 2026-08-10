# Reproduction harness. Turns a list of experiment specifications into a table
# of SDDP lower/upper bounds and certified gaps, deterministically (seeded). One
# harness covers the paper's bound tables — scalability, horizons, sensitivity,
# and single-instance convergence — which differ only in which spec fields vary:
#   * scalability : vary (nOrigins, nDestinations, nCarriers, iters)
#   * horizons    : vary tau
#   * sensitivity : vary lambda / rho_cross / spot_mult
#   * convergence : fix the instance, vary iters
#
# Because R owns all randomness (R/scenarios.R) and the engine has none, a spec
# list + seed reproduces the same numbers bit-for-bit. These are a NEW, seeded
# baseline: they are not expected to match the old unseeded rerun/*.jl oracle
# digit-for-digit (see REPRODUCE.md), only to agree with it statistically.

# Return `x` unless NULL, else `default`.
.or <- function(x, default) if (is.null(x)) default else x

#' Reproduce SDDP bound tables from a list of experiment specifications
#'
#' For each spec, generates a seed-deterministic instance, builds the copula
#' correlation matrix, trains an SDDP policy (seeded), and certifies the lower
#' bound against a fresh out-of-sample estimate of the policy cost. Returns one
#' row per spec with the lower bound, out-of-sample mean (upper bound), certified
#' gap, and validity flag.
#'
#' @param specs A list of specifications. Each element is a named list with:
#'   `tau`, `nOrigins`, `nDestinations`, `nCarriers`, `lambda`, `iters`, and
#'   optionally `label`, `trials` (default 200), `n_scenarios` (default 10),
#'   `rho_oo` (0.6), `rho_dd` (0.7), `rho_cross` (0.0), and `spot_mult` (1) —
#'   a multiplier applied to the instance spot-cost coefficients.
#' @param seed Integer base seed; the same seed reproduces the whole table.
#' @param z Standard-error margin for the bound-validity flag (default 3).
#' @return A `data.frame` with columns `label`, `tau`, `iters`, `LB`, `UB`,
#'   `gap_pct` (= 100*(UB-LB)/|LB|), `sim_se`, `margin_se`, `valid`.
#' @export
mstp_reproduce_bounds <- function(specs, seed = 1L, z = 3.0) {
  .ensure_engine()
  stopifnot(is.list(specs), length(specs) >= 1L)

  rows <- lapply(specs, function(s) {
    stopifnot(!is.null(s$tau), !is.null(s$nOrigins), !is.null(s$nDestinations),
              !is.null(s$nCarriers), !is.null(s$lambda), !is.null(s$iters))

    inst <- generate_instance(tau = as.integer(s$tau),
                              nOrigins = as.integer(s$nOrigins),
                              nDestinations = as.integer(s$nDestinations),
                              nCarriers = as.integer(s$nCarriers),
                              seed = as.integer(seed),
                              lambda = as.numeric(s$lambda))
    if (!is.null(s$spot_mult))
      inst$spot_coef <- inst$spot_coef * as.numeric(s$spot_mult)

    S <- mstp_corrmat(s$nOrigins, s$nDestinations,
                      rho_oo = .or(s$rho_oo, 0.6),
                      rho_dd = .or(s$rho_dd, 0.7),
                      rho_cross = .or(s$rho_cross, 0.0))
    cfg <- mstp_config(inst, lambda = as.numeric(s$lambda), corrmat = S,
                       n_scenarios = as.integer(.or(s$n_scenarios, 10L)))

    model <- mstp_train(cfg, iterations = as.integer(s$iters), seed = as.integer(seed))
    vb <- mstp_validate_bound(model, cfg, trials = as.integer(.or(s$trials, 200L)),
                              z = z, seed = as.integer(seed))

    data.frame(
      label     = .or(s$label, sprintf("%dx%dx%d", s$nOrigins, s$nDestinations, s$nCarriers)),
      tau       = as.integer(s$tau),
      iters     = as.integer(s$iters),
      LB        = vb$bound,
      UB        = vb$sim_mean,
      gap_pct   = 100 * (vb$sim_mean - vb$bound) / abs(vb$bound),
      sim_se    = vb$sim_se,
      margin_se = vb$margin_se,
      valid     = vb$valid,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
