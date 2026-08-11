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
                              lambda = as.numeric(s$lambda),
                              target_util = .or(s$target_util, 0.8))
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

#' Reproduce the recourse results (regret and gain of recourse)
#'
#' Generates a seed-deterministic instance, trains an SDDP policy, simulates it
#' out-of-sample, and reports regret against the perfect-foresight LP
#' (Table 2 in the paper) and gain of recourse against the myopic policy. All
#' analysis is pure R (no TLPR).
#'
#' @param spec A named list as in [mstp_reproduce_bounds()] (`tau`, `nOrigins`,
#'   `nDestinations`, `nCarriers`, `lambda`, `iters`, optionally `label`,
#'   `trials`, `n_scenarios`, `rho_oo`, `rho_dd`, `rho_cross`).
#' @param seed Integer base seed (reproduces the whole row).
#' @return A one-row `data.frame` with `label`, `mean_regret_pct`,
#'   `median_regret_pct`, `q95_regret_pct`, and `mean_gain_pct`.
#' @export
mstp_reproduce_recourse <- function(spec, seed = 1L) {
  .ensure_engine()
  stopifnot(!is.null(spec$tau), !is.null(spec$nOrigins), !is.null(spec$nDestinations),
            !is.null(spec$nCarriers), !is.null(spec$lambda), !is.null(spec$iters))

  inst <- generate_instance(tau = as.integer(spec$tau),
                            nOrigins = as.integer(spec$nOrigins),
                            nDestinations = as.integer(spec$nDestinations),
                            nCarriers = as.integer(spec$nCarriers),
                            seed = as.integer(seed), lambda = as.numeric(spec$lambda))
  S <- mstp_corrmat(spec$nOrigins, spec$nDestinations,
                    rho_oo = .or(spec$rho_oo, 0.6), rho_dd = .or(spec$rho_dd, 0.7),
                    rho_cross = .or(spec$rho_cross, 0.0))
  cfg   <- mstp_config(inst, lambda = as.numeric(spec$lambda), corrmat = S,
                       n_scenarios = as.integer(.or(spec$n_scenarios, 10L)))
  model <- mstp_train(cfg, iterations = as.integer(spec$iters), seed = as.integer(seed))
  sim   <- mstp_simulate(model, cfg, trials = as.integer(.or(spec$trials, 200L)),
                         seed = as.integer(seed))

  reg <- as.numeric(compute_regret(list(inst), list(sim)))
  gn  <- as.numeric(compute_vss(list(inst), list(sim)))

  data.frame(
    label             = .or(spec$label, sprintf("%dx%dx%d", spec$nOrigins, spec$nDestinations, spec$nCarriers)),
    mean_regret_pct   = 100 * mean(reg),
    median_regret_pct = 100 * stats::median(reg),
    q95_regret_pct    = 100 * as.numeric(stats::quantile(reg, 0.95)),
    mean_gain_pct     = 100 * mean(gn),
    stringsAsFactors  = FALSE
  )
}

#' Reproduce a capacity-optimisation regime (Tables 7-8)
#'
#' Generates a seed-deterministic instance at a given demand/utilisation regime
#' and runs [mstp_capacity_optimize()]. The reservation price is fixed at
#' `0.2 * mean(contract transport rate)`, matching the paper.
#'
#' @param spec Named list: `tau`, `nOrigins`, `nDestinations`, `nCarriers`,
#'   `lambda`, and optionally `label`, `target_util` (default 0.8), `n_scenarios`
#'   (20), `rho_cross` (0.0), `rho_oo` (0.6), `rho_dd` (0.7), and the optimisation
#'   controls `outer` (12), `cold` (100), `eval_iters` (300), `n_eval` (500),
#'   `n_capdual` (200).
#' @param seed Integer base seed.
#' @return A one-row `data.frame`: `label`, `lp_pct`, `sddp_pct`, `adv_pp`,
#'   `x0_mean`, `xlp_mean`, `xsd_mean`.
#' @export
mstp_reproduce_capopt <- function(spec, seed = 1L) {
  .ensure_engine()
  stopifnot(!is.null(spec$tau), !is.null(spec$nOrigins), !is.null(spec$nDestinations),
            !is.null(spec$nCarriers), !is.null(spec$lambda))

  inst <- generate_instance(tau = as.integer(spec$tau),
                            nOrigins = as.integer(spec$nOrigins),
                            nDestinations = as.integer(spec$nDestinations),
                            nCarriers = as.integer(spec$nCarriers),
                            seed = as.integer(seed), lambda = as.numeric(spec$lambda),
                            target_util = .or(spec$target_util, 0.8))
  S   <- mstp_corrmat(spec$nOrigins, spec$nDestinations,
                      rho_oo = .or(spec$rho_oo, 0.6), rho_dd = .or(spec$rho_dd, 0.7),
                      rho_cross = .or(spec$rho_cross, 0.0))
  cfg <- mstp_config(inst, lambda = as.numeric(spec$lambda), corrmat = S,
                     n_scenarios = as.integer(.or(spec$n_scenarios, 20L)))
  v <- 0.2 * mean(inst$transport_coef)

  r <- mstp_capacity_optimize(cfg, v = v,
                              outer = as.integer(.or(spec$outer, 12L)),
                              cold = as.integer(.or(spec$cold, 100L)),
                              eval_iters = as.integer(.or(spec$eval_iters, 300L)),
                              n_eval = as.integer(.or(spec$n_eval, 500L)),
                              n_capdual = as.integer(.or(spec$n_capdual, 200L)),
                              seed = seed)
  data.frame(
    label = .or(spec$label, sprintf("%dx%dx%d", spec$nOrigins, spec$nDestinations, spec$nCarriers)),
    lp_pct = r$lp_pct, sddp_pct = r$sddp_pct, adv_pp = r$adv_pp,
    x0_mean = r$x0_mean, xlp_mean = r$xlp_mean, xsd_mean = r$xsd_mean,
    stringsAsFactors = FALSE
  )
}
