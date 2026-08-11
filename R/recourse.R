# Recourse analyses — regret against the perfect-foresight LP and gain of
# recourse against the myopic policy. Both consume the SDDP simulation output
# (mstp_simulate: per-trajectory obj + realised noise) and the pure-R transport
# LPs in R/lp_transport.R (no TLPR, no Julia).

# Initial inventory state of an instance, as the LP builders expect it.
.init_state <- function(inst) {
  list(entry = as.numeric(inst$entry_stock_0),
       exitp = as.numeric(inst$exit_stock_0),
       exitm = as.numeric(inst$exit_short_0))
}

# Apply `f` over `items`, optionally in parallel. The closure `f` carries its
# captured data to workers automatically; `export` names the internal LP
# functions it calls, which are shipped from the environment where the analysis
# functions are defined (package namespace when installed, globalenv when
# sourced). n_cores = 1 runs serially (no cluster overhead).
.apply_maybe_parallel <- function(items, f, n_cores, export) {
  if (is.null(n_cores)) n_cores <- 1L
  if (n_cores <= 1L) return(vapply(items, f, numeric(1)))
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, export, envir = environment(.apply_maybe_parallel))
  as.numeric(parallel::parLapply(cl, items, f))
}

#' Compute out-of-sample regret against the perfect-foresight LP
#'
#' For each instance, takes the SDDP simulation objectives and realised noise
#' from [mstp_simulate()], solves the clairvoyant (perfect-foresight) LP for each
#' trajectory, and returns relative regret,
#' \eqn{(\mathrm{SDDP} - \mathrm{LP}^*) / \mathrm{LP}^*}.
#'
#' @param sims List of instances (from [generate_instance()] or
#'   [load_instances()]).
#' @param sddp_results List of [mstp_simulate()] result lists, one per instance
#'   in the same order as `sims`.
#' @param n_cores Number of parallel workers (default 1 = serial).
#' @return Numeric matrix of shape `(n_trials x n_instances)` of regret values.
#' @export
compute_regret <- function(sims, sddp_results, n_cores = 1L) {
  n_instances <- length(sims)
  stopifnot(length(sddp_results) == n_instances)
  regrets <- vector("list", n_instances)

  for (idx in seq_len(n_instances)) {
    inst  <- sims[[idx]]
    lay   <- .transport_layout(inst)
    init  <- .init_state(inst)
    obj   <- sddp_results[[idx]]$obj
    noise <- sddp_results[[idx]]$noise
    tau <- inst$tau; nI <- inst$nOrigins; nJ <- inst$nDestinations

    f <- function(i) {
      rows  <- (i - 1L) * tau + seq_len(tau)
      clair <- .clairvoyant_cost(lay,
                                 noise[rows, seq_len(nI), drop = FALSE],
                                 noise[rows, nI + seq_len(nJ), drop = FALSE], init)
      (obj[i] - clair) / clair
    }
    regrets[[idx]] <- .apply_maybe_parallel(
      seq_along(obj), f, n_cores,
      export = c(".clairvoyant_cost", ".spot_obj", "solve_lp"))
  }
  do.call(cbind, regrets)
}

#' Compute gain of recourse: SDDP vs myopic policy
#'
#' For each instance and each out-of-sample trajectory from [mstp_simulate()],
#' rolls a myopic (greedy single-stage) policy forward on the same realised noise
#' and compares its cost against the SDDP policy cost. Both costs use the SDDP
#' stage accounting (inst/julia/model.jl), so the comparison is like-for-like.
#'
#' Gain of recourse \eqn{= (\mathrm{myopic} - \mathrm{SDDP}) / \mathrm{myopic}}.
#'
#' @param sims List of instances (from [generate_instance()] or
#'   [load_instances()]).
#' @param sddp_results List of [mstp_simulate()] result lists, one per instance.
#' @param n_cores Number of parallel workers (default 1 = serial).
#' @return Numeric matrix of shape `(n_trials x n_instances)` of gain values.
#' @export
compute_vss <- function(sims, sddp_results, n_cores = 1L) {
  n_instances <- length(sims)
  stopifnot(length(sddp_results) == n_instances)
  gains <- vector("list", n_instances)

  for (idx in seq_len(n_instances)) {
    inst  <- sims[[idx]]
    lay   <- .transport_layout(inst)
    init  <- .init_state(inst)
    obj   <- sddp_results[[idx]]$obj
    noise <- sddp_results[[idx]]$noise
    tau <- inst$tau; nI <- inst$nOrigins; nJ <- inst$nDestinations

    f <- function(i) {
      rows <- (i - 1L) * tau + seq_len(tau)
      myop <- .myopic_cost(lay,
                           noise[rows, seq_len(nI), drop = FALSE],
                           noise[rows, nI + seq_len(nJ), drop = FALSE], init)
      (myop - obj[i]) / myop
    }
    gains[[idx]] <- .apply_maybe_parallel(
      seq_along(obj), f, n_cores,
      export = c(".myopic_cost", ".stage_decision", ".stage_cost",
                 ".spot_obj", "solve_lp"))
  }
  do.call(cbind, gains)
}

#' Plot regret density curves across instances and iteration counts
#'
#' @param regrets_list Named list of regret matrices (name = iteration count,
#'   e.g. `list(SDDPx500 = mat, SDDPx1000 = mat)`).
#' @param max_regret X-axis upper limit.
#' @export
plot_regret <- function(regrets_list, max_regret = NULL) {
  n <- length(regrets_list)
  old_par <- graphics::par(mfrow = c(1L, n))
  on.exit(graphics::par(old_par), add = TRUE)

  for (key in names(regrets_list)) {
    iters     <- as.numeric(gsub("SDDPx", "", key))
    densities <- apply(regrets_list[[key]], 1L, stats::density)
    ylim      <- range(sapply(densities, function(d) d$y))
    xlim      <- range(sapply(densities, function(d) d$x))
    if (!is.null(max_regret)) xlim[2L] <- min(xlim[2L], max_regret)

    graphics::plot(NA, xlim = xlim, ylim = pmin(ylim, 5e5),
         xlab = "Regret", ylab = "Density",
         main = paste0(iters, " SDDP Iterations"),
         xaxt = "n", yaxt = "n", cex.lab = 1.5, cex.main = 1.5)
    for (d in densities)
      graphics::lines(d, col = grDevices::rgb(0, 0, 0, alpha = 0.1), lwd = 2L)
  }
}
