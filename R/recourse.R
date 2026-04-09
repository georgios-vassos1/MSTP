#' Compute out-of-bag regret for a set of instances
#'
#' For each instance, takes the SDDP simulation objectives and noise realisations
#' produced by `mstp_simulate()`, solves the clairvoyant LP for each scenario,
#' and returns relative regret: (SDDP_cost - LP_optimal) / LP_optimal.
#'
#' @param sims       List of instance environments (from `load_instances()`).
#' @param sddp_results List of simulation result lists from `mstp_simulate()`,
#'                   one per instance (same order as `sims`).
#' @param n_cores    Number of parallel workers (default: all but 2).
#' @return Numeric matrix of shape `(n_trials × n_instances)` — regret values.
#' @export
compute_regret <- function(sims, sddp_results, n_cores = NULL) {
  n_instances <- length(sims)
  stopifnot(length(sddp_results) == n_instances)
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 2L)

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  regrets <- vector("list", n_instances)

  for (idx in seq_len(n_instances)) {
    env        <- build_env(sims[[idx]])
    init_state <- c(sims[[idx]]$entry_stock_0,
                    c(rbind(sims[[idx]]$exit_stock_0, sims[[idx]]$exit_short_0)))
    model      <- build_model(env, init_state,
                              sims[[idx]]$Q[-(env$tau * env$nI + env$I_)],
                              sims[[idx]]$D[-(env$tau * env$nJ + env$J_)])

    taux    <- env$tau + 1L
    sddp_obj <- sddp_results[[idx]]$obj
    sddp_ksi <- sddp_results[[idx]]$noise
    N        <- length(sddp_obj)

    pdx    <- env$nI + 2L * env$nJ + env$nCS + env$nCO
    offset <- 2L * (env$nI + env$nJ) + env$nCS + env$nCO
    pdx    <- pdx + offset * 0L:(env$tau - 1L)
    Qdx.1  <- c(outer(env$I_, pdx, "+"))
    Qdx.2  <- c(env$nI + env$nJ + Qdx.1)
    Ddx    <- Qdx.1 + env$nI

    job <- function(i) {
      Q <- c(t(sddp_ksi[(i - 1L) * taux + 1L:taux, env$I_]))
      D <- c(t(sddp_ksi[(i - 1L) * taux + 1L:taux, env$nI + env$J_]))
      model$rhs[Qdx.1] <- model$rhs[Qdx.2] <- Q[-(env$tau * env$nI + env$I_)]
      model$rhs[Ddx]   <- D[-(env$tau * env$nJ + env$J_)]
      opt <- solve_lp(model)
      (sddp_obj[i] - opt$objbound) / opt$objbound
    }

    parallel::clusterExport(cl, c("model", "env", "taux", "N",
                                   "Qdx.1", "Qdx.2", "Ddx",
                                   "sddp_obj", "sddp_ksi", "solve_lp"),
                             envir = environment())
    regrets[[idx]] <- parallel::parLapply(cl, seq(N), job) |> as.numeric()
  }

  do.call(cbind, regrets)
}

#' Plot regret density curves across instances and iteration counts
#'
#' @param regrets_list Named list of regret matrices (name = iteration count,
#'                     e.g. `list(SDDPx500 = mat, SDDPx1000 = mat, ...)`).
#' @param max_regret   X-axis upper limit.
#' @export
plot_regret <- function(regrets_list, max_regret = NULL) {
  n <- length(regrets_list)
  old_par <- par(mfrow = c(1L, n))
  on.exit(par(old_par), add = TRUE)

  for (key in names(regrets_list)) {
    iters     <- as.numeric(gsub("SDDPx", "", key))
    densities <- apply(regrets_list[[key]], 1L, density)
    ylim      <- range(sapply(densities, function(d) d$y))
    xlim      <- range(sapply(densities, function(d) d$x))
    if (!is.null(max_regret)) xlim[2L] <- min(xlim[2L], max_regret)

    plot(NA, xlim = xlim, ylim = pmin(ylim, 5e5),
         xlab = "Regret", ylab = "Density",
         main = paste0(iters, " SDDP Iterations"),
         xaxt = "n", yaxt = "n", cex.lab = 1.5, cex.main = 1.5)

    for (d in densities)
      lines(d, col = rgb(0, 0, 0, alpha = 0.1), lwd = 2L)

    rngx <- seq(xlim[1L], xlim[2L], length.out = 20L)
    axis(1L, at = rngx[1L:9L], labels = round(rngx * 1e6, 1L)[1L:9L], cex.axis = 1.5)
    mtext(expression("x" ~ 10^-6), side = 1L, line = 1L, at = par("usr")[2L])

    rngy <- seq(0, min(ylim[2L], 5e5), length.out = 6L)
    axis(2L, at = rngy, labels = rngy / 1e6, cex.axis = 1.5)
    mtext(expression("x" ~ 10^6), side = 3L, line = 1L, at = par("usr")[1L])
  }
}
