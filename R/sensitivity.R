#' Inflow sensitivity analysis
#'
#' For each instance in `sims`, solves the deterministic multi-period LP
#' under `n_samples` random inflow realisations drawn from `flow_support`.
#' Mirrors the analysis in the paper (Section 4).
#'
#' @param sims         List of instance environments (from `load_instances()`).
#' @param flow_support Integer vector — the inflow sample space.
#' @param n_samples    Number of inflow draws per instance.
#' @param n_cores      Number of parallel workers (default: all but 2).
#' @return Numeric matrix of shape `(n_samples × n_instances)` — optimal costs.
#' @export
sensitivity_inflow <- function(sims, flow_support = seq(1000L, 3000L, by = 100L),
                                n_samples = 1000L, n_cores = NULL) {
  n_instances <- length(sims)
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 2L)

  tau <- sims[[1L]]$tau
  nI  <- sims[[1L]]$nOrigins

  Q_mat <- sample(flow_support, tau * nI * n_samples, replace = TRUE) |>
    matrix(nrow = n_samples)

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  solve_lp <- solve_lp  # capture in local frame for cluster export

  optx <- numeric(n_instances * n_samples)

  for (idx in seq_len(n_instances)) {
    env        <- build_env(sims[[idx]])
    init_state <- c(sims[[idx]]$entry_stock_0,
                    c(rbind(sims[[idx]]$exit_stock_0, sims[[idx]]$exit_short_0)))
    model      <- build_model(env, init_state,
                              sims[[idx]]$Q[-(env$tau * env$nI + env$I_)],
                              sims[[idx]]$D[-(env$tau * env$nJ + env$J_)])

    pdx    <- env$nI + 2L * env$nJ + env$nCS + env$nCO
    offset <- 2L * (env$nI + env$nJ) + env$nCS + env$nCO
    pdx    <- pdx + offset * 0L:(env$tau - 1L)
    Qdx.1  <- c(outer(env$I_, pdx, "+"))
    Qdx.2  <- c(env$nI + env$nJ + Qdx.1)

    job <- function(k) {
      model$rhs[Qdx.1] <- model$rhs[Qdx.2] <- Q_mat[k, ]
      solve_lp(model)$objval
    }

    parallel::clusterExport(cl, c("model", "Qdx.1", "Qdx.2", "Q_mat",
                                   "solve_lp"), envir = environment())
    optx[(idx - 1L) * n_samples + seq(n_samples)] <-
      parallel::parLapply(cl, seq(n_samples), job) |> as.numeric()
  }

  matrix(optx, nrow = n_samples)
}

#' Spot rate sensitivity analysis
#'
#' For each instance in `sims`, solves the deterministic LP under `n_samples`
#' random spot rate realisations.
#'
#' @param sims      List of instance environments.
#' @param n_samples Number of spot rate draws per instance.
#' @param rate_min  Lower bound of uniform spot rate distribution.
#' @param rate_max  Upper bound of uniform spot rate distribution.
#' @param n_cores   Number of parallel workers.
#' @return Numeric matrix of shape `(n_samples × n_instances)`.
#' @export
sensitivity_spot <- function(sims, n_samples = 1000L,
                              rate_min = 3.0, rate_max = 9.0,
                              n_cores = NULL) {
  n_instances <- length(sims)
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 2L)

  nL_spot <- sims[[1L]]$nOrigins * sims[[1L]]$nDestinations *
             sims[[1L]]$nSpotCarriers * sims[[1L]]$tau

  spot_mat <- matrix(runif(nL_spot * n_samples, rate_min, rate_max), nrow = n_samples)

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  solve_lp <- solve_lp  # capture in local frame for cluster export

  optx <- numeric(n_instances * n_samples)

  for (idx in seq_len(n_instances)) {
    env        <- build_env(sims[[idx]])
    init_state <- c(sims[[idx]]$entry_stock_0,
                    c(rbind(sims[[idx]]$exit_stock_0, sims[[idx]]$exit_short_0)))
    model      <- build_model(env, init_state,
                              sims[[idx]]$Q[-(env$tau * env$nI + env$I_)],
                              sims[[idx]]$D[-(env$tau * env$nJ + env$J_)])

    offset   <- env$nI + 2L * env$nJ + env$nL_
    len_spot <- env$nCO * env$nL
    sdx      <- offset + outer(seq(len_spot), (offset + len_spot) * 0L:(env$tau - 1L), "+")

    job <- function(k) {
      model$obj[c(sdx)] <- spot_mat[k, ]
      solve_lp(model)$objval
    }

    parallel::clusterExport(cl, c("model", "sdx", "spot_mat", "solve_lp"),
                             envir = environment())
    optx[(idx - 1L) * n_samples + seq(n_samples)] <-
      parallel::parLapply(cl, seq(n_samples), job) |> as.numeric()
  }

  matrix(optx, nrow = n_samples)
}

#' Plot sensitivity density curves
#'
#' @param optx     Matrix from `sensitivity_inflow()` or `sensitivity_spot()`.
#' @param title    Plot title.
#' @param scale_x  Divisor for x-axis labels (default 1e6).
#' @export
plot_sensitivity <- function(optx, title = "Sensitivity", scale_x = 1e6) {
  densities <- apply(optx, 2L, density)
  ylim      <- range(sapply(densities, function(d) d$y))
  xlim      <- range(sapply(densities, function(d) d$x))

  plot(NA, xlim = xlim, ylim = ylim,
       xlab = "Total Cost", ylab = "Density", main = title,
       xaxt = "n", yaxt = "n")

  for (d in densities)
    lines(d, col = rgb(0, 0, 0, alpha = 0.2), lwd = 2L)

  rngx <- seq(xlim[1L], xlim[2L], length.out = 20L)
  axis(1L, at = rngx, labels = round(rngx / scale_x, 1L))
  rngy <- seq(0, ylim[2L], length.out = 4L)
  axis(2L, at = rngy, labels = round(rngy * scale_x, 2L))
  mtext(bquote(x ~ 10^{-6}), side = 3L, line = 1L, at = par("usr")[1L])
  mtext(bquote(x ~ 10^6),    side = 1L, line = 1L, at = par("usr")[2L])
}
