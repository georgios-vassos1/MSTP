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
  solve_lp <- solve_lp  # capture in local frame for cluster export

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

# Internal helper: solve one stage of the myopic (greedy) policy.
#
# Variables: [move(nMoves), entry_out(nI), exitp_out(nJ), exitm_out(nJ)]
# Objective: current-state holding costs (constant) + transport/spot costs
#            + next-state holding/shortage costs (entry_out, exitp_out, exitm_out).
# Including next-state costs gives the greedy policy an incentive to move
# containers and satisfy demand rather than just minimising transport spend.
#
# Constraints:
#   carrier capacity         (<=)
#   entry transition         (=)  entry_out + sum(move from i) = entry_in + Q_t
#   exit capacity            (<=) sum(move to j) <= exit_cap - exitp_in
#   exit balance             (=)  exitp_out - exitm_out - sum(move to j)
#                                   = exitp_in - exitm_in - D_t
.stage_lp <- function(t, entry_in, exitp_in, exitm_in, Q_t, D_t, env, sim) {
  nI      <- env$nI
  nJ      <- env$nJ
  nCS     <- env$nCS
  nCO     <- env$nCO
  nL      <- env$nL
  nL_     <- env$nL_
  nMoves  <- nL_ + nCO * nL
  nVars   <- nMoves + nI + 2L * nJ   # move + entry_out + exitp_out + exitm_out
  tau_sim <- env$tau + 1L            # = sim$tau; used for spot_coef indexing

  # Column offsets for state variables
  col_eo  <- nMoves + seq_len(nI)          # entry_out
  col_xp  <- nMoves + nI + seq_len(nJ)    # exitp_out
  col_xm  <- nMoves + nI + nJ + seq_len(nJ) # exitm_out

  # Spot cost: move[nL_ + (c-1)*nL + l] at stage t
  p_vec  <- seq_len(nCO * nL)
  c_vec  <- ((p_vec - 1L) %/% nL) + 1L
  l_vec  <- ((p_vec - 1L) %%  nL) + 1L
  spot_t <- sim$spot_coef[(c_vec - 1L) * nL * tau_sim + (l_vec - 1L) * tau_sim + t]

  obj <- c(sim$transport_coef, spot_t,       # move costs
           sim$entry_store_coef,              # entry_out holding
           sim$exit_store_coef,               # exitp_out holding
           sim$exit_short_coef)               # exitm_out shortage

  # Current-state holding costs (constants added to reported stage cost)
  state_cost <- sum(sim$entry_store_coef * entry_in) +
                sum(sim$exit_store_coef  * exitp_in) +
                sum(sim$exit_short_coef  * exitm_in)

  # Build sparse A: rows ordered as
  #   1..nCS           carrier cap (<=)
  #   nCS+1..+nCO      spot cap (<=)
  #   nCS+nCO+1..+nI   entry transition (=)
  #   ..+nJ            exit capacity (<=)
  #   ..+nJ            exit balance (=)
  nLc  <- env$nLc
  nCons <- nCS + nCO + nI + nJ + nJ
  is_ <- integer(0L);  js_ <- integer(0L);  xs_ <- numeric(0L)

  add <- function(i, j, x = 1.0) { is_ <<- c(is_, i); js_ <<- c(js_, j); xs_ <<- c(xs_, x) }

  # Carrier capacity
  for (k in seq_len(nCS))
    for (col in (nLc[k] + 1L):nLc[k + 1L]) add(k, col)
  for (c in seq_len(nCO))
    for (col in nL_ + (c - 1L) * nL + seq_len(nL)) add(nCS + c, col)

  # Entry transition: entry_out[i] + sum(move from i) = entry_in[i] + Q_t[i]
  r0 <- nCS + nCO
  for (i in seq_len(nI)) {
    add(r0 + i, col_eo[i])               # entry_out coefficient
    for (col in env$from_i[[i]]) add(r0 + i, col)
  }

  # Exit capacity: sum(move to j) <= exit_cap[j] - exitp_in[j]
  r0 <- nCS + nCO + nI
  for (j in seq_len(nJ))
    for (col in env$to_j[[j]]) add(r0 + j, col)

  # Exit balance: exitp_out[j] - exitm_out[j] - sum(move to j) = exitp_in[j] - exitm_in[j] - D_t[j]
  r0 <- nCS + nCO + nI + nJ
  for (j in seq_len(nJ)) {
    add(r0 + j, col_xp[j],  1.0)
    add(r0 + j, col_xm[j], -1.0)
    for (col in env$to_j[[j]]) add(r0 + j, col, -1.0)
  }

  A <- methods::as(
    Matrix::sparseMatrix(i = is_, j = js_, x = xs_, dims = c(nCons, nVars)),
    "CsparseMatrix"
  )

  sense <- c(rep("<=", nCS + nCO),       # carrier cap
             rep("=",  nI),              # entry transition
             rep("<=", nJ),              # exit capacity
             rep("=",  nJ))              # exit balance

  rhs <- c(env$Cb[t, seq_len(nCS)],
           env$Co[t, seq_len(nCO)],
           entry_in + Q_t,
           sim$exit_capacity - exitp_in,
           exitp_in - exitm_in - D_t)

  res <- solve_lp(list(obj = obj, A = A, rhs = rhs,
                        sense = sense, modelsense = "min", vtype = NULL))
  x <- pmax(0, res$x)

  list(cost      = state_cost + res$objval,
       entry_out = x[col_eo],
       exitp_out = x[col_xp],
       exitm_out = x[col_xm])
}

#' Compute gain of recourse: SDDP vs myopic policy
#'
#' For each instance and each OOB trial from `mstp_simulate()`, simulates a
#' myopic (greedy single-stage) policy that solves only the current-stage
#' carrier allocation LP at each period — no look-ahead — and compares its
#' cumulative cost against the SDDP policy cost on the same noise trajectory.
#'
#' Gain of recourse = (myopic_cost - SDDP_cost) / myopic_cost
#'
#' @param sims         List of instance lists (from `load_instances()`).
#' @param sddp_results List of simulation result lists from `mstp_simulate()`,
#'                   one per instance (same order as `sims`).
#' @param n_cores    Number of parallel workers (default: all but 2).
#' @return Numeric matrix of shape `(n_trials × n_instances)` — gain values.
#' @export
compute_vss <- function(sims, sddp_results, n_cores = NULL) {
  n_instances <- length(sims)
  stopifnot(length(sddp_results) == n_instances)
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 2L)

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  solve_lp  <- solve_lp
  .stage_lp <- .stage_lp

  gains <- vector("list", n_instances)

  for (idx in seq_len(n_instances)) {
    sim      <- sims[[idx]]
    env      <- build_env(sim)
    sddp_obj <- sddp_results[[idx]]$obj
    noise    <- sddp_results[[idx]]$noise
    N        <- length(sddp_obj)
    tau      <- env$tau
    taux     <- tau + 1L           # noise rows per trial (= sim$tau)

    job <- function(i) {
      entry_in <- as.numeric(sim$entry_stock_0)
      exitp_in <- as.numeric(sim$exit_stock_0)
      exitm_in <- as.numeric(sim$exit_short_0)
      myopic_cost <- 0.0

      for (t in seq_len(tau)) {
        row <- (i - 1L) * taux + t
        Q_t <- noise[row, env$I_]
        D_t <- noise[row, env$nI + env$J_]
        s   <- .stage_lp(t, entry_in, exitp_in, exitm_in, Q_t, D_t, env, sim)
        myopic_cost <- myopic_cost + s$cost
        entry_in    <- s$entry_out
        exitp_in    <- s$exitp_out
        exitm_in    <- s$exitm_out
      }
      # Terminal holding cost: inventory state entering the last period
      myopic_cost <- myopic_cost +
        sum(sim$entry_store_coef * entry_in) +
        sum(sim$exit_store_coef  * exitp_in) +
        sum(sim$exit_short_coef  * exitm_in)

      (myopic_cost - sddp_obj[i]) / myopic_cost
    }

    parallel::clusterExport(cl,
      c("sim", "env", "sddp_obj", "noise", "tau", "taux",
        "solve_lp", ".stage_lp"),
      envir = environment())
    gains[[idx]] <- parallel::parLapply(cl, seq(N), job) |> as.numeric()
  }

  do.call(cbind, gains)
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
