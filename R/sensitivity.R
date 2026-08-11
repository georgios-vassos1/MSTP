# Sensitivity analyses: the distribution of perfect-foresight (clairvoyant) cost
# as inflow volumes or spot rates vary. Both use the pure-R clairvoyant LP
# (R/lp_transport.R) — no TLPR. These are descriptive cost-distribution tools;
# the paper's Table-5 parameter sensitivity of the SDDP bound is produced instead
# by mstp_reproduce_bounds().

#' Inflow sensitivity of the clairvoyant cost
#'
#' For each instance, draws `n_samples` inflow trajectories from `flow_support`
#' (outflow held at `demand`) and returns the perfect-foresight LP cost of each.
#'
#' @param sims List of instances (from [generate_instance()]/[load_instances()]).
#' @param flow_support Integer vector — the per-node inflow sample space.
#' @param demand Constant per-node outflow while inflow varies. Default:
#'   `round(mean(flow_support))`.
#' @param n_samples Number of inflow draws per instance.
#' @param seed Optional integer seed (scoped; global RNG is preserved).
#' @param n_cores Parallel workers (default 1 = serial).
#' @return Numeric matrix `(n_samples x n_instances)` of clairvoyant costs.
#' @export
sensitivity_inflow <- function(sims, flow_support = seq(20L, 40L, by = 5L),
                               demand = NULL, n_samples = 200L,
                               seed = NULL, n_cores = 1L) {
  n_instances <- length(sims)
  if (is.null(demand)) demand <- round(mean(flow_support))
  cols <- vector("list", n_instances)

  for (idx in seq_len(n_instances)) {
    inst <- sims[[idx]]
    lay  <- .transport_layout(inst)
    init <- .init_state(inst)
    tau <- inst$tau; nI <- inst$nOrigins; nJ <- inst$nDestinations
    Dmat <- matrix(demand, nrow = tau, ncol = nJ)

    Qdraws <- .with_seed(seed, {
      lapply(seq_len(n_samples),
             function(k) matrix(sample(flow_support, tau * nI, replace = TRUE),
                                nrow = tau, ncol = nI))
    })
    f <- function(k) .clairvoyant_cost(lay, Qdraws[[k]], Dmat, init)
    cols[[idx]] <- .apply_maybe_parallel(
      seq_len(n_samples), f, n_cores,
      export = c(".clairvoyant_cost", ".spot_obj", "solve_lp"))
  }
  matrix(unlist(cols), nrow = n_samples)
}

#' Spot-rate sensitivity of the clairvoyant cost
#'
#' For each instance, holds inflow and outflow at `demand` and draws `n_samples`
#' spot-rate multipliers uniformly from `mult_range`, scaling the spot-cost
#' coefficients, and returns the perfect-foresight LP cost of each.
#'
#' @param sims List of instances.
#' @param demand Constant per-node inflow and outflow (numeric scalar).
#' @param mult_range Length-2 numeric — range of the spot-rate multiplier.
#' @param n_samples Number of draws per instance.
#' @param seed Optional integer seed (scoped).
#' @param n_cores Parallel workers (default 1 = serial).
#' @return Numeric matrix `(n_samples x n_instances)` of clairvoyant costs.
#' @export
sensitivity_spot <- function(sims, demand, mult_range = c(1, 4),
                             n_samples = 200L, seed = NULL, n_cores = 1L) {
  stopifnot(is.numeric(demand), length(demand) == 1L,
            is.numeric(mult_range), length(mult_range) == 2L)
  n_instances <- length(sims)
  cols <- vector("list", n_instances)

  for (idx in seq_len(n_instances)) {
    inst <- sims[[idx]]
    lay0 <- .transport_layout(inst)
    init <- .init_state(inst)
    tau <- inst$tau; nI <- inst$nOrigins; nJ <- inst$nDestinations
    Qmat <- matrix(demand, nrow = tau, ncol = nI)
    Dmat <- matrix(demand, nrow = tau, ncol = nJ)

    mults <- .with_seed(seed, stats::runif(n_samples, mult_range[1L], mult_range[2L]))
    base_spot <- lay0$spot_coef
    f <- function(k) {
      lay <- lay0; lay$spot_coef <- base_spot * mults[k]
      .clairvoyant_cost(lay, Qmat, Dmat, init)
    }
    cols[[idx]] <- .apply_maybe_parallel(
      seq_len(n_samples), f, n_cores,
      export = c(".clairvoyant_cost", ".spot_obj", "solve_lp"))
  }
  matrix(unlist(cols), nrow = n_samples)
}

#' Plot sensitivity density curves
#'
#' @param optx Matrix from [sensitivity_inflow()] or [sensitivity_spot()].
#' @param title Plot title.
#' @param scale_x Divisor for x-axis labels (default 1000).
#' @export
plot_sensitivity <- function(optx, title = "Sensitivity", scale_x = 1000) {
  densities <- apply(optx, 2L, stats::density)
  ylim <- range(sapply(densities, function(d) d$y))
  xlim <- range(sapply(densities, function(d) d$x))
  graphics::plot(NA, xlim = xlim, ylim = ylim, xlab = "Total cost", ylab = "Density",
       main = title, xaxt = "n", yaxt = "n")
  for (d in densities)
    graphics::lines(d, col = grDevices::rgb(0, 0, 0, alpha = 0.2), lwd = 2L)
  rngx <- seq(xlim[1L], xlim[2L], length.out = 6L)
  graphics::axis(1L, at = rngx, labels = round(rngx / scale_x, 1L))
}
