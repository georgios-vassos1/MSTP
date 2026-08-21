# Smoke test for the in-sample-UB extension of mstp_reproduce_bounds.
suppressMessages(library(MSTP))
root <- normalizePath("."); options(warn = 1)

tab <- mstp_reproduce_bounds(list(
  list(label = "smoke", tau = 4L, nOrigins = 2L, nDestinations = 2L, nCarriers = 4L,
       lambda = 5.5, rho_cross = 0.4, iters = 200L, trials = 200L, n_scenarios = 100L)
), seed = 42L, n_insample = 1000L)
cat("\n=== smoke: with n_insample ===\n"); print(tab, row.names = FALSE)
cat(sprintf("has ub_insample: %s | has gap_cert_pct: %s | LB<=ub_insample: %s\n",
            "ub_insample" %in% names(tab), "gap_cert_pct" %in% names(tab),
            all(tab$LB <= tab$ub_insample)))

# opt-out must be identical to the old 9-column output
tab0 <- mstp_reproduce_bounds(list(
  list(label = "smoke0", tau = 4L, nOrigins = 2L, nDestinations = 2L, nCarriers = 4L,
       lambda = 5.5, rho_cross = 0.4, iters = 50L, trials = 50L, n_scenarios = 50L)
), seed = 42L)
cat(sprintf("opt-out: %d cols, has ub_insample: %s (expect 9, FALSE)\n",
            ncol(tab0), "ub_insample" %in% names(tab0)))
