# M10: does n_scenarios=100 also fix the crash + validity at the 500-iter T3
# setting (6x6x20 tau=12 lambda=700) that failed under n_scenarios=10?

suppressMessages(library(MSTP))
root <- normalizePath("."); options(warn = 1)

t0 <- Sys.time()
tab <- mstp_reproduce_bounds(list(
  list(label = "6x6x20-500", tau = 12L, nOrigins = 6L, nDestinations = 6L,
       nCarriers = 20L, lambda = 700, n_scenarios = 100L, iters = 500L, trials = 200L)
), seed = 42L)
cat(sprintf("\n=== M10: 6x6x20 tau=12 500 iters, n_scenarios=100 (%.0fs) ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
print(tab[, c("label", "LB", "UB", "gap_pct", "margin_se", "valid")], row.names = FALSE, digits = 6)
cat(sprintf("trained without crash: TRUE ; valid: %s ; paper 6x6x20 gap ~ +1.2-1.5%%\n", tab$valid))
