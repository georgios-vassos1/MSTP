# M8 Tier A: reproduce Table 2 (perfect-foresight regret) and the gain of
# recourse on the paper's small topologies, at paper settings
# (tau=4, lambda=5.5, rho_cross=0.4, 1000 SDDP iters, 10 scenarios, 500 OOB,
# seed 42 -- draft.tex:362,368). Prints actual numbers next to the paper's.
# Numbers are a new seeded baseline, not digit-matches of the old unseeded runs.

suppressMessages(library(MSTP))
root <- normalizePath("."); options(warn = 1)

base <- list(tau = 4L, lambda = 5.5, rho_cross = 0.4, iters = 1000L,
             trials = 500L, n_scenarios = 10L)
topos <- list(
  c(label = "2x1", nOrigins = 2L, nDestinations = 1L, nCarriers = 3L),
  c(label = "1x2", nOrigins = 1L, nDestinations = 2L, nCarriers = 3L),
  c(label = "2x2", nOrigins = 2L, nDestinations = 2L, nCarriers = 4L)
)

rows <- lapply(topos, function(t) {
  spec <- c(base, list(label = t[["label"]],
                       nOrigins = as.integer(t[["nOrigins"]]),
                       nDestinations = as.integer(t[["nDestinations"]]),
                       nCarriers = as.integer(t[["nCarriers"]])))
  mstp_reproduce_recourse(spec, seed = 42L)
})
tab <- do.call(rbind, rows)
cat("\n=== M8 Tier A: Table 2 regret + gain of recourse (paper settings, seed 42) ===\n")
print(tab, row.names = FALSE, digits = 4)

cat("\nPaper (draft.tex:375-377 regret %, :479 gain):\n")
cat("  2x1  regret mean/med/q95 = 0.000/0.000/0.000\n")
cat("  1x2  regret mean/med/q95 = 0.000/0.000/0.000\n")
cat("  2x2  regret mean/med/q95 = 0.055/0.000/0.209 ; gain of recourse 52.3%\n")

reg_ok  <- all(tab$q95_regret_pct < 0.5)          # paper: all q95 < 0.25%
gain_ok <- tab$mean_gain_pct[tab$label == "2x2"] > 0
cat(sprintf("\nCHECK regret q95 < 0.5%% (near-zero, paper-consistent): %s\n", all(reg_ok)))
cat(sprintf("CHECK 2x2 gain of recourse positive: %s (value %.1f%%)\n",
            gain_ok, tab$mean_gain_pct[tab$label == "2x2"]))
quit(status = if (all(reg_ok) && gain_ok) 0L else 1L)
