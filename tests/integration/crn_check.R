# CRN-robustness check for draft.tex:549/553 -- "each advantage reproduces within
# ~1-2 pp on an independent common-random-number set". Re-runs the scarce-regime
# rho-sweep capacity optimisation at an INDEPENDENT seed (43) with the same config
# as the committed seed-42 Table 8 run (outer=15, defaults cold=100/eval=300/
# n_eval=500/n_capdual=200) and compares the SDDP advantage (pp) against seed 42.
suppressMessages(library(MSTP))
root <- normalizePath("."); options(warn = 1)

# Committed seed-42 advantages (draft Table 8 / phase4_capopt.out).
seed42 <- c("0.0" = 23.648, "0.2" = 30.757, "0.4" = 40.337, "0.6" = 61.033)

tab <- do.call(rbind, lapply(c(0.0, 0.2, 0.4, 0.6), function(rc)
  mstp_reproduce_capopt(list(label = sprintf("rho%.1f", rc), tau = 12L, nOrigins = 6L,
                             nDestinations = 6L, nCarriers = 20L, lambda = 20,
                             target_util = 0.95, rho_cross = rc, outer = 15L), seed = 43L)))

res <- data.frame(rho = c(0.0, 0.2, 0.4, 0.6),
                  seed42 = as.numeric(seed42),
                  seed43 = round(tab$adv_pp, 2),
                  abs_diff = round(abs(tab$adv_pp - as.numeric(seed42)), 2))
cat("\n=== CRN check: SDDP advantage (pp), seed 43 (independent) vs seed 42 ===\n")
print(res, row.names = FALSE)
cat(sprintf("max |diff| = %.2f pp | claim 'within ~1-2 pp' holds: %s\n",
            max(res$abs_diff), max(res$abs_diff) <= 2.0))
