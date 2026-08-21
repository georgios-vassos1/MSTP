# M8 Tier A: reproduce Table 6 (SDDP convergence on the scarce 6x6 instance,
# tau=12, rho_cross=0). Paper: gap closes 32.5% -> 2.2% over iters {5..300},
# UB flat from iter ~5 (draft.tex:499-504). Scarce regime = lambda=20, 95% util.
# Same instance, increasing iteration budget; seed 42.

suppressMessages(library(MSTP))
root <- normalizePath("."); options(warn = 1)

specs <- lapply(c(5L, 10L, 20L, 50L, 100L, 300L), function(it)
  list(label = sprintf("it%d", it), tau = 12L, nOrigins = 6L, nDestinations = 6L,
       nCarriers = 20L, lambda = 20, target_util = 0.95, rho_cross = 0.0,
       n_scenarios = 100L, iters = it, trials = 300L))

tab <- mstp_reproduce_bounds(specs, seed = 42L, n_insample = 2000L)
cat("\n=== M8 Table 6: SDDP convergence (scarce 6x6x20, tau=12, seed 42) ===\n")
print(tab[, c("iters", "LB", "UB", "gap_pct", "ub_insample", "gap_cert_pct", "valid")], row.names = FALSE, digits = 6)

# Emit Fig. 4 data (in-sample convergence certificate) for rerun/fig_make.py fig6().
fd <- file.path(root, "rerun", "figdata"); dir.create(fd, showWarnings = FALSE, recursive = TRUE)
.jarr <- function(x) paste0("[", paste(formatC(x, format = "f", digits = 4), collapse = ", "), "]")
conv <- sprintf(paste0('{\n "instance": "capreg_l20u95", "topology": "6x6x20",\n',
                       ' "iters": [%s],\n "lb": %s,\n "ub": %s,\n "se": %s,\n',
                       ' "rho_cross": 0.0, "n_insample": 2000, "tau": 12,\n',
                       ' "ub_kind": "in-sample statistical UB (mean + 1.96*SE)"\n}\n'),
                paste(as.integer(tab$iters), collapse = ", "),
                .jarr(tab$LB), .jarr(tab$ub_insample), .jarr(rep(0, nrow(tab))))
writeLines(conv, file.path(fd, "convergence.json"))
cat(sprintf("wrote %s (Fig. 4 data)\n", file.path(fd, "convergence.json")))
cat("\nPaper (draft.tex:499-504) gap%: 5->32.5, 10->19.2, 20->10.3, 50->4.3, 100->2.9, 300->2.2\n")
cat(sprintf("CHECK LB rises monotonically with iters: %s\n", all(diff(tab$LB) > -1e-6)))
cat(sprintf("CHECK gap closes to < 5%% by 300 iters: %s (%.2f%%)\n",
            tail(tab$gap_pct, 1) < 5, tail(tab$gap_pct, 1)))
