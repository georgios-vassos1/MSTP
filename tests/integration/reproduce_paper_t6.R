# M8 Tier A: reproduce Table 6 (SDDP convergence on the scarce 6x6 instance,
# tau=12, rho_cross=0). Paper: gap closes 32.5% -> 2.2% over iters {5..300},
# UB flat from iter ~5 (draft.tex:499-504). Scarce regime = lambda=20, 95% util.
# Same instance, increasing iteration budget; seed 42.

suppressMessages(library(JuliaCall))
root <- normalizePath("."); options(warn = 1)
JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf('import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
                                 file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE

specs <- lapply(c(5L, 10L, 20L, 50L, 100L, 300L), function(it)
  list(label = sprintf("it%d", it), tau = 12L, nOrigins = 6L, nDestinations = 6L,
       nCarriers = 20L, lambda = 20, target_util = 0.95, rho_cross = 0.0,
       n_scenarios = 100L, iters = it, trials = 300L))

tab <- mstp_reproduce_bounds(specs, seed = 42L)
cat("\n=== M8 Table 6: SDDP convergence (scarce 6x6x20, tau=12, seed 42) ===\n")
print(tab[, c("iters", "LB", "UB", "gap_pct", "valid")], row.names = FALSE, digits = 6)
cat("\nPaper (draft.tex:499-504) gap%: 5->32.5, 10->19.2, 20->10.3, 50->4.3, 100->2.9, 300->2.2\n")
cat(sprintf("CHECK LB rises monotonically with iters: %s\n", all(diff(tab$LB) > -1e-6)))
cat(sprintf("CHECK gap closes to < 5%% by 300 iters: %s (%.2f%%)\n",
            tail(tab$gap_pct, 1) < 5, tail(tab$gap_pct, 1)))
