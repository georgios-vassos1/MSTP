# M12: regret-based stopping table -- regret vs SDDP iteration budget on the 2x2
# instance (tau=4, lambda=5.5, rho=0.4, 500 OOB). Paper's point (G7): regret is
# small and stable across budgets, so regret-based early stopping is viable.

suppressMessages(library(JuliaCall))
root <- normalizePath("."); options(warn = 1)
JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf('import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
                                 file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE

rows <- lapply(c(100L, 250L, 500L, 1000L, 1500L), function(it)
  mstp_reproduce_recourse(list(label = sprintf("it%d", it), tau = 4L, nOrigins = 2L,
    nDestinations = 2L, nCarriers = 4L, lambda = 5.5, rho_cross = 0.4,
    n_scenarios = 10L, iters = it, trials = 500L), seed = 42L))
tab <- do.call(rbind, rows)
cat("\n=== M12: regret-based stopping (2x2, 500 OOB) ===\n")
print(tab[, c("label", "mean_regret_pct", "median_regret_pct", "q95_regret_pct")],
      row.names = FALSE, digits = 4)
cat("\nPaper (rerun/RESULTS.md:39-43) mean%: .046/.145/.135/.043/.057 ; q95%: .18/.30/.21/.17/.26\n")
cat(sprintf("CHECK regret small (all mean < 1%%) and stable across budgets: %s\n",
            all(tab$mean_regret_pct < 1)))
