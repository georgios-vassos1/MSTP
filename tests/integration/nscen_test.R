# M10 quick test: does enlarging the fixed backward support (n_scenarios) remove
# the LB overshoot seen at n_scenarios=10 on 6x6x20 tau=12? If margin_se drops
# below 3 (valid) as n_scenarios grows, the fixed-small-support SAA bias is
# confirmed as the cause and the fix direction (larger / resampled support) holds.

suppressMessages(library(JuliaCall))
root <- normalizePath("."); options(warn = 1)
JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf('import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
                                 file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE

specs <- lapply(c(10L, 50L, 100L), function(ns)
  list(label = sprintf("nscen%d", ns), tau = 12L, nOrigins = 6L, nDestinations = 6L,
       nCarriers = 20L, lambda = 700, n_scenarios = ns, iters = 100L, trials = 200L))

tab <- mstp_reproduce_bounds(specs, seed = 42L)
cat("\n=== M10: LB validity vs backward support size (6x6x20 tau=12 lambda=700, 100 iters) ===\n")
print(tab[, c("label", "LB", "UB", "gap_pct", "margin_se", "valid")], row.names = FALSE, digits = 6)
cat(sprintf("\nHYPOTHESIS (larger support -> valid): margin_se falls with n_scenarios: %s\n",
            all(diff(tab$margin_se) < 0)))
cat(sprintf("n_scenarios=100 valid: %s\n", tab$valid[tab$label == "nscen100"]))
