# Reproduce Table 5 (sensitivity, 6x6x20 tau=12, 300 iters, 200 OOB) on the
# deterministic R-first engine, seed-averaged over 42-44 (paper protocol; port of
# rerun/tables_final.jl:119-123). Unlike the legacy pipeline -- whose sweep seed
# only picked the instance file while training/simulation drew from an
# uncontrolled global RNG (so identical configs diverged) -- mstp_reproduce_bounds
# seeds every call, so lambda=700 == rho=0.0 == m=1 by construction.
#   Cost = out-of-sample mean; sd = within-seed OOS sd (= sim_se*sqrt(trials));
#   both averaged over the three seeds.
# Writes rerun/figdata/sensitivity.json for Fig. 3 (fig_make.py fig2).
suppressMessages(library(JuliaCall))
root <- normalizePath("."); options(warn = 1)
JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf('import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
                                 file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE

SEEDS <- c(42L, 43L, 44L)
base  <- list(tau = 12L, nOrigins = 6L, nDestinations = 6L, nCarriers = 20L,
              iters = 300L, trials = 200L, n_scenarios = 100L)

# One sensitivity config -> seed-averaged (cost mean, within-seed OOS sd).
sens_row <- function(label, over) {
  spec <- modifyList(base, c(list(label = label), over))
  rs <- lapply(SEEDS, function(s) mstp_reproduce_bounds(list(spec), seed = s))
  ub <- vapply(rs, function(r) r$UB, numeric(1))                          # OOS mean
  sd <- vapply(rs, function(r) r$sim_se * sqrt(spec$trials), numeric(1))  # within-seed OOS sd
  cat(sprintf("  [%s] cost=%.0f sd=%.0f  (per-seed cost: %s)\n",
              label, mean(ub), mean(sd), paste(round(ub), collapse = "/"))); flush.console()
  data.frame(label = label, cost_mean = mean(ub), cost_sd = mean(sd),
             stringsAsFactors = FALSE)
}

cfgs <- list(
  list("lam200",  list(lambda = 200,  rho_cross = 0.0)),
  list("lam700",  list(lambda = 700,  rho_cross = 0.0)),
  list("lam2000", list(lambda = 2000, rho_cross = 0.0)),
  list("rho0.0",  list(lambda = 700,  rho_cross = 0.0)),
  list("rho0.2",  list(lambda = 700,  rho_cross = 0.2)),
  list("rho0.4",  list(lambda = 700,  rho_cross = 0.4)),
  list("m1",      list(lambda = 700,  rho_cross = 0.0, spot_mult = 1)),
  list("m2",      list(lambda = 700,  rho_cross = 0.0, spot_mult = 2)),
  list("m4",      list(lambda = 700,  rho_cross = 0.0, spot_mult = 4))
)
cat("\n=== M8 Table 5: sensitivity (6x6x20 tau=12, 300 iters, 200 OOB, seed-avg 42-44, n_sc=100) ===\n")
tab <- do.call(rbind, lapply(cfgs, function(x) sens_row(x[[1]], x[[2]])))
print(tab, row.names = FALSE, digits = 7)

# Determinism check: baseline (lambda=700 / rho=0.0 / m=1) must be identical.
b <- tab$cost_mean[tab$label %in% c("lam700", "rho0.0", "m1")]
cat(sprintf("CHECK baseline identical (lam700=rho0.0=m1): %s (range %.6f)\n",
            diff(range(b)) < 1e-6, diff(range(b))))

# Emit Fig. 3 data for rerun/fig_make.py fig2().
fd <- file.path(root, "rerun", "figdata"); dir.create(fd, showWarnings = FALSE, recursive = TRUE)
g    <- function(lab) tab[tab$label == lab, ]
trip <- function(x, lab) { r <- g(lab); sprintf("[%s, %d, %d]", x, round(r$cost_mean), round(r$cost_sd)) }
sens <- sprintf(paste0('{\n "lam": [%s, %s, %s],\n "rho": [%s, %s, %s],\n',
                       ' "mm": [%s, %s, %s]\n}\n'),
  trip("200", "lam200"), trip("700", "lam700"), trip("2000", "lam2000"),
  trip("0.0", "rho0.0"), trip("0.2", "rho0.2"), trip("0.4", "rho0.4"),
  trip("1", "m1"), trip("2", "m2"), trip("4", "m4"))
writeLines(sens, file.path(fd, "sensitivity.json"))
cat(sprintf("wrote %s (Fig. 3 data)\n", file.path(fd, "sensitivity.json")))
