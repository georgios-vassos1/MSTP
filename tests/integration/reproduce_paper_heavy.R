# M8 heavy tables (run in background): T3 scalability, T4 horizons, T5
# sensitivity, T7/T8 capacity optimisation. Each experiment is wrapped so a
# failure (or the process being killed) leaves earlier results on disk.
# Results are appended to tests/integration/paper_results.txt.

suppressMessages(library(JuliaCall))
root <- normalizePath("."); options(warn = 1)
OUT <- file.path(root, "tests", "integration", "paper_results.txt")
say <- function(...) { cat(..., "\n"); cat(..., "\n", file = OUT, append = TRUE) }

JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf('import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
                                 file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE

run <- function(tag, expr) {   # expr is a lazy promise -- forced exactly once
  say(sprintf("\n### %s ###", tag))
  res <- tryCatch(expr, error = function(e) { say("ERROR:", conditionMessage(e)); NULL })
  if (!is.null(res)) {
    cap <- utils::capture.output(print(res, row.names = FALSE))
    cat(cap, sep = "\n"); cat("\n")
    cat(cap, sep = "\n", file = OUT, append = TRUE); cat("\n", file = OUT, append = TRUE)
  }
  invisible(res)
}

say(sprintf("=== M8 heavy reproduction (seed 42) ==="))

# T2 regret + gain of recourse (paper small topologies: tau=4, lambda=5.5,
# rho=0.4, 1000 iters, 500 OOB). Gain is the CORRECT model.jl-accounted value
# (~14%), not the paper's double-counted 52.3%.
run("T2 regret + gain of recourse", do.call(rbind, lapply(list(
  list(label="2x1", nOrigins=2L, nDestinations=1L, nCarriers=3L),
  list(label="1x2", nOrigins=1L, nDestinations=2L, nCarriers=3L),
  list(label="2x2", nOrigins=2L, nDestinations=2L, nCarriers=4L)
), function(t) mstp_reproduce_recourse(c(t, list(tau=4L, lambda=5.5, rho_cross=0.4,
  iters=1000L, trials=500L, n_scenarios=10L)), seed=42L))))

# n_scenarios=100 gives valid bounds (M10); 50 for the huge 40x40 to bound compute.
# T3 scalability (lambda=700, seed 42) -- paper Table 3
run("T3 scalability", mstp_reproduce_bounds(list(
  list(label="6x6x20",    tau=12L, nOrigins=6L,  nDestinations=6L,  nCarriers=20L,  lambda=700, iters=500L, trials=50L, n_scenarios=100L, n_insample=2000L),
  list(label="20x20x100", tau=12L, nOrigins=20L, nDestinations=20L, nCarriers=100L, lambda=700, iters=50L,  trials=50L, n_scenarios=100L, n_insample=500L),
  list(label="40x40x100", tau=12L, nOrigins=40L, nDestinations=40L, nCarriers=100L, lambda=700, iters=20L,  trials=50L, n_scenarios=50L,  n_insample=100L)
), seed = 42L))

# T4 horizons (6x6x20, lambda=700, 300 iters) -- paper Table 4
run("T4 horizons", mstp_reproduce_bounds(list(
  list(label="tau12", tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=700, iters=300L, trials=200L, n_scenarios=100L, n_insample=2000L),
  list(label="tau26", tau=26L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=700, iters=300L, trials=200L, n_scenarios=100L, n_insample=2000L),
  list(label="tau52", tau=52L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=700, iters=300L, trials=200L, n_scenarios=100L, n_insample=2000L)
), seed = 42L))

# T5 sensitivity (6x6x20 tau=12, 300 iters; single seed 42 -- paper averages 42-44)
run("T5 sensitivity", mstp_reproduce_bounds(list(
  list(label="lam200",  tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=200,  iters=300L, trials=200L, n_scenarios=100L),
  list(label="lam700",  tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=700,  iters=300L, trials=200L, n_scenarios=100L),
  list(label="lam2000", tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=2000, iters=300L, trials=200L, n_scenarios=100L),
  list(label="rho0.0",  tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=700, rho_cross=0.0, iters=300L, trials=200L, n_scenarios=100L),
  list(label="rho0.2",  tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=700, rho_cross=0.2, iters=300L, trials=200L, n_scenarios=100L),
  list(label="rho0.4",  tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=700, rho_cross=0.4, iters=300L, trials=200L, n_scenarios=100L),
  list(label="m1",      tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=700, spot_mult=1, iters=300L, trials=200L, n_scenarios=100L),
  list(label="m2",      tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=700, spot_mult=2, iters=300L, trials=200L, n_scenarios=100L),
  list(label="m4",      tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, lambda=700, spot_mult=4, iters=300L, trials=200L, n_scenarios=100L)
), seed = 42L))

# T7 capacity optimisation regimes (6x6x20, paper settings) -- paper Table 7
run("T7 ample (lam50 util0.80)",
    mstp_reproduce_capopt(list(label="ample", tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L,
                               lambda=50, target_util=0.80, outer=15L), seed=42L))
run("T7 scarce (lam20 util0.95)",
    mstp_reproduce_capopt(list(label="scarce", tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L,
                               lambda=20, target_util=0.95, outer=15L), seed=42L))

# T8 correlation sweep at scarce capacity -- paper Table 8
for (rc in c(0.0, 0.2, 0.4, 0.6)) {
  run(sprintf("T8 scarce rho_cross=%.1f", rc),
      mstp_reproduce_capopt(list(label=sprintf("rho%.1f", rc), tau=12L, nOrigins=6L, nDestinations=6L,
                                 nCarriers=20L, lambda=20, target_util=0.95, rho_cross=rc, outer=15L), seed=42L))
}

say("=== DONE M8 heavy ===")
