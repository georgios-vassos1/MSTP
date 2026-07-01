#!/usr/bin/env Rscript
# Generate all instance JSONs for the paper-table refresh, using the FIXED
# generator (low-overlap, differentiated rates, deterministic init_stock=0).
# Output: rerun/inst/<name>.json   Run from package root.
suppressMessages(source("R/geninst.R"))
dir.create("rerun/inst", showWarnings = FALSE, recursive = TRUE)
SEEDS <- c(42L, 43L, 44L)
emit <- function(name, ...) {
  inst <- generate_instance(...)
  jsonlite::write_json(inst, sprintf("rerun/inst/%s.json", name), auto_unbox = TRUE, digits = 12)
}
# Table 3 (scalability, tau=12): 6x6x20, 20x20x100 (seed-averaged); 40x40x100 (seed 42 only, overnight)
for (s in SEEDS) {
  emit(sprintf("scal_6_%d",  s), tau=12L, nOrigins=6L,  nDestinations=6L,  nCarriers=20L,  seed=s, lambda=700)
  emit(sprintf("scal_20_%d", s), tau=12L, nOrigins=20L, nDestinations=20L, nCarriers=100L, seed=s, lambda=700)
}
emit("scal_40_42", tau=12L, nOrigins=40L, nDestinations=40L, nCarriers=100L, seed=42L, lambda=700)
# Table 4 (horizons, 6x6x20): tau in {12,26,52}, seed-averaged
for (s in SEEDS) for (tau in c(12L,26L,52L))
  emit(sprintf("hor_%d_%d", tau, s), tau=tau, nOrigins=6L, nDestinations=6L, nCarriers=20L, seed=s, lambda=700)
# Table 5 (sensitivity, 6x6x20 tau=12): per-lambda instances (cover lambda sweep + base for rho/m), seed-averaged
for (s in SEEDS) for (lam in c(200,700,2000))
  emit(sprintf("sens_%d_%d", lam, s), tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, seed=s, lambda=lam)
# Table 7 (capacity opt, 6x6x20 lambda=50): single instance
emit("capopt_42", tau=12L, nOrigins=6L, nDestinations=6L, nCarriers=20L, seed=42L, lambda=50)
cat("instances written to rerun/inst/\n")
