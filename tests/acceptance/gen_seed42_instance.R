#!/usr/bin/env Rscript
# Pure R (no Julia): generate the exact seed-42 6x6x20 / tau=52 / lambda=700
# instance with the 1e6 storage caps used in the production tau=52 run, and
# write it to JSON for the pure-Julia acceptance check. See README.md.
#
# Usage:  Rscript tests/acceptance/gen_seed42_instance.R [out.json]
args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args) >= 1) args[[1]] else "/tmp/inst52_seed42.json"

self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(self)) normalizePath(file.path(dirname(self), "..", "..")) else "."
source(file.path(root, "R", "geninst.R"))

inst <- generate_instance(tau = 52L, nOrigins = 6L, nDestinations = 6L,
                          nCarriers = 20L, seed = 42L, lambda = 700,
                          entry_capacity = rep(1000000L, 6L),
                          exit_capacity  = rep(1000000L, 6L))
jsonlite::write_json(inst, out, auto_unbox = TRUE, digits = 12)
cat("wrote", out,
    "| carrier_capacity:", length(inst$carrier_capacity),
    "transport_coef:", length(inst$transport_coef), "\n")
