# Integration proof: the refactored engine is deterministic across the full
# R -> JuliaCall -> SDDP.jl round trip. Same seed => identical LB and identical
# out-of-sample costs; different seed => different. This is the acceptance
# criterion for M2. Run from the package root:
#   Rscript tests/integration/test_engine_determinism.R
#
# It is NOT part of the fast unit suite (it starts Julia and trains SDDP).

suppressMessages(library(JuliaCall))

root <- normalizePath(".")
options(warn = 1)

# --- wire the engine without installing the package -------------------------
JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf(
  'import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
  file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE   # engine already sourced above; skip setup_engine()

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "PASS: " else "FAIL: ", msg, "\n", sep=""); ok <<- ok && cond }

# --- reshape convention: R flatten order must match Julia reshape -----------
# Determinism alone would not catch a scrambled (stage,scenario,dim) mapping, so
# verify the flat arrays round-trip through the Julia builders unchanged.
Om_chk <- list(matrix(c(1,2,3, 4,5,6),    nrow = 2, byrow = TRUE),   # stage 1: 2 scen x 3 dim
               matrix(c(7,8,9, 10,11,12), nrow = 2, byrow = TRUE))   # stage 2
sup <- JuliaCall::julia_call("build_support", .flatten_support(Om_chk),
                             2L, 2L, 3L, need_return = "R")
check(identical(as.numeric(sup[[1]][[1]]), c(1,2,3)) &&
      identical(as.numeric(sup[[2]][[2]]), c(10,11,12)),
      "support reshape preserves (stage, scenario, dim) order")

P_chk <- list(matrix(c(1,2, 3,4), nrow = 2, byrow = TRUE),           # path 1: tau 2 x dim 2
              matrix(c(5,6, 7,8), nrow = 2, byrow = TRUE))           # path 2
hist <- JuliaCall::julia_call("build_historical", .flatten_paths(P_chk),
                              2L, 2L, 2L, need_return = "R")
# hist[[path]][[stage]] = (stage_index, noise_vector)
check(identical(as.numeric(hist[[1]][[1]][[2]]), c(1,2)) &&
      identical(as.numeric(hist[[2]][[2]][[2]]), c(7,8)),
      "trajectory reshape preserves (path, stage, dim) order")

# --- tiny instance ----------------------------------------------------------
inst <- generate_instance(tau = 4L, nOrigins = 2L, nDestinations = 2L,
                          nCarriers = 3L, seed = 42L, lambda = 30.0)
S   <- mstp_corrmat(2, 2, rho_oo = 0.6, rho_dd = 0.7, rho_cross = 0.0)
cfg <- mstp_config(inst, lambda = 30.0, corrmat = S, n_scenarios = 5L)

ITERS <- 20L; TRIALS <- 40L

run <- function(seed) {
  m  <- mstp_train(cfg, iterations = ITERS, seed = seed)
  lb <- mstp_bound(m)
  ob <- mstp_simulate(m, cfg, trials = TRIALS, seed = seed)$obj
  # Bound validity is the SE-tolerant rule (LB <= E[cost] within z SE), not a
  # raw LB <= mean comparison: on tiny instances LB can sit slightly above the
  # OOB mean within Monte-Carlo noise (a documented small-instance artefact).
  vb <- mstp_validate_bound(m, cfg, trials = 200L, z = 3.0, seed = seed)
  list(lb = lb, ub = mean(ob), obj = ob, valid = vb$valid, margin = vb$margin_se)
}

a  <- run(100L)   # baseline
b  <- run(100L)   # same seed  -> must be identical
c  <- run(200L)   # diff seed  -> must differ

cat(sprintf("seed100 run1: LB=%.6f  UB=%.6f  valid=%s (%.2f SE)\n", a$lb, a$ub, a$valid, a$margin))
cat(sprintf("seed100 run2: LB=%.6f  UB=%.6f  valid=%s (%.2f SE)\n", b$lb, b$ub, b$valid, b$margin))
cat(sprintf("seed200     : LB=%.6f  UB=%.6f  valid=%s (%.2f SE)\n", c$lb, c$ub, c$valid, c$margin))

check(identical(a$lb, b$lb),               "LB identical across same-seed runs")
check(identical(a$obj, b$obj),             "OOB cost vector identical across same-seed runs")
check(isTRUE(a$valid),                     "bound valid within 3 SE (validate_bound)")
check(!isTRUE(all.equal(a$lb, c$lb)),      "LB differs under a different seed")

cat(if (ok) ">>> M2 DETERMINISM PROOF PASSED\n" else ">>> M2 DETERMINISM PROOF FAILED\n")
quit(status = if (ok) 0L else 1L)
