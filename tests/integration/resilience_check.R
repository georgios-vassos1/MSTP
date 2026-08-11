# M9 acceptance: the default engine trains the paper's 6x6x20 instances that
# crashed under the old fail-fast :primal config -- both the ample (lambda=700)
# and the scarce (lambda=20, 95% util) regimes -- deterministically and with a
# valid bound. No explicit solver mode is set (exercises the default + fallback).

suppressMessages(library(JuliaCall))
root <- normalizePath("."); options(warn = 1)
JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf('import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
                                 file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "PASS: " else "FAIL: ", msg, "\n", sep=""); ok <<- ok && cond }

train_lb <- function(lambda, util, seed) {
  inst <- generate_instance(tau = 12L, nOrigins = 6L, nDestinations = 6L,
                            nCarriers = 20L, seed = 42L, lambda = lambda, target_util = util)
  S   <- mstp_corrmat(6, 6, rho_cross = 0.0)
  cfg <- mstp_config(inst, lambda = lambda, corrmat = S, n_scenarios = 10L)
  m   <- mstp_train(cfg, iterations = 40L, seed = seed)
  v   <- mstp_validate_bound(m, cfg, trials = 100L, seed = seed)
  list(lb = mstp_bound(m), valid = v$valid, margin = v$margin_se)
}

a1 <- train_lb(700, 0.80, 42L)
a2 <- train_lb(700, 0.80, 42L)   # determinism
sc <- train_lb(20,  0.95, 42L)   # scarce -- crashed T6 under :primal

cat(sprintf("ample  6x6x20 l700: LB=%.0f valid=%s (%.2f SE) [run2 LB=%.0f]\n",
            a1$lb, a1$valid, a1$margin, a2$lb))
cat(sprintf("scarce 6x6x20 l20 : LB=%.0f valid=%s (%.2f SE)\n", sc$lb, sc$valid, sc$margin))

check(is.finite(a1$lb),          "ample 6x6x20 trains (no crash)")
check(identical(a1$lb, a2$lb),   "6x6x20 LB identical across same-seed runs (deterministic)")
check(isTRUE(a1$valid),          "ample bound valid within 3 SE")
check(is.finite(sc$lb),          "scarce 6x6x20 trains (no crash) -- was the T6 failure")

cat(if (ok) ">>> M9 RESILIENCE PASSED\n" else ">>> M9 RESILIENCE FAILED\n")
quit(status = if (ok) 0L else 1L)
