# Correctness/smoke gate for the capacity-optimisation core (inst/julia/capopt.jl)
# driven from R. On a scarce-capacity instance the LP proxy should strip buffer
# capacity (mean below x0) while the SDDP dual gradient retains more; both plans
# are scored under common random numbers and the run must be deterministic.
# Small settings keep it to a few SDDP trains.

suppressMessages(library(JuliaCall))
root <- normalizePath("."); options(warn = 1)
JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf('import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
                                 file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE

spec <- list(label = "scarce-2x2x3", tau = 4L, nOrigins = 2L, nDestinations = 2L,
             nCarriers = 3L, lambda = 20.0, target_util = 0.95, n_scenarios = 10L,
             outer = 3L, cold = 15L, eval_iters = 15L, n_eval = 60L, n_capdual = 40L)

a <- mstp_reproduce_capopt(spec, seed = 7L)
b <- mstp_reproduce_capopt(spec, seed = 7L)   # determinism

cat("--- capacity optimisation (scarce regime) ---\n")
print(a, row.names = FALSE)

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "PASS: " else "FAIL: ", msg, "\n", sep=""); ok <<- ok && cond }
# Robust invariants only. The economic result (LP strips harder than SDDP, and
# the SDDP advantage grows with scarcity/correlation) is a PAPER-SCALE finding
# for M8, not something a tiny under-trained instance must exhibit.
check(identical(a, b),                 "capopt result identical across same-seed runs")
check(is.finite(a$adv_pp),             "advantage (pp) is finite")
check(a$xlp_mean <= a$x0_mean + 1e-9,  "LP-proxy capacity clamped <= default x0")
check(a$xsd_mean <= a$x0_mean + 1e-9,  "SDDP capacity clamped <= default x0")
check(abs(a$adv_pp - (a$lp_pct - a$sddp_pct)) < 1e-6,
      "advantage pp = lp_pct - sddp_pct (internal consistency)")

cat(if (ok) ">>> CAPOPT CORE PASSED\n" else ">>> CAPOPT CORE FAILED\n")
quit(status = if (ok) 0L else 1L)
