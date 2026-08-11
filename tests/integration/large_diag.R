# M11 diagnostic: which HiGHS mode trains the 20x20x100 instance that crashes
# (node 5, early iteration) under dual_presolve and ipm? Small n_scenarios/iters
# -- the crash is an early forward-pass primal failure, independent of support
# size -- so this is a fast probe of the solver, not a full run.

suppressMessages(library(JuliaCall))
root <- normalizePath("."); options(warn = 1)
JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf('import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
                                 file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE

inst <- generate_instance(tau = 12L, nOrigins = 20L, nDestinations = 20L,
                          nCarriers = 100L, seed = 42L, lambda = 700)
S   <- mstp_corrmat(20, 20, rho_cross = 0.0)
cfg <- mstp_config(inst, lambda = 700, corrmat = S, n_scenarios = 20L)

for (mode in c("dual_presolve", "ipm", "pdlp")) {
  JuliaCall::julia_call("set_opt_mode!", mode)
  t0 <- Sys.time()
  msg <- tryCatch({
    m <- mstp_train(cfg, iterations = 15L, seed = 42L)  # primary = `mode` via set_opt_mode!
    sprintf("OK  LB=%.0f  (%.0fs)", mstp_bound(m),
            as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }, error = function(e) paste0("CRASH: ", substr(conditionMessage(e), 1, 60)))
  cat(sprintf("%-14s -> %s\n", mode, msg))
}
JuliaCall::julia_call("set_opt_mode!", "dual_presolve")
