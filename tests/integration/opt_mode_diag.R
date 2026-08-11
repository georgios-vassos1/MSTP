# M9 diagnostic: which HiGHS configuration trains the 6x6x20 (lambda=700)
# instance that crashes under the default :primal mode? Tries each mode for a
# short run and reports success + bound (or the crash).

suppressMessages(library(JuliaCall))
root <- normalizePath("."); options(warn = 1)
JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf('import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
                                 file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE

inst <- generate_instance(tau = 12L, nOrigins = 6L, nDestinations = 6L,
                          nCarriers = 20L, seed = 42L, lambda = 700)
S   <- mstp_corrmat(6, 6, rho_cross = 0.0)
cfg <- mstp_config(inst, lambda = 700, corrmat = S, n_scenarios = 10L)

for (mode in c("primal", "presolve", "dual_presolve", "ipm")) {
  JuliaCall::julia_call("set_opt_mode!", mode)
  t0 <- Sys.time()
  msg <- tryCatch({
    m <- mstp_train(cfg, iterations = 30L, seed = 42L)
    v <- mstp_validate_bound(m, cfg, trials = 50L, seed = 42L)
    sprintf("OK  LB=%.0f UB=%.0f valid=%s  (%.0fs)",
            v$bound, v$sim_mean, v$valid, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }, error = function(e) paste0("CRASH: ", substr(conditionMessage(e), 1, 70)))
  cat(sprintf("%-14s -> %s\n", mode, msg))
}
JuliaCall::julia_call("set_opt_mode!", "primal")  # restore default
