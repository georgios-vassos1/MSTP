# M11: does the finite state bound (model.jl) let 20x20x100 train under the
# default solver, where it previously crashed under all solver modes?

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

t0 <- Sys.time()
msg <- tryCatch({
  m <- mstp_train(cfg, iterations = 20L, seed = 42L)
  v <- mstp_validate_bound(m, cfg, trials = 50L, seed = 42L)
  sprintf("OK  LB=%.0f UB=%.0f valid=%s  (%.0fs)", v$bound, v$sim_mean, v$valid,
          as.numeric(difftime(Sys.time(), t0, units = "secs")))
}, error = function(e) paste0("CRASH: ", substr(conditionMessage(e), 1, 70)))
cat("20x20x100 (state-bounded) ->", msg, "\n")
