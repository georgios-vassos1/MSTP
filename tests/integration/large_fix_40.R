# M11: does the state bound also let 40x40x100 (the largest instance) train?
suppressMessages(library(MSTP))
root <- normalizePath("."); options(warn = 1)

inst <- generate_instance(tau = 12L, nOrigins = 40L, nDestinations = 40L,
                          nCarriers = 100L, seed = 42L, lambda = 700)
S   <- mstp_corrmat(40, 40, rho_cross = 0.0)
cfg <- mstp_config(inst, lambda = 700, corrmat = S, n_scenarios = 20L)
t0 <- Sys.time()
msg <- tryCatch({
  m <- mstp_train(cfg, iterations = 20L, seed = 42L)
  v <- mstp_validate_bound(m, cfg, trials = 50L, seed = 42L)
  sprintf("OK  LB=%.0f UB=%.0f valid=%s  (%.0fs)", v$bound, v$sim_mean, v$valid,
          as.numeric(difftime(Sys.time(), t0, units = "secs")))
}, error = function(e) paste0("CRASH: ", substr(conditionMessage(e), 1, 70)))
cat("40x40x100 (state-bounded) ->", msg, "\n")
