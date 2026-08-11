# M12: Figure 5 data -- SDDP-dual projected-gradient capacity trajectory on the
# scarce 6x6x20 instance. Reduced scale (fewer outer/cold iters) to keep it fast;
# the shape is what Figure 5 shows: objective trends down, gradient norm bounded
# and noisy, capacity step shrinks (~alpha0/sqrt(k)) (draft.tex:530).

suppressMessages(library(JuliaCall))
root <- normalizePath("."); options(warn = 1)
JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf('import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
                                 file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE

inst <- generate_instance(tau = 12L, nOrigins = 6L, nDestinations = 6L,
                          nCarriers = 20L, seed = 42L, lambda = 20, target_util = 0.95)
S   <- mstp_corrmat(6, 6, rho_cross = 0.0)
cfg <- mstp_config(inst, lambda = 20, corrmat = S, n_scenarios = 20L)
v   <- 0.2 * mean(inst$transport_coef)

tr <- mstp_capacity_trajectory(cfg, v = v, outer = 8L, cold = 25L,
                               eval_iters = 25L, n_eval = 100L, n_capdual = 50L, seed = 42L)

cat("\n=== M12: Figure 5 cap-opt trajectory (scarce 6x6x20) ===\n")
print(data.frame(iter = seq_along(tr$objective),
                 objective = round(tr$objective, 1),
                 grad_norm = round(tr$grad_norm, 3),
                 step_norm = round(tr$step_norm, 3)), row.names = FALSE)
cat(sprintf("\nCHECK objective net-decreases (last < first): %s\n", tail(tr$objective,1) < tr$objective[1]))
cat(sprintf("CHECK gradient norm finite/bounded: %s\n", all(is.finite(tr$grad_norm))))
cat(sprintf("CHECK step shrinks (last < first): %s\n", tail(tr$step_norm,1) < tr$step_norm[1]))
