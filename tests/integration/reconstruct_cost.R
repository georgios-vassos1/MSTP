# Verify the R cost model matches the SDDP stage objective (inst/julia/model.jl:
# holding on incoming inventory + contract transport + spot cost). We recompute
# each trajectory's cost in pure R from the states and moves that mstp_simulate
# records, and require it to equal the SDDP objective. This validates the cost
# accounting + move/state layout with NO LP and NO TLPR — the foundation for a
# correct pure-R regret/VSS.

suppressMessages(library(JuliaCall))
root <- normalizePath(".")
options(warn = 1)
JuliaCall::julia_setup(installJulia = FALSE)
JuliaCall::julia_command(sprintf('import Pkg; Pkg.activate(raw"%s"); Pkg.offline(true)',
                                 file.path(root, "inst", "julia")))
JuliaCall::julia_source(file.path(root, "inst", "julia", "mstp.jl"))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
.mstp$loaded <- TRUE

# Non-square topology (nL = nI*nJ = 6 != nOD = nI+nJ = 5) so the spot-cost index
# distinguishes nL from nOD (they coincide on square 2x2 instances).
inst <- generate_instance(tau = 4L, nOrigins = 2L, nDestinations = 3L,
                          nCarriers = 3L, seed = 42L, lambda = 30.0)
S    <- mstp_corrmat(2, 3, rho_cross = 0.0)
cfg  <- mstp_config(inst, lambda = 30.0, corrmat = S, n_scenarios = 10L)
m    <- mstp_train(cfg, iterations = 30L, seed = 5L)
sim  <- mstp_simulate(m, cfg, trials = 25L, seed = 5L)

tau <- sim$tau; nI <- sim$nOrigins; nJ <- sim$nDestinations
nLanes <- sim$nLanes; nSpot <- sim$nSpotLanes
nL <- nI * nJ                    # spot lanes per spot carrier (model.jl / types.jl)
nSCo <- nSpot %/% nL             # number of spot carriers
trials <- sim$trials

recon <- numeric(trials)
for (i in seq_len(trials)) {
  total <- 0
  for (t in seq_len(tau)) {
    r <- (i - 1L) * tau + t
    hold <- sum(inst$entry_store_coef * sim$entry[r, ]) +
            sum(inst$exit_store_coef  * sim$exitp[r, ]) +
            sum(inst$exit_short_coef  * sim$exitm[r, ])
    contract <- sum(inst$transport_coef * sim$moves[r, seq_len(nLanes)])
    # spot-cost index (model.jl:55): carrier-major, lane-inner, lane count = nL
    kdx  <- as.vector(vapply(seq_len(nSCo),
              function(c) (c - 1L) * nL * tau + (0:(nL - 1L)) * tau + t,
              numeric(nL)))
    spot <- sum(inst$spot_coef[kdx] * sim$moves[r, nLanes + seq_len(nSpot)])
    total <- total + hold + contract + spot
  }
  recon[i] <- total
}

err <- max(abs(recon - sim$obj))
cat(sprintf("max |recon - SDDP obj| over %d trials = %.6g\n", trials, err))
ok <- err < 1e-4
cat(if (ok) "PASS: R cost model reconstructs SDDP objective\n"
    else    "FAIL: R cost model does NOT match SDDP objective\n")
quit(status = if (ok) 0L else 1L)
