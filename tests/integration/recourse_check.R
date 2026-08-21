# Correctness gate for the pure-R transport LPs (R/lp_transport.R). For each
# out-of-sample trajectory: the SDDP policy's realised cost is a FEASIBLE
# solution of the clairvoyant LP, so a correct clairvoyant optimum cannot exceed
# it (regret >= 0); and the foresight optimum cannot exceed the myopic policy
# (clairvoyant <= myopic). Both must hold on every trajectory.

suppressMessages(library(MSTP))
root <- normalizePath(".")
options(warn = 1)

inst <- generate_instance(tau = 4L, nOrigins = 2L, nDestinations = 3L,
                          nCarriers = 3L, seed = 42L, lambda = 30.0)
S    <- mstp_corrmat(2, 3, rho_cross = 0.0)
cfg  <- mstp_config(inst, lambda = 30.0, corrmat = S, n_scenarios = 10L)
m    <- mstp_train(cfg, iterations = 60L, seed = 5L)
sim  <- mstp_simulate(m, cfg, trials = 30L, seed = 5L)

lay  <- .transport_layout(inst)
tau  <- sim$tau; nI <- sim$nOrigins; nJ <- sim$nDestinations
init <- list(entry = as.numeric(inst$entry_stock_0),
             exitp = as.numeric(inst$exit_stock_0),
             exitm = as.numeric(inst$exit_short_0))

N <- sim$trials
clair <- myop <- numeric(N)
for (i in seq_len(N)) {
  rows <- (i - 1L) * tau + seq_len(tau)
  Qmat <- sim$noise[rows, seq_len(nI), drop = FALSE]
  Dmat <- sim$noise[rows, nI + seq_len(nJ), drop = FALSE]
  clair[i] <- .clairvoyant_cost(lay, Qmat, Dmat, init)
  myop[i]  <- .myopic_cost(lay, Qmat, Dmat, init)
}
sddp <- sim$obj
regret <- (sddp - clair) / clair
gain   <- (myop - sddp) / myop

cat(sprintf("mean regret = %.4f%%   (min %.4f%%)   mean gain-of-recourse = %.2f%%\n",
            100 * mean(regret), 100 * min(regret), 100 * mean(gain)))

ok <- TRUE
check <- function(cond, msg) { cat(if (cond) "PASS: " else "FAIL: ", msg, "\n", sep=""); ok <<- ok && cond }
# Relative tolerance: clairvoyant is the mathematical lower bound; any tiny
# excess over SDDP/myopic is LP solver feasibility tolerance (~1e-7) on costs ~1e4.
tol <- function(x) 1e-6 * pmax(1, abs(x))
check(all(clair <= sddp + tol(sddp)),  "clairvoyant <= SDDP cost on every trajectory (regret >= 0)")
check(all(clair <= myop  + tol(myop)), "clairvoyant <= myopic on every trajectory")
check(all(is.finite(clair)) && all(is.finite(myop)), "all LPs solved finitely")

# --- exported API (TLPR-free) end-to-end + regression lock -----------------
sr  <- list(sim)
reg <- compute_regret(list(inst), sr)
gn  <- compute_vss(list(inst), sr)
si  <- sensitivity_inflow(list(inst), flow_support = seq(20L, 40L, by = 5L),
                          n_samples = 12L, seed = 1L)
sp  <- sensitivity_spot(list(inst), demand = 30, n_samples = 12L, seed = 1L)

check(isTRUE(all.equal(as.numeric(reg), regret, tolerance = 1e-8)),
      "compute_regret matches the direct clairvoyant computation")
check(all(reg >= -1e-6) && all(is.finite(gn)), "compute_regret >= 0 and compute_vss finite")
check(all(dim(si) == c(12L, 1L)) && all(is.finite(si)) &&
      all(dim(sp) == c(12L, 1L)) && all(is.finite(sp)), "sensitivity returns finite matrices")

cat(sprintf("mean(compute_regret)=%.8f  mean(compute_vss)=%.6f\n", mean(reg), mean(gn)))
LOCK_GAIN <- 0.059147
check(isTRUE(all.equal(mean(gn), LOCK_GAIN, tolerance = 1e-4)),
      sprintf("gain-of-recourse baseline locked (%.6f)", mean(gn)))

cat(if (ok) ">>> RECOURSE LP CORRECTNESS PASSED\n" else ">>> RECOURSE LP CORRECTNESS FAILED\n")
quit(status = if (ok) 0L else 1L)
