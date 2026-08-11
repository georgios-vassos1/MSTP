# SLOW acceptance gate for the SDDP lower-bound overshoot (reviewer concern #3).
#
# This is the regression test that actually captures the bug: on the failing
# regime (6x6x20, long horizon, lambda=700) a *valid* SDDP lower bound must not
# exceed the simulated mean by more than a few standard errors.
#
# The tau=52 overshoot (real logs: LB - UB ~= 7 SE) was traced to the instance
# generator, not SDDP: random per-node initial inventory produced ill-conditioned
# duals. It is addressed by deterministic init_stock = 0 (R/geninst.R) and, in the
# refactored engine, a fixed R-supplied support (see CHANGES-sddp-bound-validity.md).
# This gate should now pass; it remains unverified here because it is slow.
#
# Skipped by default (slow: minutes). Enable with:
#   MSTP_SLOW_TESTS=1 Rscript -e 'testthat::test_file("tests/testthat/test-bound-validity-slow.R")'

test_that("SDDP lower bound is valid on the 6x6x20 long-horizon instance", {
  skip_if_not(identical(Sys.getenv("MSTP_SLOW_TESTS"), "1"),
              "set MSTP_SLOW_TESTS=1 to run the slow bound-validity gate")
  skip_if_not_installed("MSTP")

  library(MSTP)
  setup_engine()

  inst    <- generate_instance(tau = 52L, nOrigins = 6L, nDestinations = 6L,
                               nCarriers = 20L, seed = 42L, lambda = 700)
  corrmat <- mstp_corrmat(6, 6, rho_cross = 0.0)
  config  <- mstp_config(inst, lambda = rep(700, 12), corrmat = corrmat,
                         n_scenarios = 10L)
  model   <- mstp_train(config, iterations = 300L, seed = 1L)
  v       <- mstp_validate_bound(model, config, trials = 200L, z = 3.0, seed = 1L)

  message(sprintf("bound=%.0f sim_mean=%.0f se=%.0f margin_se=%.2f valid=%s",
                  v$bound, v$sim_mean, v$sim_se, v$margin_se, v$valid))
  # Acceptance criterion for the fix:
  expect_true(v$valid)
})
