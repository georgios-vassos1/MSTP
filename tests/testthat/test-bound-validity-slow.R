# SLOW acceptance gate for the SDDP lower-bound overshoot (reviewer concern #3).
#
# This is the regression test that actually captures the bug: on the failing
# regime (6x6x20, long horizon, lambda=700) a *valid* SDDP lower bound must not
# exceed the simulated mean by more than a few standard errors.
#
# STATUS: expected RED until the numerics are genuinely fixed. The storage-cap
# change in this branch reduces the overshoot margin directionally on small
# instances (1e6 -> physical caps: margin 1.27 -> 0.27 SE at 2x2/tau=52) but has
# NOT been shown to close the gap at 6x6x20 scale, where the real run logs show
# LB - UB ~= 7 standard errors (tau=52). Do not mark this green until a re-run
# confirms validity here.
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
  corrmat <- mstp_gen_corrmat(2L, 6L, cross_corr = 0.0)
  config  <- mstp_config(inst, lambda = rep(700, 12), corrmat = corrmat,
                         n_scenarios = 10L)
  model   <- mstp_train(config, iterations = 300L)
  v       <- mstp_validate_bound(model, config, trials = 200L, z = 3.0)

  message(sprintf("bound=%.0f sim_mean=%.0f se=%.0f margin_se=%.2f valid=%s",
                  v$bound, v$sim_mean, v$sim_se, v$margin_se, v$valid))
  # Acceptance criterion for the fix:
  expect_true(v$valid)
})
