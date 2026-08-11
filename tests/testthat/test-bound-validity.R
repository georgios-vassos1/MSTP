# Tests for the bound-validity decision rule (R/bound_validity.R). Ported from
# the former Julia rule test (tests/julia/test_bound_validity.jl); pure R, no
# solver, deterministic. The end-to-end SDDP validity check is exercised by
# tests/integration/test_engine_determinism.R.

test_that("bound well below the cost sample is valid and sits below the mean", {
  r <- mstp_bound_validity(90.0, rep(100.0, 50L), z = 3.0)
  expect_true(r$valid)
  expect_lt(r$margin_se, 0)
  expect_equal(r$mean, 100.0)
})

test_that("bound far above the mean with negligible noise is invalid", {
  costs <- 100.0 + 0.5 * sin(seq_len(200L))   # mean ~ 100, tiny spread
  r <- mstp_bound_validity(130.0, costs, z = 3.0)
  expect_false(r$valid)
  expect_gt(r$margin_se, 3.0)
})

test_that("bound above the mean but within Monte-Carlo noise stays valid", {
  noisy <- 100.0 + 20.0 * sin(seq_len(100L))  # SE ~ 1.4
  r <- mstp_bound_validity(101.0, noisy, z = 3.0)
  expect_true(r$valid)
  expect_gt(r$margin_se, 0)
  expect_lt(r$margin_se, 3)
})

test_that("zero-variance sample: above mean invalid, at/below mean valid", {
  expect_false(mstp_bound_validity(101.0, rep(100.0, 10L))$valid)
  expect_true(mstp_bound_validity(100.0, rep(100.0, 10L))$valid)
  expect_true(mstp_bound_validity(99.0,  rep(100.0, 10L))$valid)
  # margin_se conventions when se == 0
  expect_identical(mstp_bound_validity(101.0, rep(100.0, 10L))$margin_se, Inf)
  expect_identical(mstp_bound_validity(99.0,  rep(100.0, 10L))$margin_se, -Inf)
  expect_identical(mstp_bound_validity(100.0, rep(100.0, 10L))$margin_se, 0)
})

test_that("degenerate input is rejected, not silently mis-judged", {
  expect_error(mstp_bound_validity(1.0, numeric(0L)), "non-empty")
})

test_that("mean/se match base R moments", {
  set.seed(1)
  costs <- rnorm(500L, 100, 10)
  r <- mstp_bound_validity(95, costs, z = 3)
  expect_equal(r$mean, mean(costs))
  expect_equal(r$se, sd(costs) / sqrt(length(costs)))
})
