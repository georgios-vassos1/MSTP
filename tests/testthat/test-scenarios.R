# Tests for the deterministic scenario contract (R/scenarios.R).
# Statistical assertions use large samples with generous tolerances so they are
# not flaky; determinism assertions are exact.

test_that("mstp_corrmat has correct block structure and is symmetric PD", {
  S <- mstp_corrmat(2, 3, rho_oo = 0.6, rho_dd = 0.7, rho_cross = 0.1)
  expect_equal(dim(S), c(5L, 5L))
  expect_true(isSymmetric(S))
  expect_equal(diag(S), rep(1, 5))
  # origin block (1:2) off-diagonal
  expect_equal(S[1, 2], 0.6)
  # destination block (3:5) off-diagonal
  expect_equal(S[3, 4], 0.7)
  expect_equal(S[4, 5], 0.7)
  # cross entries
  expect_equal(S[1, 3], 0.1)
  expect_equal(S[2, 5], 0.1)
  # positive definite
  expect_true(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) > 0)
})

test_that("mstp_corrmat rejects a non-positive-definite specification", {
  # rho_cross close to 1 with unit within-block correlations breaks PD
  expect_error(mstp_corrmat(3, 3, rho_oo = 0.0, rho_dd = 0.0, rho_cross = 0.99),
               "positive definite")
})

test_that("mstp_sample_copula is deterministic under a fixed seed", {
  S <- mstp_corrmat(2, 2)
  a <- mstp_sample_copula(50, lambda = 30, corrmat = S, seed = 42)
  b <- mstp_sample_copula(50, lambda = 30, corrmat = S, seed = 42)
  c <- mstp_sample_copula(50, lambda = 30, corrmat = S, seed = 43)
  expect_identical(a, b)
  expect_false(identical(a, c))
})

test_that("mstp_sample_copula does not leak global RNG state", {
  S <- mstp_corrmat(2, 2)
  set.seed(123)
  before <- runif(1)
  set.seed(123)
  invisible(mstp_sample_copula(10, lambda = 5, corrmat = S, seed = 999))
  after <- runif(1)
  expect_identical(before, after)  # seeded call must not advance caller's stream
})

test_that("mstp_sample_copula returns correct shape, type, and support", {
  S <- mstp_corrmat(2, 3)
  x <- mstp_sample_copula(100, lambda = 10, corrmat = S, seed = 1)
  expect_equal(dim(x), c(100L, 5L))
  expect_true(is.integer(x))
  expect_true(all(x >= 0L))
})

test_that("marginals match Poisson(lambda) in mean and variance", {
  S <- mstp_corrmat(2, 2, rho_oo = 0.6, rho_dd = 0.7)
  lam <- c(50, 50, 80, 80)
  x <- mstp_sample_copula(50000, lambda = lam, corrmat = S, seed = 7)
  m <- colMeans(x)
  v <- apply(x, 2, var)
  expect_equal(unname(m), lam, tolerance = 0.03)   # ~3% relative
  expect_equal(unname(v), lam, tolerance = 0.06)   # Poisson var == mean
})

test_that("copula reproduces target correlation sign and magnitude", {
  # Large lambda: Poisson ~ Normal, so Pearson corr tracks the copula corr well.
  S <- mstp_corrmat(2, 2, rho_oo = 0.6, rho_dd = 0.7, rho_cross = 0.0)
  x <- mstp_sample_copula(60000, lambda = 700, corrmat = S, seed = 11)
  R <- cor(x)
  expect_equal(R[1, 2], 0.6, tolerance = 0.03)   # origin-origin
  expect_equal(R[3, 4], 0.7, tolerance = 0.03)   # dest-dest
  expect_lt(abs(R[1, 3]), 0.03)                  # cross-block ~ 0
  expect_lt(abs(R[2, 4]), 0.03)
})

test_that("positive cross-block correlation is recovered", {
  S <- mstp_corrmat(2, 2, rho_oo = 0.6, rho_dd = 0.7, rho_cross = 0.4)
  x <- mstp_sample_copula(60000, lambda = 700, corrmat = S, seed = 13)
  R <- cor(x)
  expect_equal(R[1, 3], 0.4, tolerance = 0.03)
  expect_equal(R[2, 4], 0.4, tolerance = 0.03)
})

test_that("mstp_scenario_support yields tau independent per-stage supports", {
  S <- mstp_corrmat(2, 2)
  Om <- mstp_scenario_support(tau = 12, n_scenarios = 10, lambda = 40,
                              corrmat = S, seed = 5)
  expect_length(Om, 12L)
  expect_true(all(vapply(Om, function(m) all(dim(m) == c(10L, 4L)), logical(1))))
  # stages are drawn independently: first two supports should differ
  expect_false(identical(Om[[1]], Om[[2]]))
  # whole structure reproducible under the same seed
  Om2 <- mstp_scenario_support(tau = 12, n_scenarios = 10, lambda = 40,
                               corrmat = S, seed = 5)
  expect_identical(Om, Om2)
})

test_that("mstp_scenario_paths yields n_paths trajectories of shape tau x d", {
  S <- mstp_corrmat(2, 3)
  P <- mstp_scenario_paths(n_paths = 200, tau = 12, lambda = 40,
                           corrmat = S, seed = 9)
  expect_length(P, 200L)
  expect_true(all(vapply(P, function(m) all(dim(m) == c(12L, 5L)), logical(1))))
  expect_identical(P, mstp_scenario_paths(200, 12, 40, S, seed = 9))
})
