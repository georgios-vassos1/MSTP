# Bound-validity decision rule (R). A valid SDDP lower bound for a minimisation
# satisfies bound <= E[policy cost]; we accept it when it does not exceed the
# out-of-sample mean by more than z standard errors (i.e. by more than Monte
# Carlo noise can explain). This is the single source of the rule — it used to
# live in Julia (inst/julia/api.jl) but R now owns evaluation, so the rule is
# here where it can be unit-tested without a solver.

#' Decide whether an SDDP lower bound is numerically valid
#'
#' Pure decision rule: given a lower `bound` and a sample of realised policy
#' `costs`, decide whether `bound <= E[cost]` within Monte-Carlo noise. A valid
#' bound must not exceed the sample mean by more than `z` standard errors.
#'
#' @param bound Numeric scalar. The SDDP lower bound.
#' @param costs Numeric vector (length >= 1) of realised out-of-sample policy
#'   costs.
#' @param z Numeric. Standard-error margin allowed before flagging (default 3).
#' @return A named list:
#'   \describe{
#'     \item{mean}{Sample mean of `costs`.}
#'     \item{se}{Standard error of the mean (sample sd / sqrt(n)).}
#'     \item{margin_se}{Standard errors the bound sits above the mean; positive
#'       means overshoot. `Inf`/`-Inf` when `se == 0` and the bound differs from
#'       the mean, `0` when they are equal.}
#'     \item{n}{Number of cost samples.}
#'     \item{valid}{`TRUE` iff `bound <= mean + z * se`.}
#'   }
#' @export
mstp_bound_validity <- function(bound, costs, z = 3.0) {
  stopifnot(
    is.numeric(bound), length(bound) == 1L, is.finite(bound),
    is.numeric(costs), is.numeric(z), length(z) == 1L
  )
  n <- length(costs)
  if (n == 0L) stop("costs must be non-empty", call. = FALSE)

  m  <- mean(costs)
  sd <- if (n > 1L) stats::sd(costs) else 0.0   # sample sd (n-1), matches Julia
  se <- sd / sqrt(n)
  margin_se <- if (se > 0) {
    (bound - m) / se
  } else if (bound > m) {
    Inf
  } else if (bound < m) {
    -Inf
  } else {
    0
  }
  list(mean = m, se = se, margin_se = margin_se, n = n,
       valid = bound <= m + z * se)
}
