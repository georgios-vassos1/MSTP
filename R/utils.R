# Solve a JuMP/Gurobi-compatible model list with HiGHS.
# Accepts model lists with fields: obj, A, rhs, sense, modelsense, vtype
# Returns list with: objval, objbound, x, status
solve_lp <- function(model) {
  n   <- ncol(model$A)
  sns <- model$sense

  is_le <- sns %in% c("<", "<=")
  is_ge <- sns %in% c(">", ">=")
  lhs <- ifelse(sns == "=", model$rhs, ifelse(is_le, -Inf, model$rhs))
  rhs <- ifelse(sns == "=", model$rhs, ifelse(is_le,  model$rhs, Inf))

  types <- if (!is.null(model$vtype)) model$vtype else rep("C", n)

  result <- highs::highs_solve(
    L       = model$obj,
    lower   = rep(0.0, n),
    upper   = rep(Inf, n),
    A       = model$A,
    lhs     = lhs,
    rhs     = rhs,
    types   = types,
    maximum = identical(model$modelsense, "max"),
    control = highs::highs_control(log_to_console = FALSE)
  )

  list(
    objval   = result$objective_value,
    objbound = result$objective_value,
    x        = result$primal_solution,
    status   = result$status
  )
}


#' Convert an MSTP instance to a TLPR-format JSON file
#'
#' Produces a JSON file that \code{TLPR::loadProblemDataCx} (and the R-side
#' \code{TLPR::dp_config}) can consume, enabling exact DP and RTDP to be run
#' on instances originally designed for SDDP.
#'
#' Because the TLPR C++ backend shares a single \code{flowSupport} vector for
#' both supply (Q) and demand (D), you must supply a single distribution object
#' \code{Q} whose \code{vals} are used for both inflow and outflow scenarios.
#' Pass a matching \code{D} for the R-side \code{dp_config} setup (defaults to
#' \code{Q} if omitted).
#'
#' \strong{Sizing note:} MSTP instances are typically generated with large
#' \code{exit_capacity} (e.g. 10,000) and carrier capacities in the hundreds.
#' For TLPR feasibility with a small \code{R}, regenerate the MSTP instance
#' with \code{exit_capacity = R} and \code{carrier_capacity} scaled to
#' \code{O(R)}, or pass a pre-scaled \code{inst}.
#'
#' @param inst  MSTP instance list (from \code{MSTP::generate_instance()}).
#' @param R     Integer storage limit for TLPR (\eqn{R = \max(\text{inventory})}).
#' @param Q     List with \code{vals} (integer supply/demand support, length
#'   \eqn{n_Q}) and \code{prob} (probability weights).
#' @param W     List with \code{vals} (numeric spot-rate support, length
#'   \eqn{n_W}) and \code{prob}.  If \code{NULL}, derived from
#'   \code{quantile(inst$spot_coef)} using \code{nW} quantiles.
#' @param path  File path for the output JSON.
#' @param nCO   Number of spot carriers to expose to TLPR (default:
#'   \code{inst$nSpotCarriers}).
#' @param nW    Number of quantile-based W support points when \code{W = NULL}.
#'   Ignored if \code{W} is supplied.
#' @param D     List with \code{vals} and \code{prob} for the demand
#'   distribution seen by \code{dp_config}.  Defaults to \code{Q}.
#' @return The assembled JSON list (invisibly).
#' @export
mstp_to_tlpr_json <- function(inst, R, Q, W = NULL, path, nCO = NULL, nW = 3L, D = NULL) {
  nI   <- inst$nOrigins
  nJ   <- inst$nDestinations
  nCS  <- inst$nCarriers
  nSC  <- if (!is.null(nCO)) as.integer(nCO) else {
    if (!is.null(inst$nSpotCarriers)) inst$nSpotCarriers else inst$nCarriers
  }
  nCO  <- nSC
  nL   <- nI * nJ
  tau  <- inst$tau
  nL_  <- length(inst$Ldx)

  # ── Normalise Q; default D to Q ──────────────────────────────────────────
  Q$prob <- Q$prob / sum(Q$prob)
  if (is.null(D)) {
    D <- Q
  } else {
    D$prob <- D$prob / sum(D$prob)
  }
  nQ <- length(Q$vals)

  # ── Derive W from spot_coef quantiles when not supplied ───────────────────
  if (is.null(W)) {
    p    <- seq(0, 1, length.out = nW + 2L)[-c(1L, nW + 2L)]
    W    <- list(
      vals = as.numeric(round(quantile(inst$spot_coef, probs = p, names = FALSE), 2)),
      prob = rep(1.0 / nW, nW)
    )
  } else {
    W$prob <- W$prob / sum(W$prob)
  }
  nW_ <- length(W$vals)

  # ── Auction structure ─────────────────────────────────────────────────────
  # winner: named list, carrier key → integer vector of bid indices
  winner    <- setNames(lapply(inst$Winners, as.integer), as.character(seq(nCS)))
  winnerKey <- as.character(seq(nCS))

  # carrierIdx: importList<int> expects {"k": [single_int]}
  # In R: named list of scalars → toJSON wraps each as [k]
  carrierIdx <- setNames(as.list(seq_len(nCS)), as.character(seq(nCS)))

  # CTb_list: per-carrier contract rates split by inst$nLc
  CTb_list <- setNames(
    lapply(seq(nCS), function(k) {
      as.numeric(inst$transport_coef[(inst$nLc[k] + 1L):inst$nLc[k + 1L]])
    }),
    as.character(seq(nCS))
  )

  # ── Routing indices ───────────────────────────────────────────────────────
  # Combined variable index: [contracted_lanes (Ldx), spot_lanes (rep(1:nL, nCO))]
  combined <- c(inst$Ldx, rep(seq(nL), nCO))

  fromx <- outer(seq(nI), (seq(nJ) - 1L) * nI, "+")  # nI × nJ
  tox   <- outer((seq(nJ) - 1L) * nI, seq(nI), "+")  # nJ × nI

  from_i <- apply(fromx, 1L, function(ldx) as.integer(which(combined %in% ldx)), simplify = FALSE)
  to_j   <- apply(tox,   1L, function(ldx) as.integer(which(combined %in% ldx)), simplify = FALSE)

  # ── Lanes: nL × 2 matrix, rows = [origin, destination] ──────────────────
  # Ordering: for j=1..nJ, i=1..nI → same as TLPR CartesianProductX(I_, J_)
  L_mat <- cbind(I = rep(seq(nI), times = nJ), J = rep(seq(nJ), each = nI))

  # ── Capacities ────────────────────────────────────────────────────────────
  Cb <- matrix(as.integer(inst$carrier_capacity[seq(nCS * tau)]),          nrow = tau)
  Co_start <- nCS * tau + 1L
  Co <- matrix(as.integer(inst$carrier_capacity[Co_start:(Co_start + nSC * tau - 1L)]), nrow = tau)

  # ── Spot rates: tau × (nCO * nL) ─────────────────────────────────────────
  CTo <- matrix(as.numeric(inst$spot_coef), nrow = tau)

  # ── State keys (mixed-radix for state indexing) ───────────────────────────
  nSI       <- R + 1L
  nSJ       <- 2L * R + 1L
  si_powers <- nSI ^ seq(nI)
  stateKeys <- if (nJ == 1L) {
    c(1L, si_powers)
  } else {
    c(1L, si_powers, (nSI ^ nI) * nSJ ^ seq(nJ - 1L))
  }

  # ── Flow keys (mixed-radix for scenario indexing) ─────────────────────────
  nQdx     <- nQ ^ nI
  nDdx     <- nQ ^ nJ   # C++ symmetry: same nQ for demand
  flowKeys <- as.integer(c(
    nQ ^ (seq(nI) - 1L),
    nQdx * nQ ^ (seq(nJ) - 1L),
    nQdx * nDdx * nW_ ^ (seq(nCO) - 1L)
  ))

  # ── Holding / shortage cost vector ───────────────────────────────────────
  alpha <- as.numeric(c(inst$entry_store_coef, inst$exit_store_coef, inst$exit_short_coef))

  # ── Assemble ──────────────────────────────────────────────────────────────
  obj <- list(
    # Scalar dimensions (C++ reads as [value][0])
    R    = as.integer(R),
    nQ   = as.integer(nQ),
    nW   = as.integer(nW_),
    nI   = as.integer(nI),
    nJ   = as.integer(nJ),
    nCS  = as.integer(nCS),
    nCO  = as.integer(nCO),
    nL_  = as.integer(nL_),
    nL   = as.integer(nL),
    # DP horizon and cost (R-side fields)
    tau   = as.integer(tau),
    alpha = alpha,
    # Distributions (R-side dp_config)
    Q = list(vals = as.integer(Q$vals), prob = as.numeric(Q$prob)),
    D = list(vals = as.integer(D$vals), prob = as.numeric(D$prob)),
    W = list(vals = as.numeric(W$vals), prob = as.numeric(W$prob)),
    # Network topology (C++)
    from_i    = from_i,
    to_j      = to_j,
    B         = lapply(inst$Bids,    as.integer),
    L         = lapply(seq(nrow(L_mat)), function(k) as.integer(L_mat[k, ])),
    winnerKey = winnerKey,
    # Capacities (C++)
    Cb  = lapply(seq(tau), function(t) as.integer(Cb[t, ])),
    Co  = lapply(seq(tau), function(t) as.integer(Co[t, ])),
    CTo = lapply(seq(tau), function(t) as.numeric(CTo[t, ])),
    # Cipher keys (C++)
    stateKeys = as.integer(stateKeys),
    flowKeys  = flowKeys,
    # Carrier data (C++)
    winner    = winner,
    CTb_list  = CTb_list,
    carrierIdx = carrierIdx,
    # Auxiliary R-side fields
    nLc = as.integer(inst$nLc)
  )

  writeLines(jsonlite::toJSON(obj, pretty = TRUE), path)
  message("TLPR JSON written to: ", path)
  invisible(obj)
}

# The former TLPR-based clairvoyant LP builder (build_model) and its env adapter
# (build_env) were removed: the multi-period transportation LP is now built in
# pure R by .clairvoyant_cost() in R/lp_transport.R (verified against the SDDP
# stage model), so the package has no TLPR dependency.
