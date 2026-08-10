# Pure-R transportation LPs used by the recourse/sensitivity analyses. These
# mirror the SDDP stage model (inst/julia/model.jl) exactly — verified by cost
# reconstruction (tests/integration/reconstruct_cost.R) and by regret >= 0
# (tests/integration/recourse_check.R) — and replace the former TLPR-based
# `build_model`, so the package has no external LP-builder dependency.
#
# Cost accounting (per stage t, model.jl:59-65): holding on the INCOMING
# inventory + contract transport + spot cost. Lanes are numbered
# lane(i,j) = (j-1)*nI + i; a spot carrier owns nL = nI*nJ lanes; the spot-cost
# coefficient index is carrier-major, lane-inner (model.jl:55).

# Derive all move-column groupings, capacities, and cost coefficients from an
# instance (no TLPR, no Julia). Returns a plain list ("layout").
.transport_layout <- function(inst) {
  nI  <- inst$nOrigins; nJ <- inst$nDestinations; tau <- inst$tau
  nSC <- inst$nCarriers
  nSCo <- if (!is.null(inst$nSpotCarriers)) inst$nSpotCarriers else nSC
  nL  <- nI * nJ
  Ldx <- inst$Ldx; nL_ <- length(Ldx)
  nLc <- inst$nLc                                  # length nSC+1, cumulative contract cols
  nMoves <- nL_ + nSCo * nL

  comb_lane <- c(Ldx, rep(seq_len(nL), nSCo))      # lane id of each move column
  origin_lanes <- lapply(seq_len(nI), function(i) (0:(nJ - 1L)) * nI + i)
  dest_lanes   <- lapply(seq_len(nJ), function(j) (j - 1L) * nI + seq_len(nI))

  list(
    nI = nI, nJ = nJ, tau = tau, nSC = nSC, nSCo = nSCo, nL = nL, nL_ = nL_,
    nMoves = nMoves,
    from_i = lapply(origin_lanes, function(ol) which(comb_lane %in% ol)),
    to_j   = lapply(dest_lanes,   function(dl) which(comb_lane %in% dl)),
    contract_cols = lapply(seq_len(nSC),  function(k) (nLc[k] + 1L):nLc[k + 1L]),
    spot_cols     = lapply(seq_len(nSCo), function(c) nL_ + (c - 1L) * nL + seq_len(nL)),
    Cb = matrix(inst$carrier_capacity[seq_len(nSC * tau)], nrow = tau),
    Co = matrix(inst$carrier_capacity[(nSC * tau + 1L):((nSC + nSCo) * tau)], nrow = tau),
    transport = inst$transport_coef, spot_coef = inst$spot_coef,
    entry_hold = inst$entry_store_coef, exit_hold = inst$exit_store_coef,
    short = inst$exit_short_coef, exit_cap = inst$exit_capacity
  )
}

# Spot-cost coefficients for stage t, ordered to match the spot move columns.
.spot_obj <- function(lay, t) {
  kdx <- as.vector(vapply(seq_len(lay$nSCo),
           function(c) (c - 1L) * lay$nL * lay$tau + (0:(lay$nL - 1L)) * lay$tau + t,
           numeric(lay$nL)))
  lay$spot_coef[kdx]
}

# Model.jl-accounted cost of one stage given incoming inventory + moves.
.stage_cost <- function(lay, t, ein, xpin, xmin, move) {
  sum(lay$entry_hold * ein) + sum(lay$exit_hold * xpin) + sum(lay$short * xmin) +
    sum(lay$transport * move[seq_len(lay$nL_)]) +
    sum(.spot_obj(lay, t) * move[lay$nL_ + seq_len(lay$nMoves - lay$nL_)])
}

# Clairvoyant (perfect-foresight) optimal cost for a single realised trajectory.
# Qmat is tau x nI inflows, Dmat is tau x nJ outflows; init is a list with
# entry/exitp/exitm. Builds and solves the full multi-period LP.
.clairvoyant_cost <- function(lay, Qmat, Dmat, init) {
  nI <- lay$nI; nJ <- lay$nJ; tau <- lay$tau; nM <- lay$nMoves
  nSt <- nI + 2L * nJ
  nVars <- tau * nM + tau * nSt
  moff <- function(t) (t - 1L) * nM
  soff <- function(t) tau * nM + (t - 1L) * nSt
  ein  <- function(t) soff(t) + seq_len(nI)
  xpin <- function(t) soff(t) + nI + seq_len(nJ)
  xmin <- function(t) soff(t) + nI + nJ + seq_len(nJ)

  obj <- numeric(nVars)
  for (t in seq_len(tau)) {
    obj[moff(t) + seq_len(lay$nL_)] <- lay$transport
    obj[moff(t) + lay$nL_ + seq_len(nM - lay$nL_)] <- .spot_obj(lay, t)
    obj[ein(t)]  <- lay$entry_hold
    obj[xpin(t)] <- lay$exit_hold
    obj[xmin(t)] <- lay$short
  }

  I <- integer(0); J <- integer(0); X <- numeric(0); rhs <- numeric(0); sns <- character(0)
  r <- 0L
  add <- function(rows, cols, x) { I <<- c(I, rows); J <<- c(J, cols); X <<- c(X, x) }
  row <- function(cols, x, sense, b) {
    r <<- r + 1L
    add(rep(r, length(cols)), cols, rep_len(x, length(cols)))
    rhs <<- c(rhs, b); sns <<- c(sns, sense)
  }

  for (t in seq_len(tau)) {
    for (k in seq_len(lay$nSC))  row(moff(t) + lay$contract_cols[[k]], 1, "<=", lay$Cb[t, k])
    for (c in seq_len(lay$nSCo)) row(moff(t) + lay$spot_cols[[c]],     1, "<=", lay$Co[t, c])
    for (i in seq_len(nI))       row(c(moff(t) + lay$from_i[[i]], ein(t)[i]),
                                     c(rep(1, length(lay$from_i[[i]])), -1), "<=", Qmat[t, i])
    for (j in seq_len(nJ))       row(c(moff(t) + lay$to_j[[j]], xpin(t)[j]),
                                     c(rep(1, length(lay$to_j[[j]])), 1), "<=", lay$exit_cap[j])
  }
  # Fix stage-1 incoming state to the initial inventory.
  for (i in seq_len(nI)) row(ein(1L)[i],  1, "=", init$entry[i])
  for (j in seq_len(nJ)) row(xpin(1L)[j], 1, "=", init$exitp[j])
  for (j in seq_len(nJ)) row(xmin(1L)[j], 1, "=", init$exitm[j])
  # Transitions define the next incoming state.
  for (t in seq_len(tau - 1L)) {
    for (i in seq_len(nI))
      row(c(ein(t + 1L)[i], ein(t)[i], moff(t) + lay$from_i[[i]]),
          c(1, -1, rep(1, length(lay$from_i[[i]]))), "=", Qmat[t, i])
    for (j in seq_len(nJ))
      row(c(xpin(t + 1L)[j], xmin(t + 1L)[j], xpin(t)[j], xmin(t)[j], moff(t) + lay$to_j[[j]]),
          c(1, -1, -1, 1, rep(-1, length(lay$to_j[[j]]))), "=", -Dmat[t, j])
  }

  A <- methods::as(Matrix::sparseMatrix(i = I, j = J, x = X, dims = c(r, nVars)), "CsparseMatrix")
  solve_lp(list(obj = obj, A = A, rhs = rhs, sense = sns, modelsense = "min", vtype = NULL))$objval
}

# Solve one greedy stage LP: choose moves minimising current transport + spot +
# holding/shortage of the resulting inventory (single-step lookahead). Returns
# the moves and the resulting (outgoing) inventory.
.stage_decision <- function(lay, t, ein, xpin, xmin, Q_t, D_t) {
  nI <- lay$nI; nJ <- lay$nJ; nM <- lay$nMoves
  n <- nM + nI + 2L * nJ
  c_eo <- nM + seq_len(nI); c_xp <- nM + nI + seq_len(nJ); c_xm <- nM + nI + nJ + seq_len(nJ)

  obj <- numeric(n)
  obj[seq_len(lay$nL_)] <- lay$transport
  obj[lay$nL_ + seq_len(nM - lay$nL_)] <- .spot_obj(lay, t)
  obj[c_eo] <- lay$entry_hold; obj[c_xp] <- lay$exit_hold; obj[c_xm] <- lay$short

  I <- integer(0); J <- integer(0); X <- numeric(0); rhs <- numeric(0); sns <- character(0)
  r <- 0L
  row <- function(cols, x, sense, b) {
    r <<- r + 1L; I <<- c(I, rep(r, length(cols))); J <<- c(J, cols)
    X <<- c(X, rep_len(x, length(cols)))
    rhs <<- c(rhs, b); sns <<- c(sns, sense)
  }
  for (k in seq_len(lay$nSC))  row(lay$contract_cols[[k]], 1, "<=", lay$Cb[t, k])
  for (c in seq_len(lay$nSCo)) row(lay$spot_cols[[c]],     1, "<=", lay$Co[t, c])
  for (i in seq_len(nI))       row(lay$from_i[[i]], 1, "<=", ein[i] + Q_t[i])
  for (j in seq_len(nJ))       row(lay$to_j[[j]],   1, "<=", lay$exit_cap[j] - xpin[j])
  for (i in seq_len(nI))       row(c(c_eo[i], lay$from_i[[i]]),
                                   c(1, rep(1, length(lay$from_i[[i]]))), "=", ein[i] + Q_t[i])
  for (j in seq_len(nJ))       row(c(c_xp[j], c_xm[j], lay$to_j[[j]]),
                                   c(1, -1, rep(-1, length(lay$to_j[[j]]))), "=", xpin[j] - xmin[j] - D_t[j])

  A <- methods::as(Matrix::sparseMatrix(i = I, j = J, x = X, dims = c(r, n)), "CsparseMatrix")
  x <- pmax(0, solve_lp(list(obj = obj, A = A, rhs = rhs, sense = sns,
                             modelsense = "min", vtype = NULL))$x)
  list(move = x[seq_len(nM)], eo = x[c_eo], xp = x[c_xp], xm = x[c_xm])
}

# Cumulative cost of the myopic (greedy, no look-ahead) policy on one trajectory,
# accounted with the SDDP cost model.
.myopic_cost <- function(lay, Qmat, Dmat, init) {
  ein <- init$entry; xpin <- init$exitp; xmin <- init$exitm
  total <- 0
  for (t in seq_len(lay$tau)) {
    s <- .stage_decision(lay, t, ein, xpin, xmin, Qmat[t, ], Dmat[t, ])
    total <- total + .stage_cost(lay, t, ein, xpin, xmin, s$move)
    ein <- s$eo; xpin <- s$xp; xmin <- s$xm
  }
  total
}
