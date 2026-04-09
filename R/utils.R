# Thin adapter: solve a Gurobi-style LP/MIP model with HiGHS.
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

# Build the model environment expected by TLPR constraint builders.
# Adapted from tsproj/R/sensitivity.R and recourse.R.
build_env <- function(sim) {
  fromx <- outer(1L:sim$nOrigins, (1L:sim$nDestinations - 1L) * sim$nOrigins, "+")
  tox   <- outer((1L:sim$nDestinations - 1L) * sim$nOrigins, 1L:sim$nOrigins, "+")

  env        <- new.env(parent = emptyenv())
  env$tau    <- sim$tau - 1L
  env$nI     <- sim$nOrigins
  env$nJ     <- sim$nDestinations
  env$I_     <- 1L:env$nI
  env$J_     <- 1L:env$nJ
  env$L      <- TLPR::CartesianProductX(env$I_, env$J_)
  env$R      <- max(sim$exit_capacity)
  env$nL     <- env$nI * env$nJ
  nSC        <- if (!is.null(sim$nSpotCarriers)) sim$nSpotCarriers else sim$nCarriers
  env$nCS    <- sim$nCarriers
  env$nCO    <- nSC
  env$Cb     <- matrix(sim$carrier_capacity[1L:(sim$nCarriers * sim$tau)],                    nrow = sim$tau)
  env$Co     <- matrix(sim$carrier_capacity[(sim$nCarriers * sim$tau + 1L):((sim$nCarriers + nSC) * sim$tau)], nrow = sim$tau)
  env$nLc    <- sim$nLc
  env$L_     <- sim$Ldx
  env$nL_    <- length(sim$Ldx)
  env$nvars  <- length(sim$Ldx) + env$nCO * env$nL
  env$CS     <- 1L:env$nCS
  env$from_i <- apply(fromx, 1L, function(ldx) which(c(sim$Ldx, rep(seq(env$nL), env$nCO)) %in% ldx), simplify = FALSE)
  env$to_j   <- apply(tox,   1L, function(ldx) which(c(sim$Ldx, rep(seq(env$nL), env$nCO)) %in% ldx), simplify = FALSE)
  env$CTb    <- sim$transport_coef
  env$CTo    <- matrix(sim$spot_coef, nrow = sim$tau)
  env$alpha  <- c(sim$entry_store_coef, c(rbind(sim$exit_store_coef, sim$exit_short_coef)))
  env
}

# Build a full multi-period LP model from an instance and realised flows.
build_model <- function(env, init_state, Q, D) {
  ccx <- TLPR::carrier_capacity_padded(env)
  tlx <- TLPR::transition_logic(env, q = Q[seq(env$nI)], d = D[seq(env$nJ)])
  slx <- TLPR::storage_limits(env, q = Q[seq(env$nI)])

  obj_ <- c(env$alpha, env$CTb, env$CTo[1L, ], env$alpha)
  A    <- rbind(ccx$A, tlx$A, slx$A)
  rhs  <- c(ccx$rhs, tlx$rhs, slx$rhs)
  sns  <- c(ccx$sense, tlx$sense, slx$sense)

  model <- TLPR::multiperiod_expansion(env, Q, D, A, obj_, rhs, sns)

  offset <- env$nI + 2L * env$nJ
  model$A <- rbind(
    Matrix::spMatrix(
      ncol = ncol(model$A), nrow = offset,
      i = 1L:offset, j = 1L:offset,
      x = rep(1L, offset)),
    model$A)
  model$sense      <- c(rep("=", offset), model$sense)
  model$rhs        <- c(init_state, model$rhs)
  model$modelsense <- "min"
  model$vtype      <- rep("I", ncol(model$A))
  model$A          <- methods::as(model$A, "CsparseMatrix")
  model
}
