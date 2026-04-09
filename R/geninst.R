#' Generate a random MSTP instance
#'
#' Produces a named list representing one problem instance: network structure
#' (bids, winners), initial stocks, capacities, and cost coefficients. The
#' list can be passed directly to `mstp_config()`.
#'
#' @param tau           Number of time periods.
#' @param nOrigins      Number of origin nodes.
#' @param nDestinations Number of destination nodes.
#' @param nCarriers     Number of contracted carriers.
#' @param nSpotCarriers Number of spot carriers (defaults to `nCarriers`).
#' @param nBids         Number of bids in the combinatorial auction.
#' @param seed          Optional integer random seed for reproducibility.
#' @return A named list suitable for `mstp_config()`.
#' @export
generate_instance <- function(
    tau           = 12L,
    nOrigins      = 6L,
    nDestinations = 6L,
    nCarriers     = 20L,
    nSpotCarriers = nCarriers,
    nBids         = 10L,
    seed          = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  nLanes <- nOrigins * nDestinations

  Bids    <- replicate(nBids,
               sample(nLanes, size = sample(6L:18L, 1L), replace = FALSE),
               simplify = FALSE)
  Winners <- replicate(nCarriers,
               sample(nBids, size = sample(1L:2L, 1L), replace = FALSE),
               simplify = FALSE)

  ordx <- unlist(Winners, use.names = FALSE)
  Ldx  <- unlist(Bids[ordx], use.names = FALSE)
  nLc  <- c(0L, cumsum(sapply(Winners, function(w) length(unlist(Bids[w])))))

  list(
    tau             = tau,
    nOrigins        = nOrigins,
    nDestinations   = nDestinations,
    nCarriers       = nCarriers,
    nSpotCarriers   = nSpotCarriers,
    Bids            = Bids,
    Winners         = Winners,
    Ldx             = Ldx,
    nLc             = nLc,
    entry_stock_0   = sample(seq(0L,  500L, by = 50L), nOrigins,      replace = TRUE),
    exit_stock_0    = sample(seq(0L, 1000L, by = 50L), nDestinations, replace = TRUE),
    exit_short_0    = rep(0L, nDestinations),
    entry_capacity  = rep(10000L, nOrigins),
    exit_capacity   = rep(10000L, nDestinations),
    entry_store_coef = rep(20.0, nOrigins),
    exit_store_coef  = rep(10.0, nDestinations),
    exit_short_coef  = rep(30.0, nDestinations),
    transport_coef   = runif(length(Ldx), 6.0, 8.0),
    spot_coef        = runif(nLanes * nSpotCarriers * tau, 3.0, 9.0),
    carrier_capacity = c(
      sample(seq(400L, 800L, by = 50L), size = nCarriers     * tau, replace = TRUE),
      rep(40L,                                  nSpotCarriers * tau)
    )
  )
}
