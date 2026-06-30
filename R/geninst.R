#' Generate a random MSTP instance
#'
#' Produces a named list representing one problem instance: network structure
#' (bids, winners), initial stocks, capacities, and cost coefficients. The
#' list can be passed directly to `mstp_config()`.
#'
#' Carrier capacities are calibrated to `lambda` so that the average contracted
#' carrier operates at roughly `target_util` utilisation (default 80%).
#' Heterogeneity is introduced by drawing per-carrier utilisation multipliers
#' from Uniform(util_lo, util_hi), producing a realistic mix of tight and slack
#' carriers.  Spot carriers are priced above contracted rates (reflecting the
#' real Chicago market: contracted median ~$680/TEU, spot ~30-50% higher) and
#' given enough capacity to cover moderate overflow.
#'
#' @param tau             Number of time periods.
#' @param nOrigins        Number of origin nodes.
#' @param nDestinations   Number of destination nodes.
#' @param nCarriers       Number of contracted carriers.
#' @param nSpotCarriers   Number of spot carriers (defaults to `nCarriers`).
#' @param nBids           Number of bids in the combinatorial auction.
#' @param seed            Optional integer random seed for reproducibility.
#' @param lambda          Expected demand per origin per period. Used to
#'   calibrate carrier capacities so that constraints are genuinely informative.
#'   Defaults to 700 (consistent with medium Illinois hub cluster).
#' @param target_util     Target mean utilisation for contracted carriers
#'   (capacity / expected per-carrier load). Default 0.8 means carriers are
#'   slightly under-provisioned on average, ensuring some constraints bind.
#' @param util_lo,util_hi Lower and upper bounds of the per-carrier utilisation
#'   multiplier drawn uniformly. Default Uniform(0.5, 1.5) produces a realistic
#'   mix: some carriers are tight (util < 1, capacity < expected load) and some
#'   have slack (util > 1).
#' @param entry_capacity  Integer vector of length `nOrigins` giving the maximum
#'   entry stock at each origin.  Defaults to 10000 at every origin.  Set to a
#'   large value (e.g. `rep(1e6L, nOrigins)`) when long horizons or high demand
#'   rates would otherwise make the constraint binding.
#' @param exit_capacity   Integer vector of length `nDestinations` giving the
#'   maximum exit stock at each destination.  Same default and guidance as
#'   `entry_capacity`.
#' @return A named list suitable for `mstp_config()`.
#' @export
generate_instance <- function(
    tau             = 12L,
    nOrigins        = 6L,
    nDestinations   = 6L,
    nCarriers       = 20L,
    nSpotCarriers   = nCarriers,
    nBids           = 10L,
    seed            = NULL,
    lambda          = 700.0,
    target_util     = 0.8,
    util_lo         = 0.5,
    util_hi         = 1.5,
    entry_capacity  = rep(10000L, nOrigins),
    exit_capacity   = rep(10000L, nDestinations)
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

  # Expected total demand per period across all origins.
  # Each contracted carrier handles (total demand / nCarriers) on average.
  # target_util < 1 means capacity < expected load on average (constraints bind).
  mean_load  <- lambda * nOrigins / nCarriers
  util_draw  <- runif(nCarriers, util_lo, util_hi)
  # Scale so that mean(util_draw) matches target_util (centre the distribution)
  util_draw  <- util_draw * (target_util / mean(util_draw))
  cap_per_cb <- pmax(1L, as.integer(round(mean_load / util_draw)))
  # Replicate across tau periods (capacity can vary by period; use same draw)
  cb_caps    <- as.integer(rep(cap_per_cb, each = tau))

  # Spot carriers: priced 30-50% above contracted rates (real market premium),
  # capacity set to cover ~20% of expected total demand as overflow relief.
  spot_cap_per <- pmax(1L, as.integer(round(lambda * nOrigins / nSpotCarriers * 0.2)))

  # Transport costs calibrated to Chicago drayage: contracted ~$680/TEU median.
  # In model units, Uniform(6, 8) represents contracted rate.
  # Spot is 30-50% higher: Uniform(8, 12).
  transport_coef <- runif(length(Ldx), 6.0, 8.0)
  spot_coef      <- runif(nLanes * nSpotCarriers * tau, 8.0, 12.0)

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
    entry_capacity  = entry_capacity,
    exit_capacity   = exit_capacity,
    entry_store_coef = rep(20.0, nOrigins),
    exit_store_coef  = rep(10.0, nDestinations),
    exit_short_coef  = rep(30.0, nDestinations),
    transport_coef   = transport_coef,
    spot_coef        = spot_coef,
    carrier_capacity = c(cb_caps, rep(spot_cap_per, nSpotCarriers * tau))
  )
}
