#' Generate a random MSTP instance
#'
#' Produces a named list representing one problem instance: network structure
#' (awarded lanes per carrier), initial stocks, capacities, and cost
#' coefficients. The list can be passed directly to `mstp_config()`.
#'
#' Carrier structure (important for numerical stability). Each origin-destination
#' lane is awarded to a SMALL number of carriers (`carriers_per_lane`, default 2)
#' with DISTINCT base rates spaced across `[rate_lo, rate_hi]`. This mirrors real
#' drayage — a lane is served by a few rate-differentiated carriers — and, unlike
#' assigning many near-identical carriers to every lane, it avoids the massive
#' primal degeneracy that produces unstable LP duals, ghost Benders-cut slopes,
#' and an invalid SDDP lower bound at long horizons / large instances. The total
#' number of carriers can still be large (scalability); only the per-lane overlap
#' is kept small.
#'
#' Carrier capacities are calibrated to `lambda` so the average contracted
#' carrier operates at roughly `target_util` utilisation (default 80%).
#' Spot carriers are priced above contracted rates and given finite capacity for
#' overflow relief.
#'
#' @param tau             Number of time periods.
#' @param nOrigins        Number of origin nodes.
#' @param nDestinations   Number of destination nodes.
#' @param nCarriers       Number of contracted carriers.
#' @param nSpotCarriers   Number of spot carriers (defaults to `nCarriers`).
#' @param carriers_per_lane Number of contracted carriers awarded each lane
#'   (default 2). Keep small (1-3) to avoid LP degeneracy; clamped to
#'   `[1, nCarriers]`.
#' @param seed            Optional integer random seed for reproducibility.
#' @param lambda          Expected demand per origin per period, used to
#'   calibrate carrier capacities. Defaults to 700.
#' @param target_util     Target mean utilisation for contracted carriers.
#' @param util_lo,util_hi Bounds of the per-carrier utilisation multiplier.
#' @param rate_lo,rate_hi Range over which distinct contracted base rates are
#'   spaced (spot rates are drawn above `rate_hi`).
#' @param store_periods   Storage capacity expressed as this many periods of
#'   demand (`ceil(store_periods * lambda)`). Default 20 keeps the buffer
#'   non-binding without the absurdly large (e.g. 1e6) values that wreck the LP
#'   right-hand-side conditioning.
#' @param init_stock     Deterministic initial inventory at every node (default
#'   0, i.e. an empty system). Must NOT be randomised: random per-node initial
#'   stock makes the long-horizon SDDP lower bound numerically invalid.
#' @param entry_capacity,exit_capacity Optional explicit storage caps (override
#'   the `store_periods` default).
#' @return A named list suitable for `mstp_config()`.
#' @export
generate_instance <- function(
    tau               = 12L,
    nOrigins          = 6L,
    nDestinations     = 6L,
    nCarriers         = 20L,
    nSpotCarriers     = nCarriers,
    carriers_per_lane = 2L,
    seed              = NULL,
    lambda            = 700.0,
    target_util       = 0.8,
    util_lo           = 0.5,
    util_hi           = 1.5,
    rate_lo           = 6.0,
    rate_hi           = 8.0,
    store_periods     = 20.0,
    init_stock        = 0L,
    entry_capacity    = NULL,
    exit_capacity     = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  nLanes <- nOrigins * nDestinations
  k      <- max(1L, min(as.integer(carriers_per_lane), nCarriers))

  # Award each lane to k distinct carriers (low overlap → well-conditioned LP).
  carrier_lanes <- vector("list", nCarriers)
  for (lane in seq_len(nLanes)) {
    for (c in sample.int(nCarriers, k)) carrier_lanes[[c]] <- c(carrier_lanes[[c]], lane)
  }
  # Guarantee every carrier has at least one lane.
  for (c in seq_len(nCarriers))
    if (length(carrier_lanes[[c]]) == 0L)
      carrier_lanes[[c]] <- sample.int(nLanes, 1L)

  # One bid per carrier = its awarded lane set; carrier c wins bid c.
  Bids    <- lapply(carrier_lanes, sort)
  Winners <- as.list(seq_len(nCarriers))

  ordx <- unlist(Winners, use.names = FALSE)
  Ldx  <- unlist(Bids[ordx], use.names = FALSE)
  nLc  <- c(0L, cumsum(sapply(Winners, function(w) length(unlist(Bids[w])))))

  # Distinct contracted base rates spaced across [rate_lo, rate_hi]; a small
  # jitter keeps lanes within a carrier from being exactly equal. Because each
  # lane sees only k carriers with well-separated base rates, the per-lane
  # routing choice is non-degenerate.
  base <- if (nCarriers == 1L) (rate_lo + rate_hi) / 2 else
            seq(rate_lo, rate_hi, length.out = nCarriers)
  transport_coef <- unlist(lapply(seq_len(nCarriers),
                     function(c) base[c] + runif(length(Bids[[c]]), -0.1, 0.1)),
                     use.names = FALSE)

  # Capacity calibration (unchanged): mean per-carrier load / utilisation draw.
  mean_load  <- lambda * nOrigins / nCarriers
  util_draw  <- runif(nCarriers, util_lo, util_hi)
  util_draw  <- util_draw * (target_util / mean(util_draw))
  cap_per_cb <- pmax(1L, as.integer(round(mean_load / util_draw)))
  cb_caps    <- as.integer(rep(cap_per_cb, each = tau))

  # Spot carriers: priced above contracted, finite overflow capacity.
  spot_cap_per <- pmax(1L, as.integer(round(lambda * nOrigins / nSpotCarriers * 0.2)))
  spot_coef    <- runif(nLanes * nSpotCarriers * tau, rate_hi, rate_hi + 4.0)

  # Storage caps sized to be non-binding but well-conditioned (not 1e6).
  hold_all <- as.integer(ceiling(store_periods * lambda))
  if (is.null(entry_capacity)) entry_capacity <- rep(hold_all, nOrigins)
  if (is.null(exit_capacity))  exit_capacity  <- rep(hold_all, nDestinations)

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
    # Deterministic initial inventory (default empty). Random per-node initial
    # stock makes the long-horizon SDDP lower bound numerically invalid
    # (confirmed by controlled isolation: same instance, random init → 28.8 SE
    # overshoot; zeroed init → valid). Standard SDDP benchmarks initialise the
    # state to a single principled value, never randomly.
    entry_stock_0   = rep(as.integer(init_stock), nOrigins),
    exit_stock_0    = rep(as.integer(init_stock), nDestinations),
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
