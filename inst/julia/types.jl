using DataStructures  # kept only for OrderedDict-free alternative below

# Immutable parameter struct — all fields use plain Julia types that map
# cleanly to/from R via JuliaCall (no OrderedDict, no Symbol-keyed Dict).
struct HyperParams
    tau              ::Int64
    nOrigins         ::Int64
    nDestinations    ::Int64
    nLanes           ::Int64                  # total contracted lanes (derived)
    nSpotLanes       ::Int64                  # total spot lanes (derived)
    nCarriers        ::Int64
    nSpotCarriers    ::Int64
    Bids             ::Vector{Vector{Int64}}  # bids[b] = lane indices covered by bid b
    Winners          ::Vector{Vector{Int64}}  # winners[k] = bid indices won by carrier k
    CarrierIdx       ::Vector{Vector{Int64}}  # move variable indices per carrier (derived)
    SpotIdx          ::Vector{Vector{Int64}}  # spot lane indices per spot carrier (derived)
    from_SG          ::Vector{Vector{Int64}}  # contracted move indices leaving origin i (derived)
    from_SP          ::Vector{Vector{Int64}}  # spot move indices leaving origin i (derived)
    to_SG            ::Vector{Vector{Int64}}  # contracted move indices arriving at dest j (derived)
    to_SP            ::Vector{Vector{Int64}}  # spot move indices arriving at dest j (derived)
    entry_stock_0    ::Vector{Int64}
    exit_stock_0     ::Vector{Int64}
    exit_short_0     ::Vector{Int64}
    entry_capacity   ::Vector{Int64}
    exit_capacity    ::Vector{Int64}
    entry_store_coef ::Vector{Float64}
    exit_store_coef  ::Vector{Float64}
    exit_short_coef  ::Vector{Float64}
    transport_coef   ::Vector{Float64}
    spot_coef        ::Vector{Float64}
    carrier_capacity ::Vector{Int64}
    lambda           ::Vector{Float64}        # Poisson intensity (broadcast to all dims)
    corrmat          ::Matrix{Float64}        # (nOrigins+nDestinations) × (nOrigins+nDestinations)
    n_scenarios      ::Int64                  # training scenarios sampled per iteration
end

# Outer constructor — user provides simple typed inputs; derived index structures are
# computed internally so the caller (R or Julia) never has to build them.
function HyperParams(;
    tau              ::Int64,
    nOrigins         ::Int64,
    nDestinations    ::Int64,
    nCarriers        ::Int64,
    nSpotCarriers    ::Int64               = nCarriers,
    Bids             ::Vector{Vector{Int64}},
    Winners          ::Vector{Vector{Int64}},
    entry_stock_0    ::Vector{Int64},
    exit_stock_0     ::Vector{Int64},
    exit_short_0     ::Vector{Int64},
    entry_capacity   ::Vector{Int64},
    exit_capacity    ::Vector{Int64},
    entry_store_coef ::Vector{Float64},
    exit_store_coef  ::Vector{Float64},
    exit_short_coef  ::Vector{Float64},
    transport_coef   ::Vector{Float64},
    spot_coef        ::Vector{Float64},
    carrier_capacity ::Vector{Int64},
    lambda           ::Vector{Float64},
    corrmat          ::Matrix{Float64},
    n_scenarios      ::Int64               = 10)

    nSpotLanes = nOrigins * nDestinations * nSpotCarriers

    # Contracted lane index (flattened from Bids in winner order)
    ordx    = vcat(Winners...)
    Ldx     = vcat(Bids[ordx]...)
    spotLdx = repeat(1:(nOrigins * nDestinations), nSpotCarriers)

    # Cumulative lane counts per carrier (for block indexing)
    nLc = [0; cumsum([length(vcat(Bids[Winners[k]]...)) for k in 1:nCarriers])]

    # Carrier → move variable index blocks
    strategic_idx = [collect(nLc[k]+1 : nLc[k+1])          for k in 1:nCarriers]
    spot_idx      = [collect((i-1)*(nOrigins*nDestinations)+1 : i*(nOrigins*nDestinations))
                     for i in 1:nSpotCarriers]

    carrier_idx = [strategic_idx; [nLc[end] .+ idx for idx in spot_idx]]

    # Origin/destination lane membership helpers
    origin_lanes(i)  = [(j-1)*nOrigins + i for j in 1:nDestinations]
    dest_lanes(j)    = [(j-1)*nOrigins + i for i in 1:nOrigins]

    from_SG = [findall(l -> l in origin_lanes(i), Ldx)      for i in 1:nOrigins]
    from_SP = [findall(l -> l in origin_lanes(i), spotLdx) .+ nLc[end] for i in 1:nOrigins]
    to_SG   = [findall(l -> l in dest_lanes(j),   Ldx)      for j in 1:nDestinations]
    to_SP   = [findall(l -> l in dest_lanes(j),   spotLdx) .+ nLc[end] for j in 1:nDestinations]

    HyperParams(
        tau, nOrigins, nDestinations,
        nLc[end],   # nLanes
        nSpotLanes,
        nCarriers, nSpotCarriers,
        Bids, Winners,
        carrier_idx, spot_idx,
        from_SG, from_SP, to_SG, to_SP,
        entry_stock_0, exit_stock_0, exit_short_0,
        entry_capacity, exit_capacity,
        entry_store_coef, exit_store_coef, exit_short_coef,
        transport_coef, spot_coef, carrier_capacity,
        lambda, corrmat, n_scenarios
    )
end
