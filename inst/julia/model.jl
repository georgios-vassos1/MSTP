using JuMP, SDDP

# One-stage subproblem for the multistage stochastic transportation model.
# Defines state variables (inventory), decision variables (carrier allocations),
# constraints (capacity, flow balance, transitions), and the stage objective.
function transportation_t(sp::Model, stage::Int64; config::HyperParams)

    ## State variables (inventory at entry and exit hubs)
    @variables(sp, begin
        0 <= entry[i = 1:config.nOrigins],      (SDDP.State, initial_value = config.entry_stock_0[i])
        0 <= exitp[j = 1:config.nDestinations], (SDDP.State, initial_value = config.exit_stock_0[j])
        0 <= exitm[j = 1:config.nDestinations], (SDDP.State, initial_value = config.exit_short_0[j])
    end)

    ## Decision variables
    @variables(sp, begin
        0 <= move[k = 1:(config.nLanes + config.nSpotLanes)]
        inflow[i  = 1:config.nOrigins]
        outflow[j = 1:config.nDestinations]
    end)

    ## Constraints

    # Carrier capacity (contracted + spot); named `cap` for dual extraction
    @constraint(sp, cap[k = 1:(config.nCarriers + config.nSpotCarriers)],
        sum(move[config.CarrierIdx[k]]) <= config.carrier_capacity[(k-1)*config.tau + stage])

    # Entry: total outbound moves <= available stock + inflow
    @constraint(sp, [i = 1:config.nOrigins],
        sum(move[config.from_SG[i]]) + sum(move[config.from_SP[i]]) <= entry[i].in + inflow[i])

    # Exit: total inbound moves + existing stock <= capacity
    @constraint(sp, [j = 1:config.nDestinations],
        exitp[j].in + sum(move[config.to_SG[j]]) + sum(move[config.to_SP[j]]) <= config.exit_capacity[j])

    # Entry inventory transition
    @constraint(sp, [i = 1:config.nOrigins],
        entry[i].out == entry[i].in + inflow[i]
                      - sum(move[config.from_SG[i]]) - sum(move[config.from_SP[i]]))

    # Exit inventory transition (net position: stock minus shortage)
    @constraint(sp, [j = 1:config.nDestinations],
        exitp[j].out - exitm[j].out == exitp[j].in - exitm[j].in
                                     + sum(move[config.to_SG[j]]) + sum(move[config.to_SP[j]])
                                     - outflow[j])

    ## Uncertainty: correlated Poisson inflows and outflows
    Ξ = sample_scenarios(config.n_scenarios, config.lambda, config.corrmat)
    SDDP.parameterize(sp, Ξ) do ξ
        JuMP.fix.(inflow,  ξ[1:config.nOrigins])
        JuMP.fix.(outflow, ξ[config.nOrigins .+ (1:config.nDestinations)])
    end

    ## Spot cost index: spot_coef is indexed [lane, carrier, stage] flattened
    kdx = vec(((1:config.nSpotCarriers)' .- 1) .* config.nOrigins .* config.nDestinations .* config.tau
              .+ (config.SpotIdx[1] .- 1) .* config.tau .+ stage)

    ## Stage objective: holding + shortage + transport + spot costs
    @stageobjective(sp,
        sum(config.entry_store_coef[i] * entry[i].in for i in 1:config.nOrigins)      +
        sum(config.exit_store_coef[j]  * exitp[j].in for j in 1:config.nDestinations) +
        sum(config.exit_short_coef[j]  * exitm[j].in for j in 1:config.nDestinations) +
        sum(config.transport_coef[1:config.nLanes] .* move[1:config.nLanes])           +
        sum(config.spot_coef[kdx] .* move[(config.nLanes+1):(config.nLanes+config.nSpotLanes)])
    )
end
