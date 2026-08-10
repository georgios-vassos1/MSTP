using JuMP, HiGHS, SDDP, LinearAlgebra, Statistics

# Capacity optimization (heavy compute; lives in the Julia core, driven from R).
# Two capacity plans are compared against the default x0:
#   * LP proxy  — projected-gradient descent on the duals of a mean-demand
#     (q_i = d_j = lambda) multi-period LP (the best a variance-ignoring planner
#     can do);
#   * SDDP dual — projected-gradient descent on the SDDP stage capacity duals,
#     which average over the demand distribution.
# All randomness (training support/trajectories, evaluation scenarios) is
# supplied by R, so results are reproducible. Ported from rerun/capopt_crn.jl
# and rerun/finalize_correlation.jl (the pre-refactor oracle).

# Copy a config with a new capacity vector (derived index structures recomputed).
function with_capacity(config::HyperParams, x::Vector{Float64})
    HyperParams(
        tau = config.tau, nOrigins = config.nOrigins, nDestinations = config.nDestinations,
        nCarriers = config.nCarriers, nSpotCarriers = config.nSpotCarriers,
        Bids = config.Bids, Winners = config.Winners,
        entry_stock_0 = config.entry_stock_0, exit_stock_0 = config.exit_stock_0,
        exit_short_0 = config.exit_short_0, entry_capacity = config.entry_capacity,
        exit_capacity = config.exit_capacity, entry_store_coef = config.entry_store_coef,
        exit_store_coef = config.exit_store_coef, exit_short_coef = config.exit_short_coef,
        transport_coef = config.transport_coef, spot_coef = config.spot_coef,
        carrier_capacity = Int64.(round.(x)), lambda = config.lambda,
        corrmat = config.corrmat, n_scenarios = config.n_scenarios)
end

# Duals of the capacity constraints of the deterministic mean-demand LP at x.
# Same stage cost/constraints as inst/julia/model.jl but with demand fixed at
# the Poisson mean (lambda). Non-positive duals (<= constraint of a Min).
function det_lp_duals(config::HyperParams, x::Vector{Float64})
    T = config.tau; nO = config.nOrigins; nD = config.nDestinations
    nK = config.nCarriers + config.nSpotCarriers; lam = config.lambda
    m = Model(HiGHS.Optimizer); JuMP.set_silent(m)
    @variable(m, mv[1:(config.nLanes + config.nSpotLanes), 1:T] >= 0)
    @variable(m, ent[1:nO, 1:T+1] >= 0)
    @variable(m, exp_[1:nD, 1:T+1] >= 0)
    @variable(m, exm[1:nD, 1:T+1] >= 0)
    for i in 1:nO; JuMP.fix(ent[i, 1], config.entry_stock_0[i]; force = true); end
    for j in 1:nD
        JuMP.fix(exp_[j, 1], config.exit_stock_0[j]; force = true)
        JuMP.fix(exm[j, 1], config.exit_short_0[j]; force = true)
    end
    @constraint(m, cap[k = 1:nK, t = 1:T], sum(mv[config.CarrierIdx[k], t]) <= x[(k-1)*T + t])
    for t in 1:T
        for i in 1:nO
            @constraint(m, sum(mv[config.from_SG[i], t]) + sum(mv[config.from_SP[i], t]) <= ent[i, t] + lam[i])
            @constraint(m, ent[i, t+1] == ent[i, t] + lam[i] - sum(mv[config.from_SG[i], t]) - sum(mv[config.from_SP[i], t]))
        end
        for j in 1:nD
            @constraint(m, exp_[j, t] + sum(mv[config.to_SG[j], t]) + sum(mv[config.to_SP[j], t]) <= config.exit_capacity[j])
            @constraint(m, exp_[j, t+1] - exm[j, t+1] == exp_[j, t] - exm[j, t] + sum(mv[config.to_SG[j], t]) + sum(mv[config.to_SP[j], t]) - lam[nO + j])
        end
    end
    obj = AffExpr(0.0)
    for t in 1:T
        for i in 1:nO; add_to_expression!(obj, config.entry_store_coef[i], ent[i, t]); end
        for j in 1:nD
            add_to_expression!(obj, config.exit_store_coef[j], exp_[j, t])
            add_to_expression!(obj, config.exit_short_coef[j], exm[j, t])
        end
        add_to_expression!(obj, sum(config.transport_coef[1:config.nLanes] .* mv[1:config.nLanes, t]))
        kdx = vec(((1:config.nSpotCarriers)' .- 1) .* nO .* nD .* T .+ (config.SpotIdx[1] .- 1) .* T .+ t)
        add_to_expression!(obj, sum(config.spot_coef[kdx] .* mv[(config.nLanes+1):(config.nLanes+config.nSpotLanes), t]))
    end
    @objective(m, Min, obj); JuMP.optimize!(m)
    d = zeros(nK * T)
    for k in 1:nK, t in 1:T; d[(k-1)*T + t] = JuMP.dual(cap[k, t]); end
    d
end

# Normalised projected-gradient descent on the total-cost gradient (grad + v),
# clamped to [1, x0]. source = :lp uses det_lp_duals; :sddp cold-trains SDDP on
# the R-supplied support/forward and reads simulate_cap_duals on cap_oob.
function pgd(config::HyperParams, v::Float64,
             support::Vector{Vector{Vector{Float64}}}, forward, cap_oob;
             source::Symbol, outer::Int64, cold::Int64, a0::Float64 = 20.0)
    x = Float64.(config.carrier_capacity); xmax = copy(x)
    for k in 1:outer
        grad = if source == :lp
            det_lp_duals(config, x) .+ v
        else
            cf = with_capacity(config, x)
            simulate_cap_duals(train_model(cf, support, forward, cold), cf, cap_oob)["duals"] .+ v
        end
        x = clamp.(x .- (a0 / sqrt(k)) .* (grad ./ (norm(grad) + 1e-12)), 1.0, xmax)
    end
    x
end

# Mean operational cost of an `iters`-iteration policy at capacity x on the
# shared out-of-sample set `oob`.
function op_on(config::HyperParams, x::Vector{Float64},
               support::Vector{Vector{Vector{Float64}}}, forward, oob, iters::Int64)
    cf = with_capacity(config, x); model = train_model(cf, support, forward, iters)
    sims = SDDP.simulate(model, length(oob), [:move];
                         sampling_scheme = SDDP.Historical(oob), parallel_scheme = SDDP.Serial())
    mean(sum(node[:stage_objective] for node in sim) for sim in sims)
end

# Full regime comparison (x0 vs LP proxy vs SDDP duals) under common random
# numbers. R supplies flat arrays; this reshapes and returns total costs and the
# SDDP advantage in percentage points. `v` is the per-unit reservation price.
function capopt_regime_flat(config::HyperParams,
                            support_flat, fwd_flat, cap_oob_flat, eval_oob_flat,
                            tau, n_scen, dim, n_fwd, n_capoob, n_evaloob;
                            v, outer, cold, eval_iters)
    tau = Int64(tau); n_scen = Int64(n_scen); dim = Int64(dim)
    support  = build_support(Vector{Float64}(support_flat), tau, n_scen, dim)
    forward  = build_historical(Vector{Float64}(fwd_flat),  Int64(n_fwd),     tau, dim)
    cap_oob  = build_historical(Vector{Float64}(cap_oob_flat),  Int64(n_capoob),  tau, dim)
    eval_oob = build_historical(Vector{Float64}(eval_oob_flat), Int64(n_evaloob), tau, dim)
    v = Float64(v); outer = Int64(outer); cold = Int64(cold); eval_iters = Int64(eval_iters)

    vvec = fill(v, length(config.carrier_capacity))
    x0  = Float64.(config.carrier_capacity)
    xlp = pgd(config, v, support, forward, cap_oob; source = :lp,   outer = outer, cold = cold)
    xsd = pgd(config, v, support, forward, cap_oob; source = :sddp, outer = outer, cold = cold)

    t0 = op_on(config, x0,  support, forward, eval_oob, eval_iters) + sum(vvec .* x0)
    tl = op_on(config, xlp, support, forward, eval_oob, eval_iters) + sum(vvec .* xlp)
    ts = op_on(config, xsd, support, forward, eval_oob, eval_iters) + sum(vvec .* xsd)

    Dict(
        "x0_total" => t0, "lp_total" => tl, "sddp_total" => ts,
        "lp_pct"   => 100 * (tl - t0) / t0,
        "sddp_pct" => 100 * (ts - t0) / t0,
        "adv_pp"   => 100 * (tl - ts) / t0,
        "x0_mean"  => mean(x0), "xlp_mean" => mean(xlp), "xsd_mean" => mean(xsd),
    )
end
