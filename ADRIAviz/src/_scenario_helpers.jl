const _SCENARIO_TYPES = [:counterfactual, :unguided, :guided]

function _CAq_off(scenarios::DataFrame)::BitVector
    return dropdims(
        sum(Matrix(scenarios[:, contains.(names(scenarios), "iv_CAq_N")]); dims=2); dims=2
    ) .== 0
end

function _counterfactual(scenarios::DataFrame)::BitVector
    iv_CAq_off = _CAq_off(scenarios)
    iv_Fog_off = scenarios.iv_Fog .== 0
    iv_Shd_off = scenarios.iv_Shd .== 0
    iv_LvM_off = scenarios.iv_LvM_N_settlers .== 0
    iv_MCB_off = scenarios.iv_MCB_duration .== 0
    return iv_CAq_off .& iv_Fog_off .& iv_Shd_off .& iv_LvM_off .& iv_MCB_off
end

function _unguided(scenarios::DataFrame)::BitVector
    iv_CAq_on = .!_CAq_off(scenarios)
    iv_Shd_on = (scenarios.iv_Fog .> 0) .| (scenarios.iv_Shd .> 0)
    iv_LvM_on = scenarios.iv_LvM_N_settlers .> 0
    iv_MCB_on = scenarios.iv_MCB_duration .> 0
    return (scenarios.guided .== 0) .& (iv_CAq_on .| iv_Shd_on .| iv_LvM_on .| iv_MCB_on)
end

function _guided(scenarios::DataFrame)::BitVector
    return .!(_counterfactual(scenarios) .| _unguided(scenarios))
end

function _scenario_rcps(scenarios::DataFrame)::Dict{Symbol,BitVector}
    rcps::Vector{Symbol} = Symbol.(:RCP, Int64.(scenarios[:, :RCP]))
    return Dict(rcp => rcps .== rcp for rcp in unique(rcps))
end

function _scenario_types(scenarios::DataFrame)::Dict{Symbol,BitVector}
    type_fns = (
        counterfactual=_counterfactual,
        unguided=_unguided,
        guided=_guided
    )
    return Dict(
        type => type_fns[type](scenarios) for
        type in _SCENARIO_TYPES if count(type_fns[type](scenarios)) != 0
    )
end

function _scenario_clusters(clusters::BitVector)::Dict{Symbol,BitVector}
    return Dict(:target => clusters, :non_target => .!clusters)
end
function _scenario_clusters(clusters::Vector{Int64})::Dict{Symbol,BitVector}
    return Dict(Symbol("Cluster_$(c)") => clusters .== c for c in unique(clusters))
end
