const _SCENARIO_TYPES = [:counterfactual, :unguided, :guided]

function _no_CAq(scenarios::DataFrame)::BitVector
    return dropdims(
        sum(Matrix(scenarios[:, contains.(names(scenarios), "iv_CAq_N")]); dims=2); dims=2
    ) .== 0
end

function _counterfactual(scenarios::DataFrame)::BitVector
    no_CAq = _no_CAq(scenarios)
    no_Fog = scenarios.iv_Fog .== 0
    no_Shd = scenarios.iv_Shd .== 0
    no_LvM = scenarios.iv_LvM_N_settlers .== 0
    no_MCB = scenarios.iv_MCB_duration .== 0
    return no_CAq .& no_Fog .& no_Shd .& no_LvM .& no_MCB
end

function _unguided(scenarios::DataFrame)::BitVector
    has_CAq = .!_no_CAq(scenarios)
    has_Shd = (scenarios.iv_Fog .> 0) .| (scenarios.iv_Shd .> 0)
    has_LvM = scenarios.iv_LvM_N_settlers .> 0
    has_MCB = scenarios.iv_MCB_duration .> 0
    return (scenarios.guided .== 0) .& (has_CAq .| has_Shd .| has_LvM .| has_MCB)
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
