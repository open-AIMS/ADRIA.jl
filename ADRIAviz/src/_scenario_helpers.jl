const _SCENARIO_TYPES = [:counterfactual, :unguided, :guided]

function _no_seed(scenarios::DataFrame)::BitVector
    return dropdims(
        sum(Matrix(scenarios[:, contains.(names(scenarios), "N_seed")]); dims=2); dims=2
    ) .== 0
end

function _counterfactual(scenarios::DataFrame)::BitVector
    no_seed = _no_seed(scenarios)
    no_fog = scenarios.fogging .== 0
    no_SRM = scenarios.SRM .== 0
    no_mc = scenarios.N_mc_settlers .== 0
    no_mcb = scenarios.mcb_duration .== 0
    return no_seed .& no_fog .& no_SRM .& no_mc .& no_mcb
end

function _unguided(scenarios::DataFrame)::BitVector
    has_seed = .!_no_seed(scenarios)
    has_shade = (scenarios.fogging .> 0) .| (scenarios.SRM .> 0)
    has_mc_corals = scenarios.N_mc_settlers .> 0
    has_mcb = scenarios.mcb_duration .> 0
    return (scenarios.guided .== 0) .& (has_seed .| has_shade .| has_mc_corals .| has_mcb)
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
