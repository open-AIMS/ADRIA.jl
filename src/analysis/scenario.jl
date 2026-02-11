const SCENARIO_TYPES = [:counterfactual, :unguided, :guided]

function scenario_clusters(clusters::BitVector)::Dict{Symbol,BitVector}
    return Dict(:target => clusters, :non_target => .!clusters)
end
function scenario_clusters(clusters::Vector{Int64})::Dict{Symbol,BitVector}
    return Dict(Symbol("Cluster_$(c)") => clusters .== c for c in unique(clusters))
end

function scenario_rcps(scenarios::DataFrame)::Dict{Symbol,BitVector}
    rcps::Vector{Symbol} = Symbol.(:RCP, Int64.(scenarios[:, :RCP]))
    return Dict(rcp => rcps .== rcp for rcp in unique(rcps))
end

function scenario_types(scenarios::DataFrame)::Dict{Symbol,BitVector}
    return Dict(
        type => eval(type)(scenarios) for
        type in SCENARIO_TYPES if count(eval(type)(scenarios)) != 0
    )
end

function counterfactual(scenarios::DataFrame)::BitVector
    no_seed = _no_seed(scenarios)
    no_fog = scenarios.fogging .== 0
    no_SRM = scenarios.SRM .== 0
    return no_seed .& no_fog .& no_SRM
end

function unguided(scenarios::DataFrame)::BitVector
    has_seed = .!_no_seed(scenarios)
    has_shade = (scenarios.fogging .> 0) .| (scenarios.SRM .> 0)
    return (scenarios.guided .== 0) .& (has_seed .| has_shade)
end

function guided(scenarios::DataFrame)::BitVector
    return .!(counterfactual(scenarios) .| unguided(scenarios))
end

function _no_seed(scenarios::DataFrame)::BitVector
    return dropdims(
        sum(Matrix(scenarios[:, contains.(names(scenarios), "N_seed")]); dims=2); dims=2
    ) .== 0
end
