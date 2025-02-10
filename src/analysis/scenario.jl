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
    has_fog_col = "fogging" in names(scenarios)
    no_fog = has_fog_col ? scenarios.fogging .== 0 : fill(true, size(scenarios, 1))
    has_SRM_col = "SRM" in names(scenarios)
    no_SRM = has_SRM_col ? scenarios.SRM .== 0 : fill(true, size(scenarios, 1))

    return no_seed .& no_fog .& no_SRM
end

function unguided(scenarios::DataFrame)::BitVector
    has_guided_col = "guided" in names(scenarios)
    # If the results set does not have a guided type column default to unguided
    is_unguided = has_guided_col ? scenarios.guided .== 0 : fill(true, size(scenarios, 1))

    has_fog_col = "fogging" in names(scenarios)
    has_fog = has_fog_col ? scenarios.fogging .> 0 : fill(true, size(scenarios, 1))

    has_SRM_col = "SRM" in names(scenarios)
    has_SRM = has_SRM_col ? scenarios.SRM .> 0 : fill(true, size(scenarios, 1))

    has_seed = .!_no_seed(scenarios)
    has_shade = (has_fog) .| (has_SRM)

    return is_unguided .& (has_seed .| has_shade)
end

function guided(scenarios::DataFrame)::BitVector
    return .!(counterfactual(scenarios) .| unguided(scenarios))
end

function _no_seed(scenarios::DataFrame)::BitVector
    seed_cols::Vector{String} = filter(
        x -> contains(x, "N_seed_"), names(scenarios)
    )

    return BitVector(
        all(collect(scen_row) .== 0) for scen_row in eachrow(scenarios[:, seed_cols])
    )
end
