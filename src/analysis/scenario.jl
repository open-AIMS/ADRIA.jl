const SCENARIO_TYPES = [:counterfactual, :unguided, :guided]

function scenario_clusters(clusters::BitVector)::Dict{Symbol,BitVector}
    return Dict(:target => clusters, :non_target => .!clusters)
end
function scenario_clusters(clusters::Vector{Int64})::Dict{Symbol,BitVector}
    return Dict(
        Symbol("Cluster_$(cluster)") => clusters .== cluster for cluster in unique(clusters)
    )
end

function scenario_rcps(
    rs_inputs::DataFrame; scenarios::Union{UnitRange,Colon,Vector{Int64},BitVector}=(:)
)::Dict{Symbol,BitVector}
    rs_rcps::Vector{Symbol} = Symbol.(:RCP, Int64.(rs_inputs[scenarios, :RCP]))
    unique_rcps = unique(rs_rcps)
    return Dict(rcp => rs_rcps .== rcp for rcp in unique_rcps)
end

function scenario_types(
    rs_inputs::DataFrame; scenarios::Union{UnitRange,Colon,Vector{Int64},BitVector}=(:)
)::Dict{Symbol,BitVector}
    return Dict(
        type => eval(type)(rs_inputs[scenarios, :]) for
        type in SCENARIO_TYPES if count(eval(type)(rs_inputs[scenarios, :])) != 0
    )
end

function counterfactual(rs_inputs::DataFrame)::BitVector
    no_seed = _no_seed(rs_inputs)
    no_fog = rs_inputs.fogging .== 0
    no_SRM = rs_inputs.SRM .== 0
    return no_seed .& no_fog .& no_SRM
end

function unguided(rs_inputs::DataFrame)::BitVector
    has_seed = .!_no_seed(rs_inputs)
    has_shade = (rs_inputs.fogging .> 0) .| (rs_inputs.SRM .> 0)
    return (rs_inputs.guided .== 0) .& (has_seed .| has_shade)
end

function guided(rs_inputs::DataFrame)::BitVector
    return .!(counterfactual(rs_inputs) .| unguided(rs_inputs))
end

function _no_seed(rs_inputs::DataFrame)::BitVector
    return (rs_inputs.N_seed_TA .== 0) .&
           (rs_inputs.N_seed_CA .== 0) .&
           (rs_inputs.N_seed_SM .== 0)
end
