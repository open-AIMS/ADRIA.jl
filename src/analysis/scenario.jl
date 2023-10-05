SCEN_TYPES = [:counterfactual, :unguided, :guided]

function scenario_types(rs::ResultSet; scenarios=(:))::Dict{Symbol,BitVector}
    return scenario_types(rs.inputs; scenarios=scenarios)
end
function scenario_types(rs_input::DataFrame; scenarios=(:))::Dict{Symbol,BitVector}
    return Dict(
        type => eval(type)(rs_input)[scenarios] for
        type in SCEN_TYPES if sum(eval(type)(rs_input)[scenarios]) != 0
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
