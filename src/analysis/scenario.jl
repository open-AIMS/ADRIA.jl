function counterfactual(rs::ResultSet)::BitVector
    inputs = rs.inputs
    no_seed =
        (inputs.N_seed_TA .== 0) .& (inputs.N_seed_CA .== 0) .& (inputs.N_seed_SM .== 0)
    no_fog = inputs.fogging .== 0
    no_SRM = inputs.SRM .== 0
    return no_seed .& no_fog .& no_SRM
end

function unguided(rs::ResultSet)::BitVector
    inputs = rs.inputs
    has_seed = (inputs.N_seed_TA .> 0) .| (inputs.N_seed_CA .> 0) .| (inputs.N_seed_SM .> 0)
    has_shade = (inputs.fogging .> 0) .| (inputs.SRM .> 0)
    return (inputs.guided .== 0) .& (has_seed .| has_shade)
end

function guided(rs::ResultSet)::BitVector
    inputs = rs.inputs
    return Bool.(ones(Int64, size(inputs, 1)) .⊻ (counterfactual(rs) .| unguided(rs)))
end
