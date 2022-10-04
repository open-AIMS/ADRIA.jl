function scenario_type(rs)
    inputs = rs.inputs

    no_seed = (inputs.seed_TA .== 0) .& (inputs.seed_CA .== 0)
    no_fog = inputs.fogging .== 0
    no_SRM = inputs.SRM .== 0
    counterfactual = no_seed .& no_fog .& no_SRM

    has_seed = (inputs.seed_TA .> 0) .| (inputs.seed_CA .> 0)
    has_shade = (inputs.fogging .> 0) .| (inputs.SRM .> 0)
    unguided = (inputs.guided .== 0) .& (has_seed .| has_shade)

    # Guided scenarios must be the inverse of unguided/counterfactual scenarios
    guided = Bool.(ones(Int64, size(inputs, 1)) .⊻ (counterfactual .| unguided))

    return (counterfactual=counterfactual, unguided=unguided, guided=guided)
end

function scenario_colors(rs, weight::Float64, hide::BitVector)
    inputs = rs.inputs
    color_map = repeat([(:blue, weight)], size(inputs, 1))
    scen_type = scenario_type(rs)
    counterfactual = scen_type.counterfactual
    unguided = scen_type.unguided

    color_map[counterfactual] .= ((:red, weight), )
    color_map[unguided] .= ((:green, weight), )
    color_map[hide] .= ((:white, 0.0), )

    return color_map
end
function scenario_colors(rs, weight::Float64)
    inputs = rs.inputs

    color_map = repeat([(:blue, weight)], size(inputs, 1))

    scen_type = scenario_type(rs)
    counterfactual = scen_type.counterfactual
    unguided = scen_type.unguided

    color_map[counterfactual] .= ((:red, weight), )
    color_map[unguided] .= ((:green, weight), )

    return color_map
end
function scenario_colors(rs)
    return scenario_colors(rs, 0.1)
end