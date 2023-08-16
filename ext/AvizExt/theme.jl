const COLORS::Dict{Symbol,Symbol} = Dict(
    :RCP45 => :darkblue,
    :RCP60 => :seagreen,
    :RCP85 => :orangered,
    :counterfactual => :red,
    :unguided => :lawngreen,
    :guided => :dodgerblue,
    :order => :dodgerblue,
    :topsis => :deepskyblue4,
    :vikor => :midnightblue
)


function scenario_type(rs; scenarios=(:))
    inputs = rs.inputs

    no_seed = (inputs.N_seed_TA .== 0) .& (inputs.N_seed_CA .== 0) .& (inputs.N_seed_SM .== 0)
    no_fog = inputs.fogging .== 0
    no_SRM = inputs.SRM .== 0
    counterfactual = no_seed .& no_fog .& no_SRM

    has_seed = (inputs.N_seed_TA .> 0) .| (inputs.N_seed_CA .> 0) .| (inputs.N_seed_SM .> 0)
    has_shade = (inputs.fogging .> 0) .| (inputs.SRM .> 0)
    unguided = (inputs.guided .== 0) .& (has_seed .| has_shade)

    # Guided scenarios must be the inverse of unguided/counterfactual scenarios
    guided = Bool.(ones(Int64, size(inputs, 1)) .⊻ (counterfactual .| unguided))

    return (counterfactual=counterfactual[scenarios], unguided=unguided[scenarios], guided=guided[scenarios])
end

function scenario_colors(rs, weight::Float64, hide::BitVector)
    color_map = fill((COLORS[:guided], weight), size(rs.inputs, 1))
    scen_type = scenario_type(rs)
    counterfactual = scen_type.counterfactual
    unguided = scen_type.unguided

    color_map[counterfactual] .= ((COLORS[:counterfactual], weight),)
    color_map[unguided] .= ((COLORS[:unguided], weight),)

    if length(hide) > 0
        color_map[hide] .= ((:white, 0.0),)
    end

    return color_map
end
function scenario_colors(rs, weight::Float64)
    color_map = fill((COLORS[:guided], weight), size(rs.inputs, 1))

    scen_type = scenario_type(rs)
    counterfactual = scen_type.counterfactual
    unguided = scen_type.unguided

    color_map[counterfactual] .= ((COLORS[:counterfactual], weight),)
    color_map[unguided] .= ((COLORS[:unguided], weight),)

    return color_map
end
function scenario_colors(rs)
    return scenario_colors(rs, 0.1)
end


"""
    scenario_colors!(obs_color, scen_types::NamedTuple, weight::Float64, hide::BitVector)

Hide selected scenarios by changing transparency.
"""
function scenario_colors!(obs_color::Observable, color_map::Vector, scen_types::NamedTuple, weight::Float64, hide::BitVector, guide_toggle_map)
    color_map .= obs_color[]
    for (t, l, c) in guide_toggle_map
        if !t.active[]
            continue
        end

        display_color = t.active[] ? c : :gray
        display_weight = t.active[] ? weight : 0.05

        scen_t = getfield(scen_types, Symbol(lowercase(l)))
        color_map[scen_t.&.!hide] .= ((display_color, display_weight),)
    end

    color_map[hide] .= ((:white, 0.0),)
    obs_color[] = color_map
end
