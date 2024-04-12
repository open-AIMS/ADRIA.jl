using Printf: @sprintf
import GeoMakie: alpha

const COLORS::Dict{Symbol,Symbol} = Dict(
    :RCP45 => :darkblue,
    :RCP60 => :seagreen,
    :RCP85 => :orangered,
    :counterfactual => :red,
    :unguided => :lawngreen,
    :guided => :dodgerblue,
    :target => :blue,
    :non_target => :red,
    :order => :dodgerblue,
    :topsis => :deepskyblue4,
    :vikor => :midnightblue,
)

function colors(
    scen_groups::Dict{Symbol,BitVector}
)::Dict{Symbol,Union{Symbol,RGBA{Float32}}}
    group_names = keys(scen_groups)
    if count(group_names .∉ [keys(COLORS)]) > 0
        colormap = categorical_colors(:seaborn_bright, length(group_names))
        return Dict(group => colormap[idx] for (idx, group) in enumerate(group_names))
    else
        return Dict(group => COLORS[group] for group in group_names)
    end
end
function colors(
    scen_groups::Dict{Symbol,BitVector},
    weight::Float64,
)::Vector{Tuple{Symbol,Float64}}
    groups = collect(keys(scen_groups))
    n_scens = length(scen_groups[groups[1]])
    _colors = Vector{Symbol}(undef, n_scens)
    scen_colors = colors(scen_groups)

    for (group, scens) in scen_groups
        _colors[scens] .= [scen_colors[group]]
    end

    return [(c, weight) for c in _colors]
end
function colors(
    scen_groups::Dict{Symbol,BitVector}, weights::Dict{Symbol,Float64}
)::Dict{Symbol,RGBA{Float32}}
    groups = collect(keys(scen_groups))
    scen_colors = colors(scen_groups)

    return Dict(
        group => RGBA{Float32}(
            scen_colors[group].r,
            scen_colors[group].g,
            scen_colors[group].b,
            weights[group],
        ) for group in groups
    )
end

function alphas(scen_groups::Dict{Symbol,BitVector})::Dict{Symbol,Float64}
    return Dict(name => alpha(scens) for (name, scens) in scen_groups)
end

function alpha(scens::BitVector)::Float64
    base_alpha::Float64 = 1.0 / (count(scens) * 0.05)
    return max(min(base_alpha, 0.6), 0.1)
end

function labels(group_names::Vector{Symbol})::Vector{String}
    return [uppercasefirst(replace(string(name), '_' => ' ')) for name in group_names]
end

"""
    scenario_colors!(obs_color, scen_types::NamedTuple, weight::Float64, hide::BitVector)

Hide selected scenarios by changing transparency.
"""
function scenario_colors!(
    obs_color::Observable,
    color_map::Vector,
    scen_types::NamedTuple,
    weight::Float64,
    hide::BitVector,
    guide_toggle_map,
)
    color_map .= obs_color[]
    for (t, l, c) in guide_toggle_map
        if !t.active[]
            continue
        end

        display_color = t.active[] ? c : :gray
        display_weight = t.active[] ? weight : 0.05

        scen_t = scen_types[Symbol(lowercase(l))]
        color_map[scen_t .& .!hide] .= ((display_color, display_weight),)
    end

    color_map[hide] .= ((:white, 0.0),)
    return obs_color[] = color_map
end
