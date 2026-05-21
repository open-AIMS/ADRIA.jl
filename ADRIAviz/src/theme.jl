const COLORS::Dict{Symbol,Union{Symbol,String}} = Dict(
    :RCP45 => :darkblue,
    :RCP60 => :seagreen,
    :RCP85 => :orangered,
    :scenarios => :steelblue,
    :interventions => :steelblue,
    :counterfactual => :red,
    :unguided => :lawngreen,
    :guided => :dodgerblue,
    :target => "#1f78b4",
    :non_target => "#ff7f00",
    :order => :dodgerblue,
    :topsis => :deepskyblue4,
    :vikor => :midnightblue
)

function labels(group_names::Vector{Symbol})::Vector{String}
    return [uppercasefirst(replace(string(name), '_' => ' ')) for name in group_names]
end
