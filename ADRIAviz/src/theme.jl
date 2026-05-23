const COLORS::Dict{Symbol,String} = Dict(
    :RCP45 => "#00008b",
    :RCP60 => "#2e8b57",
    :RCP85 => "#ff4500",
    :scenarios => "#4682b4",
    :interventions => "#4682b4",
    :counterfactual => "#ff0000",
    :unguided => "#7cfc00",
    :guided => "#1e90ff",
    :target => "#1f78b4",
    :non_target => "#ff7f00",
    :order => "#1e90ff",
    :topsis => "#00688b",
    :vikor => "#191970"
)

function labels(group_names::Vector{Symbol})::Vector{String}
    return [uppercasefirst(replace(string(name), '_' => ' ')) for name in group_names]
end
