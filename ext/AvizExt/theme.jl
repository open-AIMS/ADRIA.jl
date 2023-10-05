using Printf: @sprintf

const COLORS::Dict{Symbol,Symbol} = Dict(
    :RCP45 => :darkblue,
    :RCP60 => :seagreen,
    :RCP85 => :orangered,
    :counterfactual => :red,
    :unguided => :lawngreen,
    :guided => :dodgerblue,
    :order => :dodgerblue,
    :topsis => :deepskyblue4,
    :vikor => :midnightblue,
)

const BINARY_LABELS::Dict{Bool,String} = Dict(0 => "Non-Target", 1 => "Target")
function scenario_type(rs::ResultSet; scenarios=(:))
    counterfactual = ADRIA.analysis.counterfactual(rs)
    unguided = ADRIA.analysis.unguided(rs)
    guided = ADRIA.analysis.guided(rs)

    return (
        counterfactual=ADRIA.analysis.counterfactual(rs_input)[scenarios],
        unguided=ADRIA.analysis.unguided(rs_input)[scenarios],
        guided=ADRIA.analysis.guided(rs_input)[scenarios],
    )
function scenario_colors(rs::ResultSet, weight::Float64)::Vector{Tuple{Symbol,Float64}}
    return [(c, weight) for c in scenario_colors(rs)]
end
function scenario_colors(rs::ResultSet, weight::Float64, hide::BitVector)
    color_map::Vector{Tuple{Symbol,Float64}} = scenario_colors(rs, weight)

    if length(hide) > 0
        color_map[hide] .= ((:white, 0.0),)
    end

    return color_map
end
function scenario_colors(rs::ResultSet)::Vector{Symbol}
    return scenario_colors(rs.inputs)
end
function scenario_colors(rs_inputs::DataFrame)::Vector{Symbol}
    scen_types::Dict{Symbol,BitVector} = ADRIA.analysis.scenario_types(rs_inputs)
    colors = Vector{Symbol}(undef, size(rs_inputs, 1))

    for (key, value) in scen_types
        colors[value] .= COLORS[key]
    end

    return colors
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

        scen_t = getfield(scen_types, Symbol(lowercase(l)))
        color_map[scen_t .& .!hide] .= ((display_color, display_weight),)
    end

    color_map[hide] .= ((:white, 0.0),)
    return obs_color[] = color_map
end
"""
    cluster_colors(clusters::Vector{Int64}, unique_colors::Vector{RGBA{Float32}})::Vector{RGBA{Float32}}
    cluster_colors(clusters::Vector{Int64})::Vector{RGBA{Float32}}
    cluster_colors(clusters::BitVector)::Vector{RGBA{Float32}}

Gets colors for each clustered scenario.

# Arguments
- `clusters` : Vector with scenario cluster numbers

# Returns
Vector with one color for each clustered scenario.
"""
function cluster_colors(
    clusters::Vector{Int64}, unique_colors::Vector{RGBA{Float32}}
)::Vector{RGBA{Float32}}
    unique_clusters::Vector{Int64} = sort(unique(clusters))
    colors_map::Dict{Int64,RGBA{Float32}} = Dict(
        c => unique_colors[i] for (i, c) in enumerate(unique_clusters)
    )

    colors::Vector{RGBA{Float32}} = Vector{RGBA{Float32}}(undef, length(clusters))
    for (idx_c, cluster) in enumerate(clusters)
        colors[idx_c] = colors_map[cluster]
    end

    return colors
end
function cluster_colors(clusters::Vector{Int64})::Vector{RGBA{Float32}}
    unique_colors::Vector{RGBA{Float32}} = categorical_colors(
        :seaborn_bright, length(unique(clusters))
    )
    return cluster_colors(Int64.(clusters), unique_colors)
end
function cluster_colors(clusters::BitVector)::Vector{RGBA{Float32}}
    if unique(clusters) == [0]
        return fill(parse(Colorant, :red), length(clusters))
    elseif unique(clusters) == [1]
        return fill(parse.(Colorant, :blue), length(clusters))
    end
    unique_colors::Vector{RGBA{Float32}} = parse.(Colorant, [:red, :blue])
    return cluster_colors(Int64.(clusters), unique_colors)
end

"""
    cluster_alphas(clusters::Vector{Int64})::Vector{Float64}
    cluster_alphas(clusters::BitVector)::Vector{Float64}

Get color alphas for each cluster weighted by number of scenarios.

# Arguments
- `clusters` : Vector with scenario cluster ids

# Returns
Vector with one color alpha for each cluster.
"""
function cluster_alphas(clusters::Vector{Int64})::Dict{Int64,Float64}
    alphas::Dict{Int64,Float64} = Dict()

    for (i, cluster) in enumerate((unique(clusters)))
        n_scens::Int64 = count(clusters .== cluster)
        base_alpha::Float64 = 1.0 / (n_scens * 0.05)
        alphas[i] = max(min(base_alpha, 0.6), 0.1)
    end

    return alphas
end
function cluster_alphas(clusters::BitVector)::Dict{Int64,Float64}
    return cluster_alphas(Int64.(clusters))
end

"""
    cluster_labels(clusters::Vector{Int64})::Vector{String}
    cluster_labels(clusters::BitVector)::Vector{String}
    cluster_labels(clusters::Vector{Int64}, data::AbstractVector{<:Real})::Vector{String}
    cluster_labels(clusters::BitVector, data::AbstractVector{<:Real})::Vector{String}
    cluster_labels(cluster_names::Vector{String}, data::AbstractVector{<:Real})::Vector{String}

Get labels for each cluster.

# Arguments
- `clusters` : Vector with scenario cluster numbers

# Returns
Vector of labels for each cluster.
"""
function cluster_labels(clusters::Vector{Int64})::Vector{String}
    return "Cluster " .* string.(unique(clusters))
end
function cluster_labels(clusters::BitVector)::Vector{String}
    return [BINARY_LABELS[cluster] for cluster in unique(clusters)]
end
function cluster_labels(
    clusters::Vector{Int64}, data::AbstractVector{<:Real}
)::Vector{String}
    cluster_names::Vector{String} = "Cluster " .* string.(clusters)
    return cluster_labels(cluster_names, data)
end
function cluster_labels(clusters::BitVector, data::AbstractVector{<:Real})::Vector{String}
    cluster_names::Vector{String} = [BINARY_LABELS[cluster] for cluster in clusters]
    return cluster_labels(cluster_names, data)
end
function cluster_labels(
    cluster_names::Vector{String}, data::AbstractVector{<:Real}
)::Vector{String}
    legend_labels = Vector{String}(undef, length(unique(cluster_names)))
    for (idx, name) in enumerate(unique(cluster_names))
        # TODO Extract to external function and accept as argument
        cluster_metric = mean(data[name .== cluster_names])
        metric_formatted = @sprintf "%.1e" cluster_metric
        legend_labels[idx] = "$(name): $metric_formatted"
    end
    return legend_labels
end

function alphas(scenario_types::Dict{Symbol,BitVector})::Dict{Symbol,Float64}
    return Dict(k => alpha(v) for (k, v) in scenario_types)
end

function alpha(scenario_types::BitVector)::Float64
    base_alpha = (1 - count(scenario_types) / length(scenario_types))
    return max(min(base_alpha, 0.3), 0.1)
end
