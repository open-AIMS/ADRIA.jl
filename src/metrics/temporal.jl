"""Scenario outcomes over time and space.

Provides summary statistics across selected scenarios.
"""

import Interpolations: GriddedInterpolation


function summarize_data(data::AbstractArray{<:Real}, dims::Tuple{Vararg{Int64}}, timesteps=(:))::Dict{Symbol, AbstractArray{<:Real}}
    summarized::Dict{Symbol, AbstractArray{<:Real}} = Dict(Symbol(f) => dropdims(f(data, dims=dims), dims=dims)
                                                           for f in [mean, median, std, minimum, maximum])

    # Calculate quantiles (doesn't support `dims` so have to loop directly)
    q_series::Array{Float32} = fill(0.0, size(data, 1), 8)
    qs::Array{Float32} = Float32[0.025, 0.125, 0.25, 0.375, 0.625, 0.75, 0.875, 0.975]
    @inbounds Threads.@threads for i in 1:size(data[timesteps=timesteps], 1)
        q_series[i, :] = quantile(vec(collect(selectdim(data, 1, i))), qs)
    end

    target_keys = Symbol[:lower_95, :lower_75, :lower_50, :lower_25,
                         :upper_25, :upper_50, :upper_75, :upper_95]
    @inbounds for (i, k) in enumerate(target_keys)
        summarized[k] = q_series[:, i]
    end

    return summarized
end


function summarize_rci(rs::ResultSet, s_ids, dims::Tuple{Vararg{Int64}}=(4,3,2))::Dict{Symbol, AbstractArray{<:Real}}
    X::AbstractArray{<:Real} = rs.raw[scenarios=s_ids]

    # Divide across sites by the max possible proportional coral cover
    rc::AbstractArray{<:Real} = relative_cover(X)
    rc .= mapslices((s) -> s ./ (rs.site_max_coral_cover / 100.0), rc, dims=2)

    E::AbstractArray{<:Real} = Array(Float32[])
    SV::AbstractArray{<:Real} = shelter_volume(X, rs.inputs[s_ids, :])
    juv::AbstractArray{<:Real} = juveniles(X)

    rci::AbstractArray{<:Real} = reef_condition_index(rc, E, SV, juv)
    return summarize_data(rci, dims)
end


"""
    summarize_total_cover(data::AbstractArray{<:Real}, areas::AbstractArray{<:Real})::Dict{Symbol,AbstractArray{<:Real}}
    summarize_total_cover(rs::ResultSet, dims::Tuple=(4,3,2))::Dict{Symbol,AbstractArray{<:Real}}

Calculate summarized total cover.
"""
function summarize_total_cover(data::AbstractArray{<:Real}, areas::AbstractArray{<:Real})::Dict{Symbol,AbstractArray{<:Real}}
    rc::AbstractArray{<:Real} = total_cover(data, areas)
    return summarize_data(rc, (4,3,2))
end
function summarize_total_cover(rs::ResultSet)::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_total_cover(rs.raw, rs.site_area)
end


"""
    summarize_relative_cover(data::AbstractArray{<:Real}, dims::Tuple=(4,3,2))::Dict{Symbol,AbstractArray{<:Real}}
    summarize_relative_cover(rs::ResultSet, dims::Tuple=(4,3,2))::Dict{Symbol,AbstractArray{<:Real}}

Calculate summarized relative cover.
"""
function summarize_relative_cover(data::AbstractArray{<:Real}, dims::Tuple{Vararg{Int64}}=(4,3,2))::Dict{Symbol,AbstractArray{<:Real}}
    rc::AbstractArray{<:Real} = relative_cover(data)
    return summarize_data(rc, dims)
end
function summarize_relative_cover(rs::ResultSet, dims::Tuple{Vararg{Int64}}=(4,3,2))::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_relative_cover(rs.raw, dims)
end


"""
    summarize_coral_evenness(data::AbstractArray{<:Real}, dims::Tuple=(4,3,2))::Dict{Symbol,AbstractArray{<:Real}}
    summarize_coral_evenness(rs::ResultSet, dims::Tuple=(4,3,2))::Dict{Symbol,AbstractArray{<:Real}}

Calculate summarized coral evenness.
"""
function summarize_coral_evenness(data::AbstractArray{<:Real})::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_data(coral_evenness(data), (4,3,2))
end
function summarize_coral_evenness(rs::ResultSet)::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_coral_evenness(rs.raw)
end


function summarize_shelter_volume(rs::ResultSet, dims=(4,3,2))::Dict{Symbol, AbstractArray{<:Real}}
    return summarize_data(shelter_volume(rs.raw, rs.inputs), dims)
end


function summarize_raw(data::AbstractArray{<:Real}, dims::Tuple{Vararg{Int64}})::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_data(data, dims)
end
function summarize_raw(rs::ResultSet, dims::Tuple{Vararg{Int64}}=(4,3,2))::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_data(rs.raw, dims)
end


"""
    trajectory_heatmap(data::Matrix{Float64})::Tuple{Vector{Float64}, Vector{Float64}, Matrix{Int64}}

Estimate heatmap of trajectories from a 2D dataset.

# Arguments
- data : An N*D matrix where N is time steps and D is the scenario outcome for the given timestep in N

# Returns
OnlineStats.HeatMap
"""
function trajectory_heatmap(data::AbstractArray{<:Real})::HeatMap
    n_ts::Int64, n_scens::Int64 = size(data)
    o = HeatMap(zip(repeat(1:n_ts, n_scens), data), n_ts)

    return o
end


"""
    trajectory_heatmap_data(data::Matrix{Float64})::Tuple{Vector{Float64}, Vector{Float64}, Matrix{Int64}}

Estimate heatmap of trajectories from a 2D dataset.

# Arguments
- data : An N*D matrix where N is time steps and D is the scenario outcome for the given timestep in N

# Returns
Tuple of xedges, yedges, and bi-dimensional histogram matrix
"""
function trajectory_heatmap_data(data::AbstractArray{<:Real})::Tuple{Vector{Float64},Vector{Float64},Matrix{Int64}}
    o::HeatMap = trajectory_heatmap(data)

    return collect(o.xedges), collect(o.yedges), o.counts
end
# function temporal_histogram(rs::ResultSet)::Tuple{Vector{Float64}, Vector{Float64}, Matrix{Int64}}
