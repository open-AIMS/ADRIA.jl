"""Scenario outcomes over time and space.

Provides summary statistics across selected scenarios.
"""

import Interpolations: GriddedInterpolation


function summarize_trajectory(data::NamedDimsArray)::Dict{Symbol,AbstractArray{<:Real}}
    squash = nothing
    if :locations in NamedDims.dimnames(data)
        squash = (:scenarios, :locations)
    elseif :scenarios in NamedDims.dimnames(data) && (size(data, :scenarios) > 1)
        squash = (:scenarios,)
    end

    if !isnothing(squash)
        summarized::Dict{Symbol,AbstractArray{<:Real}} = Dict(Symbol(f) => collect(dropdims(f(data, dims=squash), dims=squash))
                                                              for f in [mean, median, std, minimum, maximum])
    else
        # Only a single scenario so don't bother doing anything
        summarized = Dict(Symbol(f) => data for f in [mean, median, minimum, maximum])
        summarized[:std] = repeat([0.0], length(data))
    end

    # Calculate quantiles (doesn't support `dims` so have to loop directly)
    q_series::Array{Float32} = fill(0.0, size(data, 1), 8)
    qs::Array{Float32} = Float32[0.025, 0.125, 0.25, 0.375, 0.625, 0.75, 0.875, 0.975]
    @inbounds for i in 1:size(data, 1)
        q_series[i, :] = quantile(vec(collect(selectdim(data, 1, i))), qs)
    end

    target_keys = Symbol[:lower_95, :lower_75, :lower_50, :lower_25,
        :upper_25, :upper_50, :upper_75, :upper_95]
    @inbounds for (i, k) in enumerate(target_keys)
        summarized[k] = q_series[:, i]
    end

    return summarized
end


"""
    summarize_raw(data::NamedDimsArray; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}

Summarize raw data, aggregating the specified dimensions (e.g., `timesteps`, `scenarios`, etc.)
and collapsing given `dims`.
"""
function summarize_raw(data::NamedDimsArray; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_trajectory(slice_results(data; kwargs...))
end
function summarize_raw(rs::ResultSet; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_raw(rs.raw; kwargs...)
end


function summarize_rci(rs::ResultSet; kwargs...)
    rc::AbstractArray{<:Real} = call_metric(relative_cover, X; kwargs...)

    # Divide across locations by the max possible proportional coral cover
    rc .= mapslices((s) -> s ./ (rs.location_max_coral_cover / 100.0), rc, dims=:locations)

    # Empty representation of evenness - currently ignored
    E::AbstractArray{<:Real} = Array(Float32[])

    s_ids = kwargs[:scenarios]
    SV::AbstractArray{<:Real} = call_metric(shelter_volume, X, rs.inputs[s_ids, :]; kwargs...)
    juv::AbstractArray{<:Real} = call_metric(juveniles, X; kwargs...)

    rci::AbstractArray{<:Real} = reef_condition_index(rc, E, SV, juv)

    return summarize_trajectory(rci)
end


"""
    summarize_total_cover(raw::AbstractArray{<:Real}, areas::AbstractArray{<:Real}; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    summarize_total_cover(rs::ResultSet; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}

Calculate summarized total absolute cover.
"""
function summarize_total_cover(raw::NamedDimsArray, areas::AbstractArray{<:Real}; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    locations = haskey(kwargs, :locations) ? kwargs[:locations] : (:)
    tac = call_metric(total_absolute_cover, raw, areas[locations]; kwargs...)
    tac = dropdims(sum(tac, dims=:locations), dims=:locations)
    return summarize_trajectory(tac)
end
function summarize_total_cover(rs::ResultSet; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    tac = dropdims(sum(slice_results(_total_absolute_cover(rs); kwargs...), dims=:locations), dims=:locations)
    return summarize_trajectory(tac)
end


"""
    summarize_relative_cover(rc::AbstractArray{<:Real}, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    summarize_relative_cover(rs::ResultSet, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}

Calculate summarized relative cover.
"""
function summarize_relative_cover(rc::NamedDimsArray; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    rc_sliced = slice_results(rc; kwargs...)
    return summarize_trajectory(rc_sliced)
end
function summarize_relative_cover(rs::ResultSet; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_relative_cover(rs.outcomes[:relative_cover]; kwargs...)
end


"""
    summarize_coral_evenness(raw::AbstractArray{<:Real}, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    summarize_coral_evenness(rs::ResultSet, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}

Calculate summarized coral evenness.
"""
function summarize_coral_evenness(raw::NamedDimsArray; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    ce::AbstractArray{<:Real} = call_metric(coral_evenness, raw; kwargs...)
    return summarize_trajectory(ce)
end
function summarize_coral_evenness(rs::ResultSet; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_coral_evenness(rs.raw; kwargs...)
end


"""
    summarize_absolute_shelter_volume(sv::AbstractArray{<:Real}, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    summarize_absolute_shelter_volume(rs::ResultSet, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}

Calculate summarized coral evenness.
"""
function summarize_absolute_shelter_volume(sv::NamedDimsArray; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_trajectory(slice_results(sv; kwargs...))
end
function summarize_absolute_shelter_volume(rs::ResultSet; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    sv_sliced = slice_results(rs.outcomes[:absolute_shelter_volume]; kwargs...)
    return summarize_trajectory(sv_sliced)
end


"""
    summarize_relative_shelter_volume(sv::AbstractArray{<:Real}, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    summarize_relative_shelter_volume(rs::ResultSet, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}

Calculate summarized coral evenness.
"""
function summarize_relative_shelter_volume(sv::NamedDimsArray; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    return summarize_trajectory(slice_results(sv; kwargs...))
end
function summarize_relative_shelter_volume(rs::ResultSet; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
    sv_sliced = slice_results(rs.outcomes[:relative_shelter_volume]; kwargs...)
    return summarize_trajectory(sv_sliced)
end



function summarize_trajectories(rs, ts_name; kwargs...)
    sliced = slice_results(rs.outcomes[ts_name]; kwargs...)
    return summarize_trajectory(sliced)
end


"""
    trajectory_heatmap(data::Matrix{Float64})::Tuple{Vector{Float64}, Vector{Float64}, Matrix{Int64}}

Estimate heatmap of trajectories from a 2D dataset.

# Arguments
- data : An N*D matrix where N is time steps and D is the scenario outcome for the given timestep in N

# Returns
OnlineStats.HeatMap
"""
function trajectory_heatmap(data::NamedDimsArray)::HeatMap
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
function trajectory_heatmap_data(data::NamedDimsArray)::Tuple{Vector{Float64},Vector{Float64},Matrix{Int64}}
    o::HeatMap = trajectory_heatmap(data)

    return collect(o.xedges), collect(o.yedges), o.counts
end
# function temporal_histogram(rs::ResultSet)::Tuple{Vector{Float64}, Vector{Float64}, Matrix{Int64}}
