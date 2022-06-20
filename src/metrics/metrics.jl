module metrics

using Interpolations, Statistics, OnlineStats

using DataFrames
import ADRIA: coral_spec, ResultSet


"""
    relative_cover(X::AbstractArray{<:Real})::AbstractArray{<:Real}

# Arguments
- X : Matrix of raw model results
"""
function relative_cover(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    return dropdims(sum(X, dims=2), dims=2)  # sum over all species and size classes
end


"""
    total_cover(X::AbstractArray{<:Real}, site_area::Vector{<:Real})::AbstractArray{<:Real}
    total_cover(rs::ResultSet)::AbstractArray{<:Real}

# Arguments
- X : Matrix of raw model results
- site_area : Vector of site areas, with sites following the same order as given indicated in X.
"""
function total_cover(X::AbstractArray{<:Real}, site_area::Vector{<:Real})::AbstractArray{<:Real}
    return relative_cover(X) .* site_area'
end
function total_cover(rs::ResultSet)::AbstractArray{<:Real}
    return total_cover(rs.raw, rs.site_area)
end


"""
    coral_cover(X::AbstractArray{<:Real})::NamedTuple
    coral_cover(rs::ResultSet)

Converts outputs from scenario runs to relative cover of the four different coral taxa.

# Returns
NamedTuple
    - relative_cover : relative coral cover
    - enhanced_tab_acr : cover of enhanced tabular acropora
    - unenhanced_tab_acr : area covered by unenhanced tabular acropora
    - enhanced_cor_acr : area covered by enhanced corymbose acropora
    - unenhanced_cor_acr : area covered by unenhanced corymbose acropora
    - tab_acr : cover of tabular acropora
    - cor_acr : cover of corymbose acropora
    - small_enc : cover of small encrusting
    - large_mass : cover of large massives
    - juveniles : area covered by juveniles
    - large : area covered by large mature corals
"""
function coral_cover(X::AbstractArray{<:Real})::NamedTuple
    # Relative total coral cover
    rc::AbstractArray{<:Real} = relative_cover(X)  # sum over all species and size classes

    _, _, cs_p::DataFrame = coral_spec()

    screen = (x, idx) -> findall(x .== idx)

    sc1::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 1), :, :, :]
    sc2::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 2), :, :, :]
    sc3::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 3), :, :, :]
    sc4::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 4), :, :, :]

    C1::AbstractArray{<:Real} = sc1 .+ sc2  # enhanced to unenhanced tabular Acropora
    C2::AbstractArray{<:Real} = sc3 .+ sc4  # enhanced to unenhanced corymbose Acropora
    C3::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 5), :, :, :]  # Encrusting and small massives
    C4::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 6), :, :, :]  # Large massives

    # Cover of juvenile corals (< 5cm diameter)
    juv_groups::AbstractArray{<:Real} = X[:, screen(cs_p.class_id, 1), :, :, :] .+ X[:, screen(cs_p.class_id, 2), :, :, :]
    juv_all::AbstractArray{<:Real} = dropdims(sum(juv_groups, dims=2), dims=2)

    large_corals::AbstractArray{<:Real} = X[:, screen(cs_p.class_id, 5), :, :, :] + X[:, screen(cs_p.class_id, 6), :, :, :]
    large_all::AbstractArray{<:Real} = dropdims(sum(large_corals, dims=2), dims=2)

    covers = (relative_cover=rc,
        enhanced_tab_acr=sc1,
        unenhanced_tab_acr=sc2,
        enhanced_cor_acr=sc3,
        unenhanced_cor_acr=sc4,
        tab_acr=C1, cor_acr=C2,
        small_enc=C3, large_mass=C4,
        juveniles=juv_all, large=large_all)

    return covers
end
function coral_cover(rs::ResultSet)::NamedTuple
    return coral_cover(rs.raw)
end


"""
    coral_evenness(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    coral_evenness(rs::ResultSet)::AbstractArray{<:Real}

Calculates evenness across functional coral groups in ADRIA.
Inverse Simpsons diversity indicator.

# Notes
Number of taxa (distinct groups with enhanced lumped with unenhanced) is hardcoded in this function.
"""
function coral_evenness(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    x::Array{<:Real} = min.(max.(X, 0.0), 1.0)
    covers::Array{<:Real} = relative_cover(x)

    # Evenness as a functional diversity metric
    n::Int64 = 4  # number of taxa
    p1::Array{<:Real} = dropdims(sum(covers.tab_acr, dims=2), dims=2) ./ covers
    p2::Array{<:Real} = dropdims(sum(covers.cor_acr, dims=2), dims=2) ./ covers
    p3::Array{<:Real} = dropdims(sum(covers.small_enc, dims=2), dims=2) ./ covers
    p4::Array{<:Real} = dropdims(sum(covers.large_mass, dims=2), dims=2) ./ covers

    sum_psqr::Array{<:Real} = p1 .^ 2 + p2 .^ 2 + p3 .^ 2 + p4 .^ 2  # functional diversity
    simpson_D::Array{<:Real} = 1 ./ sum_psqr  # Hill 1973, Ecology 54:427-432
    return simpson_D ./ n  # Group evenness
end
function coral_evenness(rs::ResultSet)::AbstractArray{<:Real}
    return coral_evenness(rs.raw)
end


"""
    shelter_volume(X::AbstractArray, inputs::DataFrame)
    shelter_volume(rs::ResultSet)

Provide indication of shelter volume.
"""
function shelter_volume(X::AbstractArray{<:Real}, inputs::DataFrame)::AbstractArray{<:Real}
    _, _, cs_p::DataFrame = coral_spec()
    n_corals::Int64 = length(unique(cs_p.taxa_id))

    colony_area_cm2::Array{Float64} = Array{Float64}(inputs[:, contains.(names(inputs), "colony_area_cm2")])'

    sheltervolume_parameters::Array{Float64} = Float64[
        -8.32 1.50;   # tabular from Urbina-Barretto 2021
        -8.32 1.50;   # tabular from Urbina-Barretto 2021
        -7.37 1.34;   # columnar from Urbina-Barretto 2021, assumed similar for corymbose Acropora
        -7.37 1.34;   # columnar from Urbina-Barretto 2021, assumed similar for corymbose Acropora
        -9.69 1.49;   # massives from Urbina-Barretto 2021, assumed similar for encrusting and small massives
        -9.69 1.49]   # massives from Urbina-Barretto 2021,  assumed similar for large massives

    sheltervolume_parameters = repeat(sheltervolume_parameters, n_corals, 1)

    ntsteps::Int64, nspecies::Int64, nsites::Int64, nint::Int64, nreps::Int64 = size(X)

    #  Estimate log colony volume (litres) based on relationship
    #  established by Urbina-Barretto 2021
    logcolony_sheltervolume = sheltervolume_parameters[:, 1] .+ sheltervolume_parameters[:, 2] .* log10.(colony_area_cm2)
    maxlogcolony_sheltervolume = sheltervolume_parameters[:, 1] .+ sheltervolume_parameters[:, 2] .* log10.(maximum(colony_area_cm2, dims=1))

    shelter_volume_colony_litres_per_cm2 = (10.0 .^ logcolony_sheltervolume)
    max_shelter_volume_colony_litres_per_cm2 = (10.0 .^ maxlogcolony_sheltervolume)

    # convert from litres per cm2 to m3 per ha
    cm2_m3::Float64 = (10^-3) * 10^4 * 10^4
    shelter_volume_colony_m3_per_ha::Array{Float64} = shelter_volume_colony_litres_per_cm2 * cm2_m3
    max_shelter_volume_colony_m3_per_ha::Array{Float64} = max_shelter_volume_colony_litres_per_cm2 * cm2_m3

    # calculate shelter volume of groups and size classes and multiply with covers
    sv::Array{Float64} = zeros(ntsteps, nspecies, nsites, nint, nreps)
    for sp::Int64 = 1:nspecies
        sv[:, sp, :, :, :] = (shelter_volume_colony_m3_per_ha[sp] / max_shelter_volume_colony_m3_per_ha[sp]) .* X[:, sp, :, :, :]
    end

    # sum over groups and size classes to estimate total shelter volume per ha
    return dropdims(sum(sv, dims=2), dims=2)
end
function shelter_volume(rs::ResultSet)::AbstractArray{<:Real}
    return shelter_volume(rs.raw, rs.inputs)
end


"""
    reef_condition_index(TC, E, SV, juveniles)
    reef_condition_index(rs)

Translates coral metrics in ADRIA to a Reef Condition Metrics.

# Arguments
- TC        : Total relative coral cover across all groups
- E         : Evenness across four coral groups
- SV        : Shelter volume based coral sizes and abundances
- juveniles : Abundance of coral juveniles < 5 cm diameter

Input dimensions: timesteps, species, sites

# Returns
Dimensions: timesteps, sites, interventions, repeats
"""
function reef_condition_index(TC::T, E::T, SV::T, juveniles::T)::T where {T<:AbstractArray{<:Real}}
    # Compare outputs against reef condition criteria provided by experts

    # These are median values for 7 experts. TODO: draw from distributions
    #  Condition        TC       E       SV      Juv
    # {'VeryGood'}      0.45     0.45    0.45    0.35
    # {'Good'    }      0.35     0.35    0.35    0.25
    # {'Fair'    }      0.25     0.25    0.30    0.25
    # {'Poor'    }      0.15     0.25    0.30    0.25
    # {'VeryPoor'}      0.05     0.15    0.18    0.15

    # Note that the scores for evenness and juveniles are slightly different
    lin_grid::Gridded{Linear{Throw{OnGrid}}} = Gridded(Linear())
    TC_func::GriddedInterpolation{Float64,1,Float64,Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.05, 0.15, 0.25, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0], lin_grid)
    # E_func::GriddedInterpolation{Float64, 1, Float64, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.15, 0.25, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.5, 0.7, 0.9, 1.0], lin_grid)
    SV_func::GriddedInterpolation{Float64,1,Float64,Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.18, 0.30, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.3, 0.5, 0.9, 1.0], lin_grid)
    juv_func::GriddedInterpolation{Float64,1,Float64,Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.15, 0.25, 0.35, 1.0],), Float64[0, 0.1, 0.5, 0.9, 1.0], lin_grid)

    TC_i::T = TC_func.(TC)
    # E_i::T = E_func.(E);
    SV_i::T = SV_func.(SV)
    juv_i::T = juv_func.(juveniles)

    # Original
    # Y = (TC_i + E_i + SV_i + juv_i) ./ 4;

    # Weighted, giving evenness 10#  weight
    # Y = (TC_i*0.3) + (E_i*0.1) + (SV_i*0.3) + (juv_i*0.3);

    # Removing evenness completely
    # Y = mean([TC_i, SV_i, juv_i])
    # Y = (TC_i .+ SV_i .+ juv_i) ./ 3;

    return mean([TC_i, SV_i, juv_i])
end
function reef_condition_index(rs::ResultSet)::Array{<:Real}
    cover::Array{<:Real} = coral_cover(rs)
    TC::Array{<:Real} = total_cover(rs)
    juv::Array{<:Real} = cover.juveniles
    E::Array{<:Real} = coral_evenness(rs)
    SV::Array{<:Real} = shelter_volume(rs)
    return reef_condition_index(TC, E, SV, juv)
end


"""
    summarize_relative_cover(data::AbstractArray{<:Real})::Dict{Symbol,Array{<:Real}}
    summarize_relative_cover(rs::ResultSet)::Dict{Symbol,Array{<:Real}}

Calculate summarized relative cover.
"""
function summarize_relative_cover(data::AbstractArray{<:Real}, dims::Tuple{<:Int}=(4,3,2))::Dict{Symbol,Array{<:Real}}
    cover::Array{<:Real} = relative_cover(data)

    summarized::Dict{Symbol,Array{<:Real}} = Dict(Symbol(f) => dropdims(f(cover, dims=dims), dims=dims)
                                                  for f in [mean, median, std, minimum, maximum])

    # Calculate quantiles (doesn't support `dims` so have to loop directly)
    q_series::Array{Float32} = fill(0.0, size(cover, 1), 8)
    qs::Array{Float32} = Float32[0.025, 0.125, 0.25, 0.375, 0.625, 0.75, 0.875, 0.975]
    @inbounds Threads.@threads for i in 1:size(cover, 1)
        q_series[i, :] = quantile(vec(collect(selectdim(cover, 1, i))), qs)
    end

    target_keys = Symbol[:lower_95, :lower_75, :lower_50, :lower_25,
                         :upper_25, :upper_50, :upper_75, :upper_95]
    @inbounds for (i, k) in enumerate(target_keys)
        summarized[k] = q_series[:, i]
    end

    return summarized
end
function summarize_relative_cover(rs::ResultSet, dims=(4,3,2))::Dict{Symbol,Array{<:Real}}
    return summarize_relative_cover(rs.raw, dims)
end


function summarize_raw(data::AbstractArray{<:Real}, dims_sum::Tuple{Int64})::Dict{Symbol,Array{<:Real}}
    cover::Array{<:Real} = relative_cover(data)

    summarized::Dict{Symbol,Array{<:Real}} = Dict(Symbol(f) => dropdims(f(cover, dims=dims_sum), dims=dims_sum)
                                                  for f in [mean, median, std, minimum, maximum])

    # Calculate quantiles (doesn't support `dims` so have to loop directly)
    q_series::Array{Float32} = fill(0.0, size(cover, 1), 8)
    qs::Array{Float32} = Float32[0.025, 0.125, 0.25, 0.375, 0.625, 0.75, 0.875, 0.975]
    @inbounds Threads.@threads for i in 1:size(cover, 1)
        q_series[i, :] = quantile(vec(collect(selectdim(cover, 1, i))), qs)
    end

    target_keys = Symbol[:lower_95, :lower_75, :lower_50, :lower_25,
        :upper_25, :upper_50, :upper_75, :upper_95]
    @inbounds for (i, k) in enumerate(target_keys)
        summarized[k] = q_series[:, i]
    end

    return summarized
end
function summarize_raw(rs::ResultSet, dims_sum::Tuple{Int64})::Dict{Symbol,Array{<:Real}}
    return summarize_raw(rs.raw, dims_sum)
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

end