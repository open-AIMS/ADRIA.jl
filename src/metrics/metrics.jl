module metrics

using Interpolations, Statistics, OnlineStats, NamedDims

using DataFrames
import ADRIA: coral_spec, ResultSet


abstract type Outcome end


struct Metric{F<:Function,T<:Tuple,S<:String} <: Outcome
    func::F
    dims::T
    unit::S
end
Metric(f, d) = Metric(f, d, "")


"""
Make Metric callable with arbitary arguments that are passed to associated function.
"""
function (f::Metric)(raw, args...; kwargs...)
    try
        return f.func(NamedDimsArray{(:timesteps, :species, :sites, :scenarios)[1:Base.ndims(raw)]}(raw), args...; kwargs...)
    catch
        raw = NamedDimsArray{(f.dims...)}(raw)

        return f.func(raw, args...; kwargs...)
    end
end
function (f::Metric)(rs::ResultSet, args...; kwargs...)
    return f.func(rs, args...; kwargs...)
end


"""
    metric_name(m::Metric)::String

Get name of metric as a string.
"""
function metric_name(m::Metric)::String
    return metric_name(m.func)
end
function metric_name(f::Function)::String
    return join(split(String(Symbol(f))[2:end], "_"), " ")
end

"""
    metric_label(m::Metric)::String
    metric_label(f::Function, unit::String)

Return name of metric in the format: "Title Case [Unit]", suitable for use as a label.

# Example
```julia
m_label = metric_label(scenario_total_cover)
# "Scenario Total Cover [m²]"
```
"""
function metric_label(m::Metric)::String
    return metric_label(m.func, m.unit)
end
function metric_label(f::Function, unit::String)::String
    n = titlecase(metric_name(f))
    if length(unit) > 0
        n *= " [$unit]"
    end

    return n
end


"""
    dims(m::Metric)::Tuple

Get dimension names for a given outcome/metric.
"""
function dims(m::Metric)::Tuple
    return m.dims
end


"""
    ndims(m::Metric)::Int64

Infer the number of dimensions for a given outcome/metric.
"""
function Base.ndims(m::Metric)::Int64
    return length(dims(m))
end


"""
    call_metric(metric, data, args...; timesteps=(:), species=(:), sites=(:), scens=(:))

Convenience method that slices the data in the specified manner.

# Arguments
- `metric` : Function, the metric function to apply to "raw" data.
- `data` : NamedDimsArray, data to pass into `metric`
- `args` : Additional positional arguments to pass into `metric`
- `dims` : dummy keyword argument, not used but defined to allow use with other methods
"""
function call_metric(metric::Function, data::NamedDimsArray{T}, args...; kwargs...) where {T<:Real}
    dims = haskey(kwargs, :dims) ? kwargs[:dims] : nothing
    if isnothing(dims)
        return metric(slice_results(data; kwargs...), args...)
    else
        return metric(slice_results(data; kwargs...), args...; dims=dims)
    end
end


"""
    slice_results(data::NamedDimsArray{T}; timesteps=(:), species=(:), sites=(:), scenarios=(:)) where {T<:Real}

Slice data as indicated.
Dimensions not found in target data are ignored.
"""
function slice_results(data::NamedDimsArray{T}; timesteps=(:), species=(:), sites=(:), scenarios=(:)) where {T<:Real}
    f_dims = (timesteps=timesteps, species=species, sites=sites, scenarios=scenarios)

    s_names = keys(f_dims)
    d_names = NamedDims.dimnames(data)
    common_dims = intersect(s_names, d_names)

    selected_slice = (; zip(common_dims, [getfield(f_dims, k) for k in common_dims])...)
    return data[selected_slice...]
end


""" 
    relative_cover(X::AbstractArray{T}, k_area::Vector{T})::AbstractArray{T} where {T}
    relative_cover(rs::ResultSet)::AbstractArray

# Arguments
- `X` : Matrix of raw model results
"""
function _relative_cover(X::AbstractArray{T}, k_area::Vector{T})::AbstractArray{T} where {T<:Real}
    # sum over all species and size classes
    return (dropdims(sum(X, dims=:species), dims=:species) ./ replace(k_area, 0.0 => 1.0)')
end
function _relative_cover(rs::ResultSet)::AbstractArray
    denom = replace(((rs.site_max_coral_cover ./ 100.0) .* rs.site_area), 0.0 => 1.0)'
    return (rs.outcomes[:total_absolute_cover] ./ denom)
end
relative_cover = Metric(_relative_cover, (:timesteps, :sites, :scenarios))


"""
    total_absolute_cover(X::AbstractArray{T}, site_area::Vector{T})::AbstractArray{T} where {T}
    total_absolute_cover(rs::ResultSet)::AbstractArray{T} where {T}

The Total Absolute Coral Cover.
Sum of proportional area taken up by all corals, multiplied by total site area.

# Arguments
- `X` : Matrix of raw model results
- `site_area` : Vector of site areas, with sites following the same order as given indicated in X.
"""
function _total_absolute_cover(X::AbstractArray{T}, site_area::Vector{T})::AbstractArray{T} where {T<:Real}
    return dropdims(sum(X, dims=:species), dims=:species) .* site_area'
end
function _total_absolute_cover(rs::ResultSet)::AbstractArray
    return rs.outcomes[:total_absolute_cover]
end
total_absolute_cover = Metric(_total_absolute_cover, (:timesteps, :sites, :scenarios), "m²")


"""
    relative_taxa_cover(X::AbstractArray{T,3}) where {T<:Real}

Results grouped by taxa/species.

TODO: Uses hardcoded index values, to be replaced by something more generic.

# Arguments
- `X` : Raw model results for a single scenario

# Returns
Coral cover, grouped by taxa for the given scenario.
"""
function _relative_taxa_cover(X::AbstractArray{T,3}) where {T<:Real}
    nsteps, nspecies, _ = size(X)

    taxa_cover = zeros(nsteps, 6)
    for (taxa_id, grp) in enumerate([i:i+5 for i in 1:6:nspecies])
        # Sum over groups
        taxa_cover[:, taxa_id] = dropdims(mean(sum(X[:, grp, :], dims=2), dims=3), dims=3)
    end

    return taxa_cover
end
function _relative_taxa_cover(rs)::AbstractArray
    return rs.outcomes[:relative_taxa_cover]
end
relative_taxa_cover = Metric(_relative_taxa_cover, (:timesteps, :taxa, :scenarios))


# """
#     coral_cover(X::AbstractArray{T})::NamedTuple where {T<:Real}
#     coral_cover(rs::ResultSet)

# Converts outputs from scenario runs to relative cover of the four different coral taxa.

# # Returns
# NamedTuple
#     - relative_cover : relative coral cover
#     - enhanced_tab_acr : cover of enhanced tabular acropora
#     - unenhanced_tab_acr : area covered by unenhanced tabular acropora
#     - enhanced_cor_acr : area covered by enhanced corymbose acropora
#     - unenhanced_cor_acr : area covered by unenhanced corymbose acropora
#     - tab_acr : cover of tabular acropora
#     - cor_acr : cover of corymbose acropora
#     - small_enc : cover of small encrusting
#     - large_mass : cover of large massives
#     - juveniles : area covered by juveniles
#     - large : area covered by large mature corals
# """
# function coral_cover(X::AbstractArray{<:Real})::NamedTuple
#     rc::AbstractArray{<:Real} = _relative_cover(X)  # sum over all species and size classes

#     _, _, cs_p::DataFrame = coral_spec()

#     screen = (x, idx) -> findall(x .== idx)

#     sc1::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 1), :, :, :]
#     sc2::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 2), :, :, :]
#     sc3::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 3), :, :, :]
#     sc4::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 4), :, :, :]

#     C1::AbstractArray{<:Real} = sc1 .+ sc2  # enhanced to unenhanced tabular Acropora
#     C2::AbstractArray{<:Real} = sc3 .+ sc4  # enhanced to unenhanced corymbose Acropora
#     C3::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 5), :, :, :]  # Encrusting and small massives
#     C4::AbstractArray{<:Real} = X[:, screen(cs_p.taxa_id, 6), :, :, :]  # Large massives

#     # Cover of juvenile corals (< 5cm diameter)
#     juv_all = _absolute_juveniles(X)

#     large_corals::AbstractArray{<:Real} = X[:, screen(cs_p.class_id, 5), :, :, :] + X[:, screen(cs_p.class_id, 6), :, :, :]
#     large_all::AbstractArray{<:Real} = dropdims(sum(large_corals, dims=2), dims=2)

#     covers = (relative_cover=rc,
#         enhanced_tab_acr=sc1,
#         unenhanced_tab_acr=sc2,
#         enhanced_cor_acr=sc3,
#         unenhanced_cor_acr=sc4,
#         tab_acr=C1, cor_acr=C2,
#         small_enc=C3, large_mass=C4,
#         juveniles=juv_all, large=large_all)

#     return covers
# end
# function coral_cover(rs::ResultSet)::NamedTuple
#     return coral_cover(rs.raw)
# end


"""
    relative_juveniles(X::AbstractArray{T})::AbstractArray{T} where {T}
    relative_juveniles(rs::ResultSet)::AbstractArray

Juvenile coral cover relative to total site area.
"""
function _relative_juveniles(X::AbstractArray{T})::AbstractArray{T} where {T<:Real}
    _, _, cs_p::DataFrame = coral_spec()

    # Cover of juvenile corals (< 5cm diameter)
    juv_groups::AbstractArray{<:Real} = X[species=cs_p.class_id .== 1] .+ X[species=cs_p.class_id .== 2]

    return dropdims(sum(juv_groups, dims=:species), dims=:species)
end
function _relative_juveniles(rs::ResultSet)::AbstractArray
    return rs.outcomes[:relative_juveniles]
end
relative_juveniles = Metric(_relative_juveniles, (:timesteps, :sites, :scenarios))


"""
    absolute_juveniles(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    absolute_juveniles(rs::ResultSet)::AbstractArray{<:Real}

Juvenile coral cover in m².
"""
function _absolute_juveniles(X::AbstractArray{T}, area::AbstractVector{T})::AbstractArray{T} where {T<:Real}
    return _relative_juveniles(X) .* area'
end
function _absolute_juveniles(rs::ResultSet)::AbstractArray
    return rs.outcomes[:relative_juveniles] .* rs.site_area'
end
absolute_juveniles = Metric(_absolute_juveniles, (:timesteps, :sites, :scenarios), "m²")


function _max_juvenile_density(max_juv_density=51.8)
    cs_p::DataFrame = coral_spec().params

    max_size_m² = maximum(cs_p[cs_p.class_id.==2, :colony_area_cm2]) ./ 10^4
    return max_juv_density * max_size_m²
end


"""
    juvenile_indicator(X::AbstractArray{T}) where {T<:Real}
    juvenile_indicator(rs::ResultSet)

Indicator for juvenile density (0 - 1), where 1 indicates the maximum theoretical density 
for juveniles have been achieved.

Maximum density is 51.8 juveniles / m², where juveniles are defined as < 5cm diameter.
"""
function _juvenile_indicator(X::AbstractArray{T})::AbstractArray{T} where {T<:Real}
    return _relative_juveniles(X) ./ _max_juvenile_density()
end
function _juvenile_indicator(rs::ResultSet)::AbstractArray
    return rs.outcomes[:relative_juveniles] ./ _max_juvenile_density()
end
juvenile_indicator = Metric(_juvenile_indicator, (:timesteps, :sites, :scenarios))


"""
    coral_evenness(X::AbstractArray{T})::AbstractArray{T} where {T<:Real}
    coral_evenness(rs::ResultSet)::AbstractArray{T} where {T}

Calculates evenness across functional coral groups in ADRIA.
Inverse Simpsons diversity indicator.

# Notes
Number of taxa (distinct groups with enhanced lumped with unenhanced) is hardcoded in this function.

# References
1. Hill, M. O. (1973). 
    Diversity and Evenness: A Unifying Notation and Its Consequences. 
   Ecology, 54(2), 427-432.
   https://doi.org/10.2307/1934352

"""
function _coral_evenness(X::AbstractArray{T})::AbstractArray{T} where {T<:Real}
    x::AbstractArray{<:Real} = min.(max.(X, 0.0), 1.0)
    covers::NamedTuple = coral_cover(x)

    # Evenness as a functional diversity metric
    n::Int64 = 4  # number of taxa
    p1::AbstractArray{<:Real} = dropdims(sum(covers.tab_acr, dims=2), dims=2) ./ covers.relative_cover
    p2::AbstractArray{<:Real} = dropdims(sum(covers.cor_acr, dims=2), dims=2) ./ covers.relative_cover
    p3::AbstractArray{<:Real} = dropdims(sum(covers.small_enc, dims=2), dims=2) ./ covers.relative_cover
    p4::AbstractArray{<:Real} = dropdims(sum(covers.large_mass, dims=2), dims=2) ./ covers.relative_cover

    sum_psqr::AbstractArray{<:Real} = p1 .^ 2 + p2 .^ 2 + p3 .^ 2 + p4 .^ 2  # functional diversity
    simpson_D::AbstractArray{<:Real} = 1 ./ sum_psqr  # Hill 1973, Ecology 54:427-432
    return simpson_D ./ n  # Group evenness
end
function _coral_evenness(rs::ResultSet)::AbstractArray
    return _coral_evenness(rs.raw)
end
coral_evenness = Metric(_coral_evenness, (:timesteps, :sites, :scenarios))


"""
    _colony_Lcm2_to_m3_m2(inputs::DataFrame)::Tuple

Helper function to convert coral colony values from Litres/cm² to m³/m²

# Arguments
- `inputs` : Scenario values for the simulation

# Returns
Tuple : Assumed colony volume (m³/m²) for each species/size class, theoretical maximum for each species/size class

# References
1. Urbina-Barreto, I., Chiroleu, F., Pinel, R., Fréchon, L., Mahamadaly, V., Elise, S., Kulbicki, M., Quod, J.-P., 
     Dutrieux, E., Garnier, R., Henrich Bruggemann, J., Penin, L., & Adjeroud, M. (2021). 
   Quantifying the shelter capacity of coral reefs using photogrammetric 3D modeling: 
     From colonies to reefscapes. 
   Ecological Indicators, 121, 107151. 
   https://doi.org/10.1016/j.ecolind.2020.107151

"""
function _colony_Lcm2_to_m3m2(inputs::Union{DataFrame,DataFrameRow})::Tuple
    _, _, cs_p::DataFrame = coral_spec()
    n_corals::Int64 = length(unique(cs_p.taxa_id))
    n_species::Int64 = length(unique(cs_p.coral_id))

    # Extract assumed colony area (in cm^2) for each taxa/size class from scenario inputs
    # Have to be careful to extract data in the correct order, matching coral id
    colony_area_cm2 = nothing
    try
        colony_area_cm2 = Array{Float64}(inputs[:, cs_p.coral_id.*"_colony_area_cm2"])'
    catch
        # Get from DataFrameRow instead
        colony_area_cm2 = Array{Float64}(inputs[cs_p.coral_id.*"_colony_area_cm2"])
    end

    # Colony planar area parameters (see second column of Table 1 in Urbina-Barreto et al., [1])
    # First column is `b`, second column is `a`
    # log(S) = b + a * log(x)
    pa_params::Array{Float64,2} = Array{Float64,2}([
        -8.32 1.50   # tabular from Urbina-Barretto 2021
        -8.32 1.50   # tabular from Urbina-Barretto 2021
        -7.37 1.34   # columnar from Urbina-Barretto 2021, assumed similar for Corymbose Acropora
        -7.37 1.34   # columnar from Urbina-Barretto 2021, assumed similar for Corymbose Acropora
        -9.69 1.49   # massives from Urbina-Barretto 2021, assumed similar for encrusting and small massives
        -9.69 1.49   # massives from Urbina-Barretto 2021,  assumed similar for large massives
    ])

    log_colony = zeros(n_species)

    # Estimate log colony volume (litres) based on relationship
    # established by Urbina-Barretto 2021, for each taxa/size class and scenario
    # Excuse the dodgy hardcoded indexing

    # Do acropora first ...
    coral_subsets = [i:i+n_corals-1 for i in 1:n_corals:n_species-12]
    for (taxa_id, c_ss) in zip(1:4, coral_subsets)
        c_areas = cs_p[cs_p.taxa_id.==taxa_id, :colony_area_cm2]
        log_colony[c_ss] = pa_params[:, 1] .+ pa_params[:, 2] .* log.(c_areas)
    end

    # ... then massives
    coral_subsets = [i:i+n_corals-1 for i in 25:n_corals:n_species]
    for (taxa_id, c_ss) in zip(5:6, coral_subsets)
        c_areas = cs_p[cs_p.taxa_id.==taxa_id, :colony_area_cm2]
        log_colony[c_ss] = pa_params[:, 1] .+ pa_params[:, 2] .* log.(c_areas)
    end

    # Maximum colony area for each species and scenario, using largest size class
    if ndims(colony_area_cm2) == 1
        max_log_colony = pa_params[6:6:end, 1] .+ pa_params[6:6:end, 2] .* log.(colony_area_cm2[6:6:end])
    else
        max_log_colony = vec(pa_params[6:6:end, 1] .+ pa_params[6:6:end, 2] .* log.(colony_area_cm2[6:6:end, :]))
    end

    colony_litres_per_cm2 = 10.0 .^ log_colony
    max_colony_litres_per_cm2 = 10.0 .^ max_log_colony

    # Convert from dm^3 to m^3
    cm2_m3_per_m2::Float64 = 10^-3
    colony_vol_m3_per_m2::Array{Float64} = colony_litres_per_cm2 * cm2_m3_per_m2
    max_colony_vol_m3_per_m2::Array{Float64} = max_colony_litres_per_cm2 * cm2_m3_per_m2

    return colony_vol_m3_per_m2, max_colony_vol_m3_per_m2
end


"""
    _shelter_species_loop(X, nspecies::Int64, scen::Int64, colony_vol_m3_per_m2, max_colony_vol_m3_per_m2, site_area)

Helper method to calculate relative shelter volume metric across each species/size class for a given scenario.

Note: Species dimension is an amalgamation of taxa and size class.
e.g., X[species=1:6] is Taxa 1, size classes 1-6; X[species=7:12] is Taxa 2, size class 1-6, etc.

# Arguments
- `X` : raw results (proportional coral cover relative to full site area)
- `nspecies` : number of species (taxa and size classes) considered
- `scen` : scenario number to calculate metric for
- `colony_vol_m3_per_m2` : estimated cubic volume per m² of coverage for each species/size class (36)
- `max_colony_vol_m3_per_m2` : theoretical maximum volume per m² of coverage for each taxa (6)
- `site_area` : total area of site in m²
- `k_area` : habitable area of site in m² (i.e., `k` area)
"""
function _shelter_species_loop(X::AbstractArray{T1,3}, nspecies::Int64, colony_vol_m3_per_m2::Array{F}, max_colony_vol_m3_per_m2::Array{F}, site_area::Array{F}, k_area::Array{F}) where {T1<:Real,F<:Float64}
    # Calculate absolute shelter volumes first
    ASV = NamedDimsArray{(:timesteps, :species, :sites)}(zeros(size(X)...))
    _shelter_species_loop!(X, ASV, nspecies, colony_vol_m3_per_m2, site_area)

    MSV::Matrix{Float64} = k_area' .* max_colony_vol_m3_per_m2  # in m³
    # Ensure zero division does not occur
    # ASV should be 0.0 where MSV is 0.0 so the end result is 0.0 / 1.0
    MSV[MSV.==0.0] .= 1.0

    # Loop over each taxa group
    RSV = NamedDimsArray{(:timesteps, :species, :sites)}(zeros(size(X[species=1:6])...))
    taxa_max_map = zip([i:i+5 for i in 1:6:36], 1:6)  # map maximum SV for each group

    # Work out RSV for each taxa
    @inbounds for (sp, sq) in taxa_max_map
        for site in 1:size(ASV, :sites)
            RSV[species=sq, sites=site] .= dropdims(sum(ASV[species=sp, sites=site], dims=:species), dims=:species) ./ MSV[sq, site]
        end
    end

    return RSV
end


"""
    _shelter_species_loop!(X::AbstractArray{T1,3}, ASV::AbstractArray{T1,3}, nspecies::Int64, colony_vol_m3_per_m2, site_area) where {T1}

Helper method to calculate absolute shelter volume metric across each species/size class for a given scenario.

# Arguments
- `X` : raw results (proportional coral cover relative to full site area)
- `ASV` : matrix to hold shelter volume results
- `nspecies` : number of species (taxa and size classes) considered
- `scen` : scenario number to calculate metric for
- `colony_vol_m3_per_m2` : estimated cubic volume per m² of coverage for each species/size class (36)
- `site_area` : area of site in m²
- `k_area` : habitable area of site in m²
"""
function _shelter_species_loop!(X::AbstractArray{T1,3}, ASV::AbstractArray{T1,3}, nspecies::Int64, colony_vol_m3_per_m2, site_area) where {T1<:Real}
    covered_area = nothing

    @inbounds for sp::Int64 in 1:nspecies
        covered_area = X[species=sp] .* site_area'

        # SV represents absolute shelter volume in cubic meters
        ASV[species=sp] .= covered_area .* colony_vol_m3_per_m2[sp]
    end

    clamp!(ASV, 0.0, maximum(ASV))
end


"""
    absolute_shelter_volume(X::NamedDimsArray, site_area::Vector{<:Real}, inputs::Union{DataFrame,DataFrameRow})
    absolute_shelter_volume(rs::ResultSet)

Provide indication of shelter volume in volume of cubic meters.

The metric applies log-log linear models developed by Urbina-Barreto et al., [1]
which uses colony diameter and planar area (2D metrics) to estimate
shelter volume (a 3D metric).


# Arguments
- `X` : raw results
- `site_area` : area in m^2 for each site
- `max_cover` : maximum possible coral cover for each site (in percentage of site_area)
- `inputs` : DataFrame of scenario inputs

# References
1. Urbina-Barreto, I., Chiroleu, F., Pinel, R., Fréchon, L., Mahamadaly, V.,
     Elise, S., Kulbicki, M., Quod, J.-P., Dutrieux, E., Garnier, R.,
     Henrich Bruggemann, J., Penin, L., & Adjeroud, M. (2021).
   Quantifying the shelter capacity of coral reefs using photogrammetric
     3D modeling: From colonies to reefscapes.
   Ecological Indicators, 121, 107151.
   https://doi.org/10.1016/j.ecolind.2020.107151
"""
function _absolute_shelter_volume(X::AbstractArray{T,4}, site_area::Vector{T}, inputs::DataFrame)::AbstractArray{T} where {T<:Real}
    nspecies::Int64 = size(X, :species)

    # Calculate shelter volume of groups and size classes and multiply with area covered
    nscens::Int64 = size(X, :scenarios)
    ASV = NamedDimsArray{(:timesteps, :species, :sites, :scenarios)}(zeros(size(X)...))
    for scen::Int64 in 1:nscens
        colony_vol, _ = _colony_Lcm2_to_m3m2(inputs[scen, :])
        _shelter_species_loop!(X[scenarios=scen], ASV, nspecies, colony_vol, site_area)
    end

    # Sum over groups and size classes to estimate total shelter volume per site
    return dropdims(sum(ASV, dims=:species), dims=:species)
end
function _absolute_shelter_volume(X::AbstractArray{T,3}, site_area::Vector{T}, inputs::Union{DataFrame,DataFrameRow})::AbstractArray{T} where {T<:Real}
    # Collate for a single scenario
    nspecies::Int64 = size(X, :species)

    # Calculate shelter volume of groups and size classes and multiply with area covered
    ASV = NamedDimsArray{(:timesteps, :species, :sites)}(zeros(size(X)...))
    colony_vol, _ = _colony_Lcm2_to_m3m2(inputs)
    _shelter_species_loop!(X, ASV, nspecies, colony_vol, site_area)

    # Sum over groups and size classes to estimate total shelter volume per site
    return dropdims(sum(ASV, dims=:species), dims=:species)
end
function _absolute_shelter_volume(rs::ResultSet)::AbstractArray
    return rs.outcomes[:absolute_shelter_volume]
end
absolute_shelter_volume = Metric(_absolute_shelter_volume, (:timesteps, :sites, :scenarios))


"""
    relative_shelter_volume(X::AbstractArray, site_area::Vector{<:Real}, inputs::DataFrame)
    relative_shelter_volume(rs::ResultSet)

Provide indication of shelter volume relative to theoretical maximum volume for
the area covered by coral.

The metric applies log-log linear models developed by Urbina-Barreto et al., [1]
which uses colony diameter and planar area (2D metrics) to estimate
shelter volume (a 3D metric).

```math
RSV = \\begin{cases}
TASV / MSV & TASV > 0, \\\\
0 & \\text{otherwise}
\\end{cases}
```

where ``TASV`` represents Total Absolute Shelter Volume and ``MSV`` represents the
maximum shelter volume possible.


# Arguments
- `X` : raw results
- `site_area` : area in m^2 for each site
- `inputs` : DataFrame of scenario inputs

# References
1. Urbina-Barreto, I., Chiroleu, F., Pinel, R., Fréchon, L., Mahamadaly, V.,
     Elise, S., Kulbicki, M., Quod, J.-P., Dutrieux, E., Garnier, R.,
     Henrich Bruggemann, J., Penin, L., & Adjeroud, M. (2021).
   Quantifying the shelter capacity of coral reefs using photogrammetric
     3D modeling: From colonies to reefscapes.
   Ecological Indicators, 121, 107151.
   https://doi.org/10.1016/j.ecolind.2020.107151
"""
function _relative_shelter_volume(X::AbstractArray{T,3}, site_area::Vector{T}, k_area::Vector{T}, inputs::Union{DataFrame,DataFrameRow})::AbstractArray{T} where {T<:Real}
    # Collate for a single scenario
    nspecies::Int64 = size(X, :species)

    # Calculate shelter volume of groups and size classes and multiply with covers
    colony_vol, max_colony_vol = _colony_Lcm2_to_m3m2(inputs)
    RSV = _shelter_species_loop(X, nspecies, colony_vol, max_colony_vol, site_area, k_area)

    # @assert !any(RSV .> 1.1)  # Error out in cases where RSV significantly .> 1.0

    # Sum over groups and size classes to estimate total shelter volume
    # proportional to the theoretical maximum (per site)
    RSV = dropdims(sum(RSV, dims=:species), dims=:species)

    clamp!(RSV, 0.0, 1.0)
    return RSV
end
function _relative_shelter_volume(X::AbstractArray{T,4}, site_area::Vector{T}, k_area::Vector{T}, inputs::Union{DataFrame,DataFrameRow})::AbstractArray{T} where {T<:Real}
    @assert nrow(inputs) == size(X, :scenarios)  # Number of results should match number of scenarios

    nspecies::Int64 = size(X, :species)
    nscens::Int64 = size(X, :scenarios)

    # Result template - six entries, one for each taxa
    RSV = NamedDimsArray{(:timesteps, :species, :sites, :scenarios)}(zeros(size(X[:, 1:6, :, :])...))
    for scen::Int64 in 1:nscens
        colony_vol, max_colony_vol = _colony_Lcm2_to_m3m2(inputs[scen, :])
        RSV[scenarios=scen] .= _shelter_species_loop(X[scenarios=scen], nspecies, colony_vol, max_colony_vol, site_area, k_area)
    end

    @assert !any(RSV .> 1.1)  # Error out in cases where RSV significantly .> 1.0

    # Sum over groups and size classes to estimate total shelter volume
    # proportional to the theoretical maximum (per site)
    RSV = dropdims(sum(RSV, dims=:species), dims=:species)

    clamp!(RSV, 0.0, 1.0)
    return RSV
end
function _relative_shelter_volume(rs::ResultSet)::AbstractArray
    return rs.outcomes[:relative_shelter_volume]
end
relative_shelter_volume = Metric(_relative_shelter_volume, (:timesteps, :sites, :scenarios))


"""
    reef_condition_index(TC, E, SV, juveniles)
    reef_condition_index(rs)

Translates coral metrics in ADRIA to a Reef Condition Metrics.

# Notes
Juveniles are made relative to maximum observed juvenile density (51.8/m²)
See email correspondence 
from: Dr. A Thompson; to: Dr. K. Anthony
Subject: RE: Max density of juvenile corals on the GBR
Sent: Friday, 14 October 2022 2:58 PM

# Arguments
- `TC`        : Total relative coral cover across all groups
- `E`         : Evenness across four coral groups
- `SV`        : Shelter volume based coral sizes and abundances
- `juveniles` : Abundance of coral juveniles < 5 cm diameter

Input dimensions: timesteps, species, sites, repeats, scenarios

# Returns
Dimensions: timesteps, sites, repeats, scenarios
"""
function _reef_condition_index(rc::AbstractArray{T}, E::AbstractArray{T}, SV::AbstractArray{T}, juveniles::AbstractArray{T})::AbstractArray{T} where {T<:Real}
    # Compare outputs against reef condition criteria provided by experts

    # These are median values for 7 experts. TODO: draw from distributions
    #  Condition        RC       E       SV      Juv
    # {'VeryGood'}      0.45     0.45    0.45    0.35
    # {'Good'    }      0.35     0.35    0.35    0.25
    # {'Fair'    }      0.25     0.25    0.30    0.25
    # {'Poor'    }      0.15     0.25    0.30    0.25
    # {'VeryPoor'}      0.05     0.15    0.18    0.15

    # Ignoring evenness for now

    # Note that the scores for evenness and juveniles are slightly different
    lin_grid::Gridded{Linear{Throw{OnGrid}}} = Gridded(Linear())
    TC_func::GriddedInterpolation{Float64,1,Vector{Float64},Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.05, 0.15, 0.25, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0], lin_grid)
    # E_func::GriddedInterpolation{Float64, 1, Vector{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.15, 0.25, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.5, 0.7, 0.9, 1.0], lin_grid)
    SV_func::GriddedInterpolation{Float64,1,Vector{Float64},Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.18, 0.30, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.3, 0.5, 0.9, 1.0], lin_grid)
    juv_func::GriddedInterpolation{Float64,1,Vector{Float64},Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.15, 0.25, 0.35, 1.0],), Float64[0, 0.1, 0.5, 0.9, 1.0], lin_grid)

    rc_i::AbstractArray{<:Real} = TC_func.(rc)
    # E_i::T = E_func.(E)
    SV_i::AbstractArray{<:Real} = SV_func.(SV)
    juv_i::AbstractArray{<:Real} = juv_func.(juveniles)

    return mean([rc_i, SV_i, juv_i])
end
function _reef_condition_index(rs::ResultSet)::AbstractArray{<:Real}
    rc::AbstractArray{<:Real} = _relative_cover(rs)

    # Divide across sites by the max possible proportional coral cover
    rc = mapslices((s) -> s ./ (rs.site_max_coral_cover ./ 100.0), rc, dims=2)

    juv::AbstractArray{<:Real} = _juvenile_indicator(rs)
    # E::AbstractArray{<:Real} = _coral_evenness(rs)
    E = Array(Float32[])
    SV::AbstractArray{<:Real} = _relative_shelter_volume(rs)

    return _reef_condition_index(rc, E, SV, juv)
end
reef_condition_index = Metric(_reef_condition_index, (:timesteps, :sites, :scenarios))


include("temporal.jl")
include("site_level.jl")
include("scenario.jl")
include("ranks.jl")
include("pareto.jl")


# """
#     @extend_metric(name, m, args)

# Macro to extend a given metric with additional argument values that are made constant.

# # Arguments
# - name : arbitrary name for generated function
# - m : metric function
# - args : additional arguments whose values will be made constant

# # Example

# ```julia
# import ADRIA.metrics: total_absolute_cover
# extended_metric = @extend_metric(example_func, total_absolute_cover, [site_area(domain)])

# Y = extended_metric(raw_results)  # Equivalent to total_absolute_cover(raw_results, site_area(domain))
# ```
# """
# macro extend_metric(name, m, args)
#     eval(:(($name)(X) = ($m.func)(X, $args...)))
#     return :(Metric(eval($name), $m.dims))
# end


end
