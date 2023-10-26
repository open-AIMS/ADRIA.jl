module metrics

using Distributions, Statistics, Interpolations
using DataFrames, OnlineStats, NamedDims, AxisKeys, JuliennedArrays

using ADRIA: coral_spec, colony_mean_area, ResultSet, timesteps, site_k_area, site_area


abstract type Outcome end


struct Metric{F<:Function,T<:Tuple,S<:String} <: Outcome
    func::F
    dims::T
    unit::S
end
Metric(f, d) = Metric(f, d, "")


"""
    (f::Metric)(raw, args...; kwargs...)
    (f::Metric)(rs::ResultSet, args...; kwargs...)

Makes Metric types callable with arbitary arguments that are passed to associated function.
"""
function (f::Metric)(raw, args...; kwargs...)::NamedDimsArray
    local res

    res = NamedDimsArray{(:timesteps, :species, :sites, :scenarios)[1:Base.ndims(raw)]}(raw)
    res = f.func(res, args...; kwargs...)

    try
        res = NamedDims.rename(res, f.dims[1:Base.ndims(res)])
    catch err
        if !(err isa MethodError)
            rethrow(err)
        end

        res = NamedDimsArray{(f.dims[1:Base.ndims(res)])}(res)
    end

    return res
end
function (f::Metric)(rs::ResultSet, args...; kwargs...)::NamedDimsArray
    met = f.func(rs, args...; kwargs...)
    try
        return NamedDimsArray(met, f.dims)
    catch
        return NamedDims.rename(met, f.dims)
    end
end

"""
    to_string(m::Metric)::String

Get name of metric as a string.
"""
function to_string(m::Metric)::String
    return to_string(m.func)
end
function to_string(f::Function)::String
    return join(split(String(Symbol(f))[2:end], "_"), " ")
end

"""
    to_symbol(m::Metric)::String

Get name of metric as a symbol.
"""
function to_symbol(m::Metric)::Symbol
    return Symbol(replace(to_string(m), ' ' => '_'))
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
    n = titlecase(to_string(f))
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
function call_metric(metric::Union{Function,Metric}, data::NamedDimsArray, args...; kwargs...)
    dims = haskey(kwargs, :dims) ? kwargs[:dims] : nothing
    if isnothing(dims)
        return metric(slice_results(data; kwargs...), args...)
    else
        return metric(slice_results(data; kwargs...), args...; dims=dims)
    end
end


"""
    slice_results(data::NamedDimsArray; timesteps=(:), species=(:), sites=(:), scenarios=(:))

Slice data as indicated.
Dimensions not found in target data are ignored.
"""
function slice_results(data::NamedDimsArray; timesteps=(:), species=(:), sites=(:), scenarios=(:))
    f_dims = (timesteps=timesteps, species=species, sites=sites, scenarios=scenarios)

    s_names = keys(f_dims)
    d_names = NamedDims.dimnames(data)
    common_dims = intersect(s_names, d_names)

    selected_slice = (; zip(common_dims, [getfield(f_dims, k) for k in common_dims])...)
    return data[selected_slice...]
end


function _relative_cover(X::AbstractArray{U}, total_area::Vector{T}, k_area::Vector{T})::AbstractArray{T} where {T<:Real,U<:Real}
    # sum over all species and size classes
    return (dropdims(sum(X, dims=2), dims=2) .* total_area') ./ replace(k_area, 0.0 => 1.0)'
end
function _relative_cover(rs::ResultSet)::AbstractArray
    k_area = site_k_area(rs)'
    rc = copy(rs.outcomes[:total_absolute_cover])
    non_zero_locs = findall(k_area' .> 0.0)
    rc[:, non_zero_locs, :] ./= k_area[non_zero_locs]'

    return rc
end

"""
    relative_cover(X::AbstractArray{T}, k_area::Vector{T})::AbstractArray{T} where {T}
    relative_cover(rs::ResultSet)::AbstractArray

Indicate coral cover relative to available hard substrate (\$k\$ area).

# Arguments
- `X` : Matrix of raw model results

# Returns
Coral cover [0 - 1], relative to available \$k\$ area for a given location.
"""
relative_cover = Metric(_relative_cover, (:timesteps, :sites, :scenarios))


function _total_absolute_cover(X::AbstractArray{T}, k_area::Vector{T})::AbstractArray{T} where {T<:Real}
    return dropdims(sum(X, dims=:species), dims=:species) .* k_area'
end
function _total_absolute_cover(rs::ResultSet)::AbstractArray
    return rs.outcomes[:total_absolute_cover]
end

"""
    total_absolute_cover(X::AbstractArray{T}, k_area::Vector{T})::AbstractArray{T} where {T}
    total_absolute_cover(rs::ResultSet)::AbstractArray{T} where {T}

The Total Absolute Coral Cover.
Sum of proportional area taken up by all corals, multiplied by total site area.

# Arguments
- `X` : Matrix of raw model results
- `k_area` : Vector of site areas, with sites following the same order as given indicated in X.

# Returns
Absolute coral cover for a given location in m².
"""
total_absolute_cover = Metric(_total_absolute_cover, (:timesteps, :sites, :scenarios), "m²")


function _relative_taxa_cover(X::AbstractArray{T}, k_area::Vector{T})::AbstractArray where {T<:Real}
    n_steps, n_species, n_locs = size(X)
    n_sc = 6

    taxa_cover = zeros(n_steps, n_sc)
    k_cover = zeros(n_steps, n_sc, n_locs)
    for (taxa_id, grp) in enumerate([i:i+(n_sc-1) for i in 1:n_sc:n_species])
        for (loc, a) in enumerate(k_area)
            k_cover[:, :, loc] .= X[:, grp, loc] .* a
        end

        # Sum over size class groups
        taxa_cover[:, taxa_id] = vec(sum(k_cover, dims=(2, 3))) ./ sum(k_area)
    end

    return taxa_cover
end
function _relative_taxa_cover(rs::ResultSet)::AbstractArray
    return rs.outcomes[:relative_taxa_cover]
end

"""
    relative_taxa_cover(X::AbstractArray{T}, k_area::Vector{T}) where {T<:Real}
    relative_taxa_cover(rs::ResultSet)

Results grouped by taxa/species.

TODO: Number of size classes is hard coded.

# Arguments
- `X` : Raw model results for a single scenario
- `k_area` : the coral habitable area
- `area` : total location area

# Returns
Coral cover, grouped by taxa for the given scenario, relative to location k area.
"""
relative_taxa_cover = Metric(_relative_taxa_cover, (:timesteps, :taxa, :scenarios))


function _relative_loc_taxa_cover(X::AbstractArray{T}, k_area::Vector{T})::AbstractArray where {T<:Real}
    n_steps, n_species, n_locs = size(X)
    n_sc = 6

    taxa_cover = zeros(n_steps, n_sc, n_locs)
    k_cover = zeros(n_steps, n_sc)
    for (taxa_id, grp) in enumerate([i:i+(n_sc-1) for i in 1:n_sc:n_species])
        @floop for (loc, a) in enumerate(k_area)
            k_cover .= X[:, grp, loc] .* a

            # Sum over size class groups
            taxa_cover[:, taxa_id, loc] = vec(sum(k_cover, dims=2)) ./ a
        end
    end

    return replace!(taxa_cover, NaN => 0.0)
end
# function _relative_loc_taxa_cover(rs::ResultSet)::AbstractArray
#     return rs.outcomes[:relative_loc_taxa_cover]
# end

relative_loc_taxa_cover = Metric(_relative_loc_taxa_cover, (:timesteps, :taxa, :location, :scenarios))


function _relative_juveniles(X::AbstractArray{T}, coral_spec::DataFrame)::AbstractArray{T} where {T<:Real}
    # Cover of juvenile corals (< 5cm diameter)
    juv_groups::AbstractArray{<:Real} = X[species=coral_spec.class_id .== 1] .+ X[species=coral_spec.class_id .== 2]

    return dropdims(sum(juv_groups, dims=:species), dims=:species)
end
function _relative_juveniles(rs::ResultSet)::AbstractArray
    return rs.outcomes[:relative_juveniles]
end

"""
    relative_juveniles(X::AbstractArray{T})::AbstractArray{T} where {T}
    relative_juveniles(rs::ResultSet)::AbstractArray

Juvenile coral cover relative to total site area.
"""
relative_juveniles = Metric(_relative_juveniles, (:timesteps, :sites, :scenarios))


function _absolute_juveniles(X::AbstractArray{T}, coral_spec::DataFrame, area::AbstractVector{T})::AbstractArray{T} where {T<:Real}
    return _relative_juveniles(X, coral_spec) .* area'
end
function _absolute_juveniles(rs::ResultSet)::AbstractArray
    return rs.outcomes[:relative_juveniles] .* site_k_area(rs)'
end

"""
    absolute_juveniles(X::AbstractArray{T}, coral_spec::DataFrame, area::AbstractVector{T})::AbstractArray{T} where {T<:Real}
    absolute_juveniles(rs::ResultSet)::AbstractArray

Juvenile coral cover in m².
"""
absolute_juveniles = Metric(_absolute_juveniles, (:timesteps, :sites, :scenarios), "m²")


"""
    _max_juvenile_area(coral_params::DataFrame, max_juv_density::Float64=51.8)

Calculate the maximum possible area that can be covered by juveniles for a given m².
"""
function _max_juvenile_area(coral_params::DataFrame, max_juv_density::Float64=51.8)
    max_size_m² = maximum(colony_mean_area(coral_params[coral_params.class_id.==2, :mean_colony_diameter_m]))
    return max_juv_density * max_size_m²
end

function _juvenile_indicator(
    X::AbstractArray{T},
    coral_params::DataFrame,
    k_area::V
)::AbstractArray{T} where {T<:Real,V<:Vector{Float64}}
    # Replace 0 k areas with 1.0 to avoid zero-division error
    usable_k_area = Float64[k > 0.0 ? k : 1.0 for k in k_area]'

    return _absolute_juveniles(X, coral_params, k_area) ./ (_max_juvenile_area(coral_params) .* usable_k_area)
end
function _juvenile_indicator(rs::ResultSet)::AbstractArray
    return rs.outcomes[:juvenile_indicator]
end

"""
    juvenile_indicator(X::AbstractArray{T}, coral_params::DataFrame, area::Vector{Float64}, k_area::Vector{Float64}) where {T<:Real}
    juvenile_indicator(rs::ResultSet)

Indicator for juvenile density (0 - 1), where 1 indicates the maximum theoretical density
for juveniles have been achieved.

# Notes
Maximum density is 51.8 juveniles / m², where juveniles are defined as < 5cm diameter.
See email correspondence
from: Dr. A Thompson; to: Dr. K. Anthony
Subject: RE: Max density of juvenile corals on the GBR
Sent: Friday, 14 October 2022 2:58 PM
"""
juvenile_indicator = Metric(_juvenile_indicator, (:timesteps, :sites, :scenarios))

function _coral_evenness(r_taxa_cover::AbstractArray{T})::Array{T} where {T<:Real}
    # Evenness as a functional diversity metric
    n_steps, n_grps, n_locs = size(r_taxa_cover)

    # Sum across groups represents functional diversity
    # Group evenness (Hill 1973, Ecology 54:427-432)
    loc_cover = dropdims(sum(r_taxa_cover, dims=2), dims=2)
    simpsons_diversity = zeros(n_steps, n_locs)
    for loc in axes(loc_cover, 2)
        simpsons_diversity[:, loc] = 1.0 ./ sum((r_taxa_cover[:, :, loc] ./ loc_cover[:, loc]) .^ 2, dims=2)
    end

    return replace!(simpsons_diversity, NaN => 0.0, Inf => 0.0) ./ n_grps
end
function _coral_evenness(rs::ResultSet)::AbstractArray
    return rs.outcomes[:coral_evenness]
end

"""
    coral_evenness(r_taxa_cover::AbstractArray{T})::Array{T} where {T<:Real}
    coral_evenness(rs::ResultSet)::AbstractArray{T} where {T}

Calculates evenness across functional coral groups in ADRIA as a diversity metric.
Inverse Simpsons diversity indicator.

# References
1. Hill, M. O. (1973).
    Diversity and Evenness: A Unifying Notation and Its Consequences.
   Ecology, 54(2), 427-432.
   https://doi.org/10.2307/1934352
"""
coral_evenness = Metric(_coral_evenness, (:timesteps, :sites, :scenarios))

"""
    _colony_Lcm2_to_m3_m2(inputs::DataFrame)::Tuple

Helper function to convert coral colony values from Litres/cm² to m³/m²

# Arguments
- `inputs` : Scenario values for the simulation

# Returns
Tuple : Assumed colony volume (m³/m²) for each species/size class, theoretical maximum for each species/size class

# References
1. Aston Eoghan A., Duce Stephanie, Hoey Andrew S., Ferrari Renata (2022).
    A Protocol for Extracting Structural Metrics From 3D Reconstructions of Corals.
    Frontiers in Marine Science, 9.
    https://doi.org/10.3389/fmars.2022.854395

"""
function _colony_Lcm2_to_m3m2(inputs::NamedDimsArray)::Tuple{Vector{Float64},Vector{Float64}}
    _, _, cs_p::DataFrame = coral_spec()
    n_corals::Int64 = length(unique(cs_p.taxa_id))
    n_species::Int64 = length(unique(cs_p.coral_id))
    n_sizes::Int64 = length(unique(cs_p.class_id))

    # Extract colony diameter (in cm) for each taxa/size class from scenario inputs
    # Have to be careful to extract data in the correct order, matching coral id
    colony_mean_diams_cm::Vector{Float64} = vec(inputs(cs_p.coral_id .* "_mean_colony_diameter_m")) .* 100.0

    # Colony planar area parameters (see Fig 2B in Aston et al., [1])
    # First column is `b`, second column is `a`
    # log(S) = b + a * log(x)
    pa_params::Array{Float64,2} = Array{Float64,2}([
        -8.97 3.14   # Abhorescent Acropora (using branching porites parameters as similar method of growing ever expanding colonies).
        -8.95 2.80   # Tabular Acropora
        -9.13 2.94   # Corymbose Acropora
        -8.90 2.94   # Corymbose non-Acropora (using branching pocillopora values from fig2B)
        -8.87 2.30   # Small massives
        -8.87 2.30   # Large massives
    ])

    # Repeat each entry `n_sizes` times to cover the number size classes represented
    pa_params = repeat(pa_params, inner=(n_sizes, 1))

    # Estimate colony volume (litres) based on relationship
    # established by Aston et al. 2022, for each taxa/size class and scenario
    # Aston et. al. log-log relationship so we apply `exp()` to transform back to dm³
    colony_litres_per_cm2::Vector{Float64} = exp.(pa_params[:, 1] .+ pa_params[:, 2] .* log.(colony_mean_diams_cm))

    # Convert from dm^3 to m^3
    cm2_to_m3_per_m2::Float64 = 10^-3
    colony_vol_m3_per_m2::Vector{Float64} = colony_litres_per_cm2 * cm2_to_m3_per_m2

    # Assumed maximum colony area for each species and scenario, using largest size class
    max_colony_vol_m3_per_m2::Vector{Float64} = colony_vol_m3_per_m2[n_sizes:n_sizes:end]

    return colony_vol_m3_per_m2, max_colony_vol_m3_per_m2
end
function _colony_Lcm2_to_m3m2(inputs::DataFrame)::Tuple
    nd = NamedDimsArray(Matrix(inputs), scenarios=1:nrow(inputs), factors=names(inputs))
    return _colony_Lcm2_to_m3m2(nd)
end


"""
    _shelter_species_loop(X, nspecies::Int64, scen::Int64, colony_vol_m3_per_m2, max_colony_vol_m3_per_m2, site_area)

Helper method to calculate relative shelter volume metric across each species/size class for a given scenario.

Note: Species dimension is an amalgamation of taxa and size class.
e.g., X[species=1:6] is Taxa 1, size classes 1-6; X[species=7:12] is Taxa 2, size class 1-6, etc.

# Arguments
- `X` : raw results (proportional coral cover relative to full site area)
- `n_species` : number of species (taxa and size classes) considered
- `scen` : scenario number to calculate metric for
- `colony_vol_m3_per_m2` : estimated cubic volume per m² of coverage for each species/size class (36)
- `max_colony_vol_m3_per_m2` : theoretical maximum volume per m² of coverage for each taxa (6)
- `site_area` : total area of site in m²
- `k_area` : habitable area of site in m² (i.e., `k` area)
"""
function _shelter_species_loop(X::AbstractArray{T1,3}, n_species::Int64, colony_vol_m3_per_m2::Array{F}, max_colony_vol_m3_per_m2::Array{F}, k_area::Array{F})::NamedDimsArray where {T1<:Real,F<:Float64}
    # Calculate absolute shelter volumes first
    ASV = NamedDimsArray{(:timesteps, :species, :sites)}(zeros(size(X)...))
    _shelter_species_loop!(X, ASV, n_species, colony_vol_m3_per_m2, k_area)

    MSV::Matrix{Float64} = k_area' .* max_colony_vol_m3_per_m2  # in m³
    # Ensure zero division does not occur
    # ASV should be 0.0 where MSV is 0.0 so the end result is 0.0 / 1.0
    MSV[MSV.==0.0] .= 1.0

    # Loop over each taxa group
    RSV = NamedDimsArray{(:timesteps, :species, :sites)}(zeros(size(X[species=1:6])...))
    taxa_max_map = zip([i:i+5 for i in 1:6:n_species], 1:6)  # map maximum SV for each group

    # Work out RSV for each taxa
    for (sp, sq) in taxa_max_map
        @floop for site in 1:size(ASV, :sites)
            RSV[species=sq, sites=site] .= dropdims(sum(ASV[species=sp, sites=site], dims=:species), dims=:species) ./ MSV[sq, site]
        end
    end

    return RSV
end


"""
    _shelter_species_loop!(X::T1, ASV::T1, nspecies::Int64, colony_vol_m3_per_m2::V, k_area::V) where {T1<:NamedDims.NamedDimsArray{(:timesteps, :species, :sites),Float64,3,Array{Float64,3}},V<:AbstractVector{<:Float64}}

Helper method to calculate absolute shelter volume metric across each species/size class for a given scenario.

# Arguments
- `X` : raw results (proportional coral cover relative to full site area)
- `ASV` : matrix to hold shelter volume results
- `nspecies` : number of species (taxa and size classes) considered
- `scen` : scenario number to calculate metric for
- `colony_vol_m3_per_m2` : estimated cubic volume per m² of coverage for each species/size class (36)
- `k_area` : habitable area of site in m²
"""
function _shelter_species_loop!(X::T1, ASV::T1, nspecies::Int64, colony_vol_m3_per_m2::V, k_area::V) where {T1<:NamedDims.NamedDimsArray{(:timesteps, :species, :sites),Float64,3,Array{Float64,3}},V<:AbstractVector{<:Float64}}
    @floop for sp::Int64 in 1:nspecies
        # SV represents absolute shelter volume in cubic meters
        @inbounds ASV[species=sp] .= (X[species=sp] .* k_area') .* colony_vol_m3_per_m2[sp]
    end

    clamp!(ASV, 0.0, maximum(ASV))
end

function _absolute_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::NamedDimsArray)::AbstractArray{T} where {T<:Real}
    nspecies::Int64 = size(X, :species)

    # Calculate shelter volume of groups and size classes and multiply with area covered
    nscens::Int64 = size(X, :scenarios)
    ASV = NamedDimsArray{(:timesteps, :species, :sites, :scenarios)}(zeros(size(X)...))
    for scen::Int64 in 1:nscens
        colony_vol, _ = _colony_Lcm2_to_m3m2(inputs[scen, :])
        _shelter_species_loop!(X[scenarios=scen], ASV, nspecies, colony_vol, k_area)
    end

    # Sum over groups and size classes to estimate total shelter volume per site
    return dropdims(sum(ASV, dims=:species), dims=:species)
end
function _absolute_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::DataFrame)::AbstractArray{T} where {T<:Real}
    ins = NamedDimsArray(Matrix(inputs), scenarios=1:size(inputs, 1), factors=names(inputs))
    return _absolute_shelter_volume(X, k_area, ins)
end
function _absolute_shelter_volume(X::AbstractArray{T,3}, k_area::Vector{T}, inputs::DataFrameRow)::AbstractArray{T} where {T<:Real}
    ins = NamedDimsArray(Matrix(Vector(inputs)'), scenarios=1:1, params=names(inputs))
    return _absolute_shelter_volume(X, k_area, ins)
end
function _absolute_shelter_volume(X::AbstractArray{T,3}, k_area::Vector{T}, inputs::NamedDimsArray)::AbstractArray{T} where {T<:Real}
    # Collate for a single scenario
    nspecies::Int64 = size(X, :species)

    # Calculate shelter volume of groups and size classes and multiply with area covered
    ASV = NamedDimsArray{(:timesteps, :species, :sites)}(zeros(size(X)...))
    colony_vol, _ = _colony_Lcm2_to_m3m2(inputs)
    _shelter_species_loop!(X, ASV, nspecies, colony_vol, k_area)

    # Sum over groups and size classes to estimate total shelter volume per site
    return dropdims(sum(ASV, dims=:species), dims=:species)
end
function _absolute_shelter_volume(rs::ResultSet)::AbstractArray
    return rs.outcomes[:absolute_shelter_volume]
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
absolute_shelter_volume = Metric(_absolute_shelter_volume, (:timesteps, :sites, :scenarios))


function _relative_shelter_volume(X::AbstractArray{T,3}, k_area::Vector{T}, inputs::NamedDimsArray)::AbstractArray{T} where {T<:Real}
    # Collate for a single scenario
    nspecies::Int64 = size(X, :species)

    # Calculate shelter volume of groups and size classes and multiply with covers
    colony_vol::Array{Float64}, max_colony_vol::Array{Float64} = _colony_Lcm2_to_m3m2(inputs)
    RSV::NamedDimsArray = _shelter_species_loop(X, nspecies, colony_vol, max_colony_vol, k_area)

    # @assert !any(RSV .> 1.1)  # Error out in cases where RSV significantly .> 1.0

    # Sum over groups and size classes to estimate total shelter volume
    # proportional to the theoretical maximum (per site)
    RSV = dropdims(sum(RSV, dims=:species), dims=:species)

    clamp!(RSV, 0.0, 1.0)
    return RSV
end
function _relative_shelter_volume(X::AbstractArray{T,3}, k_area::Vector{T}, inputs::DataFrame)::AbstractArray{T} where {T<:Real}
    # Collate for a single scenario
    nscens = size(inputs, 1)
    ins = NamedDimsArray(Matrix(inputs), scenarios=1:nscens, factors=names(inputs))
    return _relative_shelter_volume(X, k_area, ins)
end
function _relative_shelter_volume(X::AbstractArray{T,3}, k_area::Vector{T}, inputs::DataFrameRow)::AbstractArray{T} where {T<:Real}
    # Collate for a single scenario
    ins = NamedDimsArray(Matrix(Vector(inputs)'), scenarios=1, factors=names(inputs))
    return _relative_shelter_volume(X, k_area, ins)
end
function _relative_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::NamedDimsArray)::NamedDimsArray where {T<:Real}
    @assert size(inputs, :scenarios) == size(X, :scenarios)  # Number of results should match number of scenarios

    nspecies::Int64 = size(X, :species)
    nscens::Int64 = size(X, :scenarios)

    # Result template - six entries, one for each taxa
    RSV = NamedDimsArray{(:timesteps, :species, :sites, :scenarios)}(zeros(size(X[:, 1:6, :, :])...))
    for scen::Int64 in 1:nscens
        colony_vol, max_colony_vol = _colony_Lcm2_to_m3m2(inputs[scen, :])
        RSV[scenarios=scen] .= _shelter_species_loop(X[scenarios=scen], nspecies, colony_vol, max_colony_vol, k_area)
    end

    @assert !any(RSV .> 1.1)  # Error out in cases where RSV significantly .> 1.0

    # Sum over groups and size classes to estimate total shelter volume
    # proportional to the theoretical maximum (per site)
    RSV = dropdims(sum(RSV, dims=:species), dims=:species)

    clamp!(RSV, 0.0, 1.0)
    return RSV
end
function _relative_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::DataFrame)::NamedDimsArray where {T<:Real}
    nscens = size(inputs, 1)
    ins = NamedDimsArray(Matrix(inputs), scenarios=1:nscens, factors=names(inputs))
    return _relative_shelter_volume(X, k_area, ins)
end
function _relative_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::DataFrameRow)::NamedDimsArray where {T<:Real}
    ins = NamedDimsArray(Vector(inputs), scenarios=1, factors=names(inputs))
    return _relative_shelter_volume(X, k_area, ins)
end
function _relative_shelter_volume(rs::ResultSet)::NamedDimsArray
    return rs.outcomes[:relative_shelter_volume]
end

"""
    _relative_shelter_volume(X::AbstractArray{T,3}, site_area::Vector{T}, k_area::Vector{T}, inputs::DataFrame)::AbstractArray{T} where {T<:Real}
    _relative_shelter_volume(X::AbstractArray{T,3}, site_area::Vector{T}, k_area::Vector{T}, inputs::DataFrameRow)::AbstractArray{T} where {T<:Real}
    _relative_shelter_volume(X::AbstractArray{T,3}, site_area::Vector{T}, k_area::Vector{T}, inputs::NamedDimsArray)::NamedDimsArray where {T<:Real}
    _relative_shelter_volume(X::AbstractArray{T,4}, site_area::Vector{T}, k_area::Vector{T}, inputs::DataFrame)::NamedDimsArray where {T<:Real}
    _relative_shelter_volume(X::AbstractArray{T,4}, site_area::Vector{T}, k_area::Vector{T}, inputs::DataFrameRow)::NamedDimsArray where {T<:Real}
    _relative_shelter_volume(X::AbstractArray{T,4}, site_area::Vector{T}, k_area::Vector{T}, inputs::NamedDimsArray)::NamedDimsArray where {T<:Real}
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

# Returns
Shelter volume relative to a theoretical maximum volume for the available \$k\$ area.

# References
1. Urbina-Barreto, I., Chiroleu, F., Pinel, R., Fréchon, L., Mahamadaly, V.,
     Elise, S., Kulbicki, M., Quod, J.-P., Dutrieux, E., Garnier, R.,
     Henrich Bruggemann, J., Penin, L., & Adjeroud, M. (2021).
   Quantifying the shelter capacity of coral reefs using photogrammetric
     3D modeling: From colonies to reefscapes.
   Ecological Indicators, 121, 107151.
   https://doi.org/10.1016/j.ecolind.2020.107151
"""
relative_shelter_volume = Metric(_relative_shelter_volume, (:timesteps, :sites, :scenarios))


include("reef_indices.jl")
include("temporal.jl")
include("site_level.jl")
include("scenario.jl")
include("ranks.jl")
include("pareto.jl")


if ccall(:jl_generating_output, Cint, ()) == 1
    Base.precompile(Tuple{Metric{typeof(_relative_shelter_volume),NTuple{4,Symbol},String},Array{Float64,3},Vector{Float64},Vararg{Any}})   # time: 1.6421927
    Base.precompile(Tuple{Metric{typeof(_relative_juveniles),Tuple{Symbol,Symbol,Symbol},String},Array{Float64,3}})   # time: 0.4322211
    Base.precompile(Tuple{Metric{typeof(_absolute_shelter_volume),Tuple{Symbol,Symbol,Symbol},String},Array{Float64,3},Vector{Float64},Vararg{Any}})   # time: 0.3056465
    Base.precompile(Tuple{Metric{typeof(_total_absolute_cover),Tuple{Symbol,Symbol,Symbol},String},Array{Float64,3},Vector{Float64}})   # time: 0.2397704
    Base.precompile(Tuple{Metric{typeof(_relative_taxa_cover),Tuple{Symbol,Symbol,Symbol},String},Array{Float64,3}})   # time: 0.1171238
    Base.precompile(Tuple{typeof(_shelter_species_loop),NamedDimsArray{(:timesteps, :species, :sites),Float64,3,Array{Float64,3}},Int64,Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64}})   # time: 0.0916976
    Base.precompile(Tuple{typeof(_absolute_shelter_volume),NamedDimsArray{(:timesteps, :species, :sites),Float64,3,Array{Float64,3}},Vector{Float64},DataFrame})   # time: 0.0089155
    Base.precompile(Tuple{typeof(_relative_shelter_volume),NamedDimsArray{(:timesteps, :species, :sites),Float64,3,Array{Float64,3}},Vector{Float64},Vector{Float64},DataFrame})   # time: 0.0080667
end

# """
#     @extend_metric(name, m, args)

# Macro to extend a given metric with additional argument values that are made constant.

# # Arguments
# - name : arbitrary name for generated function
# - m : metric function
# - args : additional arguments whose values will be made constant

# # Example

# ```julia
# using ADRIA.metrics: total_absolute_cover
# extended_metric = @extend_metric(example_func, total_absolute_cover, [site_area(domain)])

# Y = extended_metric(raw_results)  # Equivalent to total_absolute_cover(raw_results, site_area(domain))
# ```
# """
# macro extend_metric(name, m, args)
#     eval(:(($name)(X) = ($m.func)(X, $args...)))
#     return :(Metric(eval($name), $m.dims))
# end


end
