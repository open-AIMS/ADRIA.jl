module metrics

using
    Distributions,
    Interpolations,
    JuliennedArrays,
    OnlineStats,
    Statistics

using YAXArrays
using ADRIA:
    DataCube, ZeroDataCube, axes_names, axis_labels, axis_index, default_coral_params,
    default_coral_spec
using ADRIA: n_sizes, group_indices

using FLoops
using DataFrames

using ADRIA: coral_spec, colony_mean_area, ResultSet, timesteps, site_k_area, loc_area,
    planar_area_params

abstract type Outcome end

const UNIT_VOLUME = "m³"
const UNIT_AREA = "m²"
const UNIT_AREA_INVERSE = "m⁻²"
const IS_RELATIVE = true
const IS_NOT_RELATIVE = false

struct Metric{F<:Function,T<:Tuple,S<:String,B<:Bool} <: Outcome
    func::F
    dims::T     # output dimension axes ?
    feature::S
    is_relative::B
    unit::S
end

Metric(func, dims, feature, is_relative) = Metric(func, dims, feature, is_relative, "")

"""
    (f::Metric)(raw, args...; kwargs...)
    (f::Metric)(rs::ResultSet, args...; kwargs...)

Makes Metric types callable with arbitary arguments that are passed to associated function.
"""
function (f::Metric)(raw, args...; kwargs...)::YAXArray
    axes::Tuple = (:timesteps, :species, :locations, :scenarios)[1:ndims(raw)]
    return fill_metadata!(f.func(DataCube(raw, axes), args...; kwargs...), f)
end
function (f::Metric)(rs::ResultSet, args...; kwargs...)::YAXArray
    return fill_metadata!(f.func(rs, args...; kwargs...), f)
end

"""
    relative_cover(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    relative_cover(rs::ResultSet)::AbstractArray{<:Real}

Indicate coral cover relative to available hard substrate (\$k\$ area).

# Arguments
- `X` : Matrix with dimensions (n_timesteps, n_functional_groups * n_size_classes,
n_locations) of raw model results (coral cover relative to available space)

# Returns
Coral cover [0 - 1], relative to available \$k\$ area for a given location.
"""
function _relative_cover(X::YAXArray{<:Real})::YAXArray{<:Real}
    # Sum over all species and size classes
    result::YAXArray = dropdims(sum(X; dims=2); dims=2)
    return result
end
function _relative_cover(rs::ResultSet)::YAXArray{<:Real}
    return rs.outcomes[:relative_cover]
end
relative_cover = Metric(
    _relative_cover, (:timesteps, :locations, :scenarios), "Relative Cover", IS_RELATIVE
)

"""
    total_absolute_cover(X::AbstractArray{<:Real}, k_area::Vector{<:Real})::AbstractArray{<:Real}
    total_absolute_cover(rs::ResultSet)::AbstractArray{<:Real}

The Total Absolute Coral Cover.
Sum of proportional area taken up by all corals, multiplied by total site area.

# Arguments
- `relative_cover` : Array with relative_cover
- `k_area` : Site areas, with sites following the same order as given indicated in X.

# Returns
Absolute coral cover for a given location in $UNIT_AREA.
"""
function _total_absolute_cover(
    relative_cover::AbstractArray{<:Real},
    k_area::Vector{<:Real}
)::AbstractArray{<:Real}
    return relative_cover .* k_area'
end
function _total_absolute_cover(rs::ResultSet)::AbstractArray{<:Real}
    return _total_absolute_cover(rs.outcomes[:relative_cover], site_k_area(rs))
end
total_absolute_cover = Metric(
    _total_absolute_cover,
    (:timesteps, :locations, :scenarios),
    "Cover",
    IS_NOT_RELATIVE,
    UNIT_AREA
)

"""
    relative_taxa_cover(X::AbstractArray{<:Real}, k_area::Vector{<:Real}, n_groups::Int64)::AbstractArray{<:Real,2}
    relative_taxa_cover(rs::ResultSet)::AbstractArray{<:Real,2}

Relative coral cover grouped by taxa/species sumed up across all locations.

# Arguments
- `X` : Raw model results for a single scenario. Dimensions (n_timesteps, n_group_sizes,
n_locations).
- `k_area` : The coral habitable area.
- `n_groups` : Number of function coral groups.

# Returns
Coral cover, grouped by taxa for the given scenario, summed up across all locations,
relative to total k area.
"""
function _relative_taxa_cover(
    X::AbstractArray{<:Real},
    k_area::Vector{<:Real},
    n_groups::Int64
)::AbstractArray{<:Real,2}
    n_timesteps, n_group_sizes, n_locs = size(X)
    _n_sizes::Int64 = n_sizes(n_groups, n_group_sizes)

    taxa_cover::YAXArray = ZeroDataCube(
        (:timesteps, :species), (n_timesteps, n_groups), X.properties
    )
    k_cover = zeros(n_timesteps, _n_sizes, n_locs)
    _group_indices::Vector{UnitRange{Int64}} = group_indices(_n_sizes, n_group_sizes)
    for (idx_group, group) in enumerate(_group_indices)
        # Fill k_cover with absolute cover for each timestep and location and fixed group
        for (idx_loc, loc_k_area) in enumerate(k_area)
            k_cover[:, :, idx_loc] .= X[:, group, idx_loc] .* loc_k_area
        end

        # Sum over size classes and locations  and divide by total k area
        taxa_cover[:, idx_group] = vec(sum(k_cover; dims=(2, 3))) ./ sum(k_area)
    end

    return taxa_cover
end
function _relative_taxa_cover(rs::ResultSet)::AbstractArray{<:Real,3}
    return rs.outcomes[:relative_taxa_cover]
end
relative_taxa_cover = Metric(
    _relative_taxa_cover, (:timesteps, :species, :scenarios), "Cover", IS_RELATIVE
)

"""
    relative_loc_taxa_cover(X::AbstractArray{T}, k_area::Vector{T}, n_groups::Int64)::AbstractArray{T,3} where {T<:Real}

# Arguments
- `X` : Raw model results for a single scenario. Dimensions (n_timesteps, n_group_sizes,
n_locations)
- `k_area` : The coral habitable area.
- `n_groups` : Number of function coral groups.

# Returns
Coral cover, grouped by taxa for the given scenario, for each timestep and location,
relative to location k area.
"""
function _relative_loc_taxa_cover(
    X::AbstractArray{T}, k_area::Vector{T}, n_groups::Int64
)::AbstractArray{T,3} where {T<:Real}
    n_timesteps, n_group_sizes, n_locs = size(X)
    _n_sizes::Int64 = n_sizes(n_groups, n_group_sizes)

    taxa_cover::YAXArray = ZeroDataCube(
        (:timesteps, :species, :locations), (n_timesteps, n_groups, n_locs), X.properties
    )
    k_cover = zeros(n_timesteps, _n_sizes)
    _group_indices::Vector{UnitRange{Int64}} = group_indices(_n_sizes, n_group_sizes)
    for (idx_group, group) in enumerate(_group_indices)
        for (idx_loc, loc_k_area) in enumerate(k_area)
            k_cover .= X[:, group, idx_loc] .* loc_k_area

            # Sum over size class groups
            taxa_cover[:, idx_group, idx_loc] = vec(sum(k_cover; dims=2)) ./ loc_k_area
        end
    end

    return replace!(taxa_cover, NaN => 0.0)
end
relative_loc_taxa_cover = Metric(
    _relative_loc_taxa_cover,
    (:timesteps, :species, :locations, :scenarios),
    "Relative Cover",
    IS_RELATIVE
)

"""
    relative_juveniles(X::AbstractArray{T,3}, coral_spec::DataFrame)::AbstractArray{T,2} where {T<:Real}
    relative_juveniles(rs::ResultSet)::AbstractArray{<:Real,2}

Juvenile coral cover relative to total site area.

# Arguments
- `X` : Raw model results for a single scenario. Dimensions (n_timesteps, n_group_sizes,
n_locations)
- `coral_spec` : Coral spec DataFrame
"""
function _relative_juveniles(
    X::AbstractArray{T,3}, coral_spec::DataFrame
)::AbstractArray{T,2} where {T<:Real}
    # Cover of juvenile corals (< 5cm diameter)
    juv_groups =
        X[species=(coral_spec.class_id .== 1)] .+ X[species=(coral_spec.class_id .== 2)]

    return dropdims(sum(juv_groups; dims=:species); dims=:species)
end
function _relative_juveniles(rs::ResultSet)::AbstractArray{<:Real,3}
    return rs.outcomes[:relative_juveniles]
end
relative_juveniles = Metric(
    _relative_juveniles, (:timesteps, :locations, :scenarios), "Relative Cover", IS_RELATIVE
)

"""
    absolute_juveniles(X::AbstractArray{T,3}, coral_spec::DataFrame, k_area::AbstractVector{T})::AbstractArray{T,2} where {T<:Real}
    absolute_juveniles(rs::ResultSet)::AbstractArray{<:Real,2}

Juvenile coral cover in m².

# Arguments
- `X` : Raw model results for a single scenario. Dimensions (n_timesteps, n_group_sizes,
n_locations)
- `coral_spec` : Coral spec DataFrame
- `k_area` : The coral habitable area.
"""
function _absolute_juveniles(
    X::AbstractArray{T,3}, coral_spec::DataFrame, k_area::AbstractVector{T}
)::AbstractArray{T,2} where {T<:Real}
    return _relative_juveniles(X, coral_spec) .* k_area'
end
function _absolute_juveniles(rs::ResultSet)::AbstractArray{<:Real,3}
    return rs.outcomes[:relative_juveniles] .* site_k_area(rs)'
end
absolute_juveniles = Metric(
    _absolute_juveniles,
    (:timesteps, :locations, :scenarios),
    "Cover",
    IS_NOT_RELATIVE,
    UNIT_AREA
)

"""
    _max_juvenile_area(coral_params::DataFrame, max_juv_density::Float64=51.8)

Calculate the maximum possible area that can be covered by juveniles for a given m².
"""
function _max_juvenile_area(coral_params::DataFrame, max_juv_density::Float64=51.8)
    max_size_m² = maximum(
        colony_mean_area(coral_params[coral_params.class_id .== 2, :mean_colony_diameter_m])
    )
    return max_juv_density * max_size_m²
end

"""
    juvenile_indicator(X::AbstractArray{T,3}, coral_spec::DataFrame, k_area::Vector{Float64})::AbstractArray{T,2} where {T<:Real}
    juvenile_indicator(rs::ResultSet)::AbstractArray{<:Real,2}

Indicator for juvenile density (0 - 1), where 1 indicates the maximum theoretical density
for juveniles have been achieved.

# Arguments
- `X` : Raw model results for a single scenario. Dimensions (n_timesteps, n_group_sizes,
n_locations).
- `coral_spec` : Coral spec DataFrame.
- `k_area` : The coral habitable area.

# Notes
Maximum density is 51.8 juveniles / m², where juveniles are defined as < 5cm diameter.
See email correspondence
from: Dr. A Thompson; to: Dr. K. Anthony
Subject: RE: Max density of juvenile corals on the GBR
Sent: Friday, 14 October 2022 2:58 PM
"""
function _juvenile_indicator(
    X::AbstractArray{T,3}, coral_spec::DataFrame, k_area::Vector{Float64}
)::AbstractArray{T,2} where {T<:Real}
    # Replace 0 k areas with 1.0 to avoid zero-division error
    usable_k_area = Float64[k > 0.0 ? k : 1.0 for k in k_area]'

    return _absolute_juveniles(X, coral_spec, k_area) ./
           (_max_juvenile_area(coral_spec) .* usable_k_area)
end
function _juvenile_indicator(rs::ResultSet)::AbstractArray{<:Real,3}
    return rs.outcomes[:juvenile_indicator]
end
juvenile_indicator = Metric(
    _juvenile_indicator,
    (:timesteps, :locations, :scenarios),
    "Density Indicator",
    IS_NOT_RELATIVE,
    UNIT_AREA_INVERSE)

"""
    coral_evenness(r_taxa_cover::AbstractArray{T})::AbstractArray{T} where {T<:Real}
    coral_evenness(rs::ResultSet)::AbstractArray{T} where {T}

Calculates evenness across functional coral groups in ADRIA as a diversity metric.
Inverse Simpsons diversity indicator.

# References
1. Hill, M. O. (1973).
Diversity and Evenness: A Unifying Notation and Its Consequences.
Ecology, 54(2), 427-432.
https://doi.org/10.2307/1934352
"""
function _coral_evenness(
    r_taxa_cover::AbstractArray{T,3}
)::AbstractArray{T,2} where {T<:Real}
    # Evenness as a functional diversity metric
    n_steps, n_grps, n_locs = size(r_taxa_cover)

    # Sum across groups represents functional diversity
    # Group evenness (Hill 1973, Ecology 54:427-432)
    simpsons_diversity::YAXArray = ZeroDataCube(
        (:timesteps, :locations), (n_steps, n_locs), r_taxa_cover.properties
    )
    loc_cover = dropdims(sum(r_taxa_cover; dims=2); dims=2)
    for loc in axes(loc_cover, 2)
        simpsons_diversity[:, loc] =
            1.0 ./ sum((r_taxa_cover[:, :, loc] ./ loc_cover[:, loc]) .^ 2; dims=2)
    end

    return replace!(
        simpsons_diversity, NaN => 0.0, Inf => 0.0
    ) ./ n_grps
end
function _coral_evenness(rs::ResultSet)::AbstractArray{<:Real,3}
    return rs.outcomes[:coral_evenness]
end
coral_evenness = Metric(
    _coral_evenness,
    (:timesteps, :locations, :scenarios),
    "Evenness Indicator",
    IS_NOT_RELATIVE
)

"""
    _colony_Lcm2_to_m3m2(inputs::DataFrame)::Tuple
    _colony_Lcm2_to_m3m2(inputs::YAXArray)::Tuple{Vector{Float64},Vector{Float64}}

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
function _colony_Lcm2_to_m3m2(inputs::DataFrame, coral_params::Dict)::Tuple
    _inputs = DataCube(Matrix(inputs); scenarios=1:nrow(inputs), factors=names(inputs))
    return _colony_Lcm2_to_m3m2(_inputs, coral_params)
end
function _colony_Lcm2_to_m3m2(
    inputs::YAXArray, coral_params::Dict
)::Tuple{Vector{Float64},Vector{Float64}}
    _, _, cs_p::DataFrame = coral_spec(coral_params)
    n_sizes::Int64 = length(unique(cs_p.class_id))

    # Extract colony diameter (in cm) for each taxa/size class from scenario inputs
    # Have to be careful to extract data in the correct order, matching coral id
    n_groups_sizes::Int64 = size(cs_p, 1)

    colony_mean_diams_cm::Vector{Float64} = reshape(
        (inputs[factors=At(cs_p.coral_id .* "_mean_colony_diameter_m")] .* 100.0).data,
        n_groups_sizes
    )

    # Colony planar area parameters (see Fig 2B in Aston et al., [1])
    # First column is `b`, second column is `a`
    # log(S) = b + a * log(x)
    pa_params::Array{Float64,2} = coral_params[:planar_area_params]

    # Repeat each entry `n_sizes` times to cover the number size classes represented
    pa_params = repeat(pa_params; inner=(n_sizes, 1))

    # Estimate colony volume (litres) based on relationship
    # established by Aston et al. 2022, for each taxa/size class and scenario
    # Aston et. al. log-log relationship so we apply `exp()` to transform back to dm³
    colony_litres_per_cm2::Vector{Float64} =
        exp.(pa_params[:, 1] .+ pa_params[:, 2] .* log.(colony_mean_diams_cm))

    # Convert from dm^3 to m^3
    cm2_to_m3_per_m2::Float64 = 10^-3
    colony_vol_m3_per_m2::Vector{Float64} = colony_litres_per_cm2 * cm2_to_m3_per_m2

    # Assumed maximum colony area for each species and scenario, using largest size class
    max_colony_vol_m3_per_m2::Vector{Float64} = colony_vol_m3_per_m2[n_sizes:n_sizes:end]

    return colony_vol_m3_per_m2, max_colony_vol_m3_per_m2
end

"""
    _shelter_species_loop(X::AbstractArray{T1,3}, n_species::Int64, colony_vol_m3_per_m2::Array{F}, max_colony_vol_m3_per_m2::Array{F}, k_area::Array{F})::YAXArray where {T1<:Real,F<:Float64}

Helper method to calculate relative shelter volume metric across each species/size class for a given scenario.

Note: Species dimension is an amalgamation of taxa and size class.
e.g., X[species=1:6] is Taxa 1, size classes 1-6; X[species=7:12] is Taxa 2, size class 1-6, etc.

# Arguments
- `X` : raw results (proportional coral cover relative to full site area)
- `n_group_and_size` : number of species (taxa and size classes) considered
- `colony_vol_m3_per_m2` : estimated cubic volume per m² of coverage for each species/size class
- `max_colony_vol_m3_per_m2` : theoretical maximum volume per m² of coverage for each taxa
- `k_area` : habitable area of site in m² (i.e., `k` area)
"""
function _shelter_species_loop(
    X::AbstractArray{T1,3},
    n_group_and_size::Int64,
    colony_vol_m3_per_m2::Array{F},
    max_colony_vol_m3_per_m2::Array{F},
    k_area::Array{F}
)::YAXArray where {T1<:Real,F<:Float64}
    # Calculate absolute shelter volumes first
    ASV::YAXArray = ZeroDataCube((:timesteps, :species, :locations), size(X))

    _shelter_species_loop!(X, ASV, n_group_and_size, colony_vol_m3_per_m2, k_area)

    # Maximum shelter volume
    MSV::Matrix{Float64} = k_area' .* max_colony_vol_m3_per_m2  # in m³
    # Ensure zero division does not occur
    # ASV should be 0.0 where MSV is 0.0 so the end result is 0.0 / 1.0
    MSV[MSV .== 0.0] .= 1.0

    # Number of functional groups
    n_groups::Int64 = size(MSV, 1)
    # Number of size classes
    n_sizes::Int64 = Int64(n_group_and_size / n_groups)
    # Loop over each taxa group

    RSV::YAXArray = ZeroDataCube(
        (:timesteps, :species, :locations), size(X[species=1:n_groups]), X.properties
    )
    taxa_max_map = zip(
        [i:(i + n_sizes - 1) for i in 1:n_sizes:n_group_and_size], 1:n_groups
    )  # map maximum SV for each group

    # Work out RSV for each taxa
    for (sp, sq) in taxa_max_map
        for site in 1:size(ASV, :locations)
            RSV[species=At(sq), locations=At(site)] .=
                dropdims(
                    sum(ASV[species=At(sp), locations=At(site)]; dims=:species);
                    dims=:species
                ) ./ MSV[sq, site]
        end
    end

    return RSV
end

"""
_shelter_species_loop!(X::YAXArray, ASV::YAXArray, nspecies::Int64, colony_vol_m3_per_m2::V, k_area::V) where {V<:AbstractVector{<:Float64}}

Helper method to calculate absolute shelter volume metric across each species/size class for a given scenario.

# Arguments
- `X` : raw results (proportional coral cover relative to full site area)
- `ASV` : matrix to hold shelter volume results
- `nspecies` : number of species (taxa and size classes) considered
- `scen` : scenario number to calculate metric for
- `colony_vol_m3_per_m2` : estimated cubic volume per m² of coverage for each species/size class
- `k_area` : habitable area of site in m²
"""
function _shelter_species_loop!(
    X::YAXArray,
    ASV::YAXArray,
    nspecies::Int64,
    colony_vol_m3_per_m2::V,
    k_area::V
) where {V<:AbstractVector{<:Float64}}
    for sp::Int64 in 1:nspecies
        # SV represents absolute shelter volume in cubic meters
        ASV[species=At(sp)] = (X[species=At(sp)] .* k_area') .* colony_vol_m3_per_m2[sp]
    end
end

"""
    absolute_shelter_volume(X::YAXArray{T,3}, k_area::Vector{T}, inputs::DataFrameRow)::AbstractArray{T} where {T<:Real}
    absolute_shelter_volume(X::YAXArray{T,4}, k_area::Vector{T}, inputs::DataFrame)::AbstractArray{T} where {T<:Real}
    absolute_shelter_volume(X::YAXArray{T,3}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
    absolute_shelter_volume(X::YAXArray{T,4}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
    absolute_shelter_volume(rs::ResultSet)

Provide indication of shelter volume in volume of cubic meters.

The metric applies log-log linear models developed by Urbina-Barreto et al., [1]
which uses colony diameter and planar area (2D metrics) to estimate
shelter volume (a 3D metric).


# Arguments
- `X` : raw results
- `k_area` : area in m^2 for each site
- `max_cover` : maximum possible coral cover for each site (in percentage of loc_area)
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
function _absolute_shelter_volume(
    X::YAXArray{T,3},
    k_area::Vector{T},
    inputs::DataFrameRow;
    coral_params=default_coral_params()
)::AbstractArray{T} where {T<:Real}
    _inputs::YAXArray = DataCube(
        Matrix(Vector(inputs)'); scenarios=1:1, factors=names(inputs)
    )

    return _absolute_shelter_volume(X, k_area, _inputs; coral_params=coral_params)
end
function _absolute_shelter_volume(
    X::YAXArray{T,4},
    k_area::Vector{T},
    inputs::DataFrame;
    coral_params=default_coral_params()
)::AbstractArray{T} where {T<:Real}
    _inputs::YAXArray = DataCube(
        Matrix(inputs); scenarios=1:size(inputs, 1), factors=names(inputs)
    )

    return _absolute_shelter_volume(X, k_area, _inputs; coral_params=coral_params)
end
function _absolute_shelter_volume(
    X::YAXArray{T,3},
    k_area::Vector{T},
    inputs::YAXArray;
    coral_params=default_coral_params()
)::AbstractArray{T} where {T<:Real}
    # Collate for a single scenario
    nspecies::Int64 = size(X, :species)

    # Calculate shelter volume of groups and size classes and multiply with area covered
    ASV::YAXArray = ZeroDataCube(
        (:timesteps, :species, :locations), size(X), X.properties
    )

    colony_vol, _ = _colony_Lcm2_to_m3m2(inputs, coral_params)
    _shelter_species_loop!(X, ASV, nspecies, colony_vol, k_area)

    # Sum over groups and size classes to estimate total shelter volume per site
    return dropdims(sum(ASV; dims=:species); dims=:species)
end
function _absolute_shelter_volume(
    X::YAXArray{T,4},
    k_area::Vector{T},
    inputs::YAXArray;
    coral_params=default_coral_params()
)::AbstractArray{T} where {T<:Real}
    nspecies::Int64 = size(X, :species)

    # Calculate shelter volume of groups and size classes and multiply with area covered
    nscens::Int64 = size(X, :scenarios)
    ASV::YAXArray = ZeroDataCube(
        (:timesteps, :species, :locations, :scenarios), size(X), X.properties
    )

    for scen::Int64 in 1:nscens
        colony_vol, _ = _colony_Lcm2_to_m3m2(inputs[scen, :], coral_params)
        _shelter_species_loop!(
            X[scenarios=scen], ASV[scenarios=scen], nspecies, colony_vol, k_area
        )
    end

    # Sum over groups and size classes to estimate total shelter volume per site
    return dropdims(sum(ASV; dims=:species); dims=:species)
end
function _absolute_shelter_volume(rs::ResultSet)::AbstractArray
    return rs.outcomes[:absolute_shelter_volume]
end
absolute_shelter_volume = Metric(
    _absolute_shelter_volume,
    (:timesteps, :locations, :scenarios),
    "Volume",
    IS_NOT_RELATIVE,
    UNIT_VOLUME
)

"""
    relative_shelter_volume(X::AbstractArray{T,3}, k_area::Vector{T}, inputs::DataFrameRow)::AbstractArray{T} where {T<:Real}
    relative_shelter_volume(X::AbstractArray{T,3}, k_area::Vector{T}, inputs::DataFrame)::AbstractArray{T} where {T<:Real}
    relative_shelter_volume(X::AbstractArray{T,3}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
    relative_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::DataFrameRow)::AbstractArray{T} where {T<:Real}
    relative_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::DataFrame)::AbstractArray{T} where {T<:Real}
    relative_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
    relative_shelter_volume(rs::ResultSet)

Provide indication of shelter volume relative to theoretical maximum volume for the area
covered by coral.

The metric applies log-log linear models developed by Urbina-Barreto et al., [1] which uses
colony diameter and planar area (2D metrics) to estimate shelter volume (a 3D metric).

```math
RSV = \\begin{cases}
TASV / MSV & MSV > 0, \\\\
0 & \\text{otherwise}
\\end{cases}
```

where ``TASV`` represents Total Absolute Shelter Volume and ``MSV`` represents the
maximum shelter volume possible.

# Arguments
- `X` : raw results
- `k_area` : area in m^2 for each site
- `scens` : DataFrame of scenario inputs

# Returns
Shelter volume relative to a theoretical maximum volume for the available \$k\$ area.

# References
1. Urbina-Barreto, I., Chiroleu, F., Pinel, R., Fréchon, L., Mahamadaly, V., Elise, S.,
Kulbicki, M., Quod, J.-P., Dutrieux, E., Garnier, R., Henrich Bruggemann, J., Penin, L.,
& Adjeroud, M. (2021). Quantifying the shelter capacity of coral reefs using photogrammetric
3D modeling: From colonies to reefscapes. Ecological Indicators, 121, 107151.
https://doi.org/10.1016/j.ecolind.2020.107151
"""
function _relative_shelter_volume(
    X::AbstractArray{T,3},
    k_area::Vector{T},
    scens::DataFrameRow;
    coral_params=default_coral_params()
)::AbstractArray{T} where {T<:Real}
    return _relative_shelter_volume(X, k_area, DataFrame(scens); coral_params=coral_params)
end
function _relative_shelter_volume(
    X::AbstractArray{T,3},
    k_area::Vector{T},
    scens::DataFrame;
    coral_params=default_coral_params()
)::AbstractArray{T} where {T<:Real}
    @assert size(scens, 1) == 1 "Scens DataFrame should have only one line"
    _inputs::YAXArray = DataCube(
        Matrix(scens); scenarios=axes(scens, 1), factors=names(scens)
    )
    return _relative_shelter_volume(X, k_area, _inputs; coral_params=coral_params)
end
function _relative_shelter_volume(
    X::AbstractArray{T,3},
    k_area::Vector{T},
    scens::YAXArray;
    coral_params=default_coral_params()
)::AbstractArray{T} where {T<:Real}
    # Collate for a single scenario
    nspecies::Int64 = size(X, :species)

    # Calculate shelter volume of groups and size classes and multiply with covers
    colony_vol::Array{Float64}, max_colony_vol::Array{Float64} = _colony_Lcm2_to_m3m2(
        scens, coral_params
    )
    RSV::YAXArray = _shelter_species_loop(X, nspecies, colony_vol, max_colony_vol, k_area)

    # Sum over groups and size classes to estimate total shelter volume
    # proportional to the theoretical maximum (per site)
    RSV = dropdims(sum(RSV; dims=:species); dims=:species)

    clamp!(RSV, 0.0, 1.0)
    return RSV
end
function _relative_shelter_volume(
    X::AbstractArray{T,4},
    k_area::Vector{T},
    scens::DataFrameRow;
    coral_params=default_coral_params()
)::AbstractArray{T} where {T<:Real}
    return _relative_shelter_volume(X, k_area, DataFrame(scens); coral_params=coral_params)
end
function _relative_shelter_volume(
    X::AbstractArray{T,4},
    k_area::Vector{T},
    scens::DataFrame;
    coral_params=default_coral_params()
)::AbstractArray{T} where {T<:Real}
    _inputs::YAXArray = DataCube(
        Matrix(scens); scenarios=axes(scens, 1), factors=names(scens)
    )
    return _relative_shelter_volume(X, k_area, _inputs; coral_params=coral_params)
end
function _relative_shelter_volume(
    X::AbstractArray{T,4},
    k_area::Vector{T},
    scens::YAXArray;
    coral_params=default_coral_params()
)::YAXArray where {T<:Real}
    @assert size(scens, :scenarios) == size(X, :scenarios)  # Number of results should match number of scenarios

    nspecies::Int64 = size(X, :species)
    nscens::Int64 = size(X, :scenarios)

    # Result template - six entries, one for each taxa
    n_groups::Int64 = length(default_coral_spec().taxa_names)
    RSV::YAXArray = ZeroDataCube(
        (:timesteps, :species, :locations, :scenarios),
        size(X[:, 1:n_groups, :, :]),
        X.properties
    )
    for scen::Int64 in 1:nscens
        colony_vol, max_colony_vol = _colony_Lcm2_to_m3m2(scens[scen, :], coral_params)
        RSV[scenarios=scen] .= _shelter_species_loop(
            X[scenarios=scen], nspecies, colony_vol, max_colony_vol, k_area
        )
    end

    # Sum over groups and size classes to estimate total shelter volume
    # proportional to the theoretical maximum (per site)
    RSV = dropdims(sum(RSV; dims=:species); dims=:species)

    clamp!(RSV, 0.0, 1.0)
    return RSV
end
function _relative_shelter_volume(rs::ResultSet)::YAXArray
    return rs.outcomes[:relative_shelter_volume]
end
relative_shelter_volume = Metric(
    _relative_shelter_volume,
    (:timesteps, :locations, :scenarios),
    "Relative Volume",
    IS_RELATIVE
)

include("metadata.jl")
include("pareto.jl")
include("ranks.jl")
include("reef_indices.jl")
include("scenario.jl")
include("spatial.jl")
include("temporal.jl")
include("utils.jl")

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
# extended_metric = @extend_metric(example_func, total_absolute_cover, [loc_area(domain)])

# Y = extended_metric(raw_results)  # Equivalent to total_absolute_cover(raw_results, loc_area(domain))
# ```
# """
# macro extend_metric(name, m, args)
#     eval(:(($name)(X) = ($m.func)(X, $args...)))
#     return :(Metric(eval($name), $m.dims))
# end

end
