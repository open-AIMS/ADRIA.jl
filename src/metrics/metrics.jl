module metrics

using ADRIAIndicators

using
    Distributions,
    Interpolations,
    JuliennedArrays,
    OnlineStats,
    Statistics

using YAXArrays
using ADRIA:
    DataCube, ZeroDataCube, axes_names, axis_labels, axis_index
using ADRIA: n_sizes, group_indices

using FLoops
using DataFrames

using ADRIA: coral_spec, colony_mean_area, ResultSet, timesteps, loc_k_area, loc_area,
    planar_area_params

abstract type Outcome end

const UNIT_VOLUME = "m³"
const UNIT_AREA = "m²"
const UNIT_AREA_INVERSE = "m⁻²"
const IS_RELATIVE = true
const IS_NOT_RELATIVE = false

struct Metric{F<:Function,T<:Tuple,U<:Tuple,S<:String,B<:Bool} <: Outcome
    func::F
    in_dims::U
    dims::T
    feature::S
    is_relative::B
    unit::S
end

Metric(func, in_dims, dims, feature, is_relative) =
    Metric(func, in_dims, dims, feature, is_relative, "")

"""
    (f::Metric)(raw, args...; kwargs...)
    (f::Metric)(rs::ResultSet, args...; kwargs...)

Makes Metric types callable with arbitary arguments that are passed to associated function.
"""
function (f::Metric)(raw, args...; kwargs...)::YAXArray
    if :scenarios in f.in_dims
        axes = f.in_dims[1:ndims(raw)]
    else
        axes = f.in_dims
    end
    return fill_metadata!(f.func(DataCube(raw, axes), args...; kwargs...), f)
end
function (f::Metric)(rs::ResultSet, args...; kwargs...)::YAXArray
    return fill_metadata!(f.func(rs, args...; kwargs...), f)
end

"""
    relative_cover(X::AbstractArray{<:Real}, loc_area::AbstractVector{<:Real})::AbstractArray{<:Real}
    relative_cover(rs::ResultSet)::AbstractArray{<:Real}

Indicate coral cover relative to available hard substrate (\$k\$ area).

# Arguments
- `X` : Matrix with dimensions (n_timesteps, n_functional_groups * n_size_classes,
n_locations) of raw model results (coral cover relative to available space)

# Returns
Coral cover [0 - 1], relative to available \$k\$ area for the entire study area.
"""
function _relative_cover(
    X::YAXArray{<:Real,4}
)::YAXArray{<:Real}
    dims = (:timesteps, :locations)
    return DataCube(
        ADRIAIndicators.relative_cover(X.data), dims
    )
end
function _relative_cover(X::YAXArray{<:Real,5})
    return dropdims(sum(X; dims=(:groups, :sizes)); dims=(:groups, :sizes))
end
function _relative_cover(rs::ResultSet)::YAXArray{<:Real}
    return rs.outcomes[:relative_cover]
end
relative_cover = Metric(
    _relative_cover,
    (:timesteps, :groups, :sizes, :locations, :scenarios),
    (:timesteps, :locations, :scenarios),
    "Relative Cover",
    IS_RELATIVE
)

"""
    total_absolute_cover(relative_cover::AbstractArray{<:Real}, k_area::Vector{<:Real})::AbstractArray{<:Real}
    total_absolute_cover(rs::ResultSet)::AbstractArray{<:Real}

The Total Absolute Coral Cover.
Sum of proportional area taken up by all corals, multiplied by the location area.

# Arguments
- `relative_cover` : Array with relative_cover
- `k_area` : Proportional area, with locations following the same order as given indicated in `relative_cover`.

# Returns
Absolute coral cover for a given location in $UNIT_AREA.
"""
function _total_absolute_cover(
    relative_cover::YAXArray{<:Real,3},
    k_area::Vector{<:Real}
)::YAXArray
    return DataCube(
        relative_cover.data .* k_area', parentmodule(metrics).axes_names(relative_cover)
    )
end
function _total_absolute_cover(rs::ResultSet)::AbstractArray{<:Real}
    return _total_absolute_cover(relative_cover(rs), loc_k_area(rs))
end
total_absolute_cover = Metric(
    _total_absolute_cover,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :locations, :scenarios),
    "Cover",
    IS_NOT_RELATIVE,
    UNIT_AREA
)

"""
    relative_taxa_cover(X::AbstractArray{<:Real}, k_area::Vector{<:Real}, n_groups::Int64)::AbstractArray{<:Real,2}
    relative_taxa_cover(rs::ResultSet)::AbstractArray{<:Real,2}

Relative coral cover grouped by groups summed up across all locations.

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
    X::AbstractArray{<:Real,4},
    k_area::Vector{<:Real}
)::AbstractArray{<:Real,2}
    n_timesteps = size(X, 1)
    n_groups = size(X, 2)

    taxa_cover::YAXArray = ZeroDataCube(
        (:timesteps, :groups), (n_timesteps, n_groups), X.properties
    )
    ADRIAIndicators.relative_taxa_cover!(X, k_area, taxa_cover.data)

    return taxa_cover
end
function _relative_taxa_cover(rs::ResultSet)::AbstractArray{<:Real,3}
    return rs.outcomes[:relative_taxa_cover]
end
relative_taxa_cover = Metric(
    _relative_taxa_cover,
    (:timesteps, :groups, :sizes, :locations),
    (:timesteps, :groups, :scenarios),
    "Cover",
    IS_RELATIVE
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
    X::AbstractArray{T,4}
)::AbstractArray{T,3} where {T<:Real}
    n_timesteps, n_groups, _, n_locs = size(X)

    taxa_cover::YAXArray = ZeroDataCube(
        (:timesteps, :groups, :locations), (n_timesteps, n_groups, n_locs), X.properties
    )
    ADRIAIndicators.relative_loc_taxa_cover!(X, taxa_cover.data)

    return replace!(taxa_cover, NaN => 0.0)
end
relative_loc_taxa_cover = Metric(
    _relative_loc_taxa_cover,
    (:timesteps, :groups, :sizes, :locations),
    (:timesteps, :groups, :locations, :scenarios),
    "Relative Cover",
    IS_RELATIVE
)

"""
    relative_juveniles(X::AbstractArray{T,3}, coral_spec::DataFrame)::AbstractArray{T,2} where {T<:Real}
    relative_juveniles(rs::ResultSet)::AbstractArray{<:Real,2}

Juvenile coral cover relative to the location's area.

# Arguments
- `X` : Raw model results for a single scenario. Dimensions (n_timesteps, n_group_sizes,
n_locations)
- `coral_spec` : Coral spec DataFrame
"""
function _relative_juveniles(
    X::YAXArray{T,4}
)::YAXArray{T,2} where {T<:Real}
    dims = (:timesteps, :locations)
    is_juvenile::BitVector = falses(size(X, 3))
    is_juvenile[1:2] .= true
    return DataCube(
        ADRIAIndicators.relative_juveniles(X.data, is_juvenile), dims
    )
end
function _relative_juveniles(
    X::YAXArray{T,5}
)::YAXArray{T,2} where {T<:Real}
    dims = (
        :timesteps,
        :locations,
        :scenarios
    )
    is_juvenile::BitVector = falses(size(X, 3))
    is_juvenile[1:2] .= true
    return DataCube(
        ADRIAIndicators.relative_juveniles(X.data, is_juvenile), dims
    )
end
function _relative_juveniles(rs::ResultSet)::AbstractArray{<:Real,3}
    return rs.outcomes[:relative_juveniles]
end
relative_juveniles = Metric(
    _relative_juveniles,
    (:timesteps, :groups, :sizes, :locations, :scenarios),
    (:timesteps, :locations, :scenarios),
    "Relative Cover",
    IS_RELATIVE
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
    X::YAXArray{T}, k_area::AbstractVector{T}
)::YAXArray where {T<:Real}
    dims = ndims(X) == 4 ? (:timesteps, :locations) : (:timesteps, :locations, :scenarios)
    is_juvenile::BitVector = falses(size(X, 3))
    is_juvenile[1:2] .= true
    return DataCube(
        ADRIAIndicators.absolute_juveniles(X.data, is_juvenile, k_area), dims
    )
end
function _absolute_juveniles(rs::ResultSet)::YAXArray
    rel_juv = relative_juveniles(rs)
    return DataCube(rel_juv .* loc_k_area(rs)', parentmodule(metrics).axes_names(rel_juv))
end
absolute_juveniles = Metric(
    _absolute_juveniles,
    (:timesteps, :groups, :sizes, :locations, :scenarios),
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
    juvenile_indicator(X::AbstractArray{T}, coral_spec::DataFrame, k_area::Vector{Float64})::AbstractArray{T,2} where {T<:Real}
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
    X::AbstractArray{T}, coral_spec::DataFrame, k_area::Vector{Float64}
)::AbstractArray{T} where {T<:Real}
    n_groups, n_sizes = size(X)[2:3]
    is_juvenile::BitVector = falses(size(X, 3))
    is_juvenile[1:2] .= true

    mean_diams = permutedims(
        reshape(coral_spec.mean_colony_diameter_m, (n_sizes, n_groups)), (2, 1)
    )
    MAX_DENSITY = 51.8

    # If calculating over multiple scenario include scenarios
    dims = ndims(X) == 4 ? (:timesteps, :locations) : (:timesteps, :locations, :scenarios)
    return DataCube(
        ADRIAIndicators.juvenile_indicator(
            X.data, is_juvenile, k_area, mean_diams, MAX_DENSITY
        ), dims
    )
end
function _juvenile_indicator(rs::ResultSet)::AbstractArray{<:Real,3}
    return rs.outcomes[:juvenile_indicator]
end
juvenile_indicator = Metric(
    _juvenile_indicator,
    (:timesteps, :groups, :sizes, :locations, :scenarios),
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
    n_steps, _, n_locs = size(r_taxa_cover)

    # Sum across groups represents functional diversity
    # Group evenness (Hill 1973, Ecology 54:427-432)
    simpsons_diversity::YAXArray = ZeroDataCube(
        (:timesteps, :locations), (n_steps, n_locs), r_taxa_cover.properties
    )
    ADRIAIndicators.coral_evenness!(r_taxa_cover, simpsons_diversity.data)

    return simpsons_diversity
end
function _coral_evenness(rs::ResultSet)::AbstractArray{<:Real,3}
    return rs.outcomes[:coral_evenness]
end
coral_evenness = Metric(
    _coral_evenness,
    (:timesteps, :groups, :locations),
    (:timesteps, :locations, :scenarios),
    "Evenness Indicator",
    IS_NOT_RELATIVE
)

"""
    coral_diversity(ce::AbstractArray{T})::AbstractArray{T} where {T}
    coral_diversity(rs::ResultSet)::AbstractArray{T} where {T}

Calculates coral diversity metric as the Gini-Simpson index.
This is calculated from coral evenness (which is the inverse Simpson's index, `1/D`)
as `1 - 1/evenness`, which is equivalent to `1 - D`.

# Arguments
- `ce` : Coral evenness (inverse Simpson's index).
- `rs` : A ResultSet object.
"""
function _coral_diversity(
    ce::YAXArray{T}
)::YAXArray{T} where {T<:Real}
    # cd = 1 - (1 / ce)
    # Replace NaNs and Infs with 0.0
    cd = 1.0 .- (1.0 ./ ce)
    replace!(cd.data, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
    return cd
end
function _coral_diversity(rs::ResultSet)::AbstractArray{<:Real}
    ce = coral_evenness(rs)
    return _coral_diversity(ce)
end
coral_diversity = Metric(
    _coral_diversity,
    (:timesteps, :locations, :scenarios), # Input from coral_evenness
    (:timesteps, :locations, :scenarios), # Output dims
    "Diversity",
    IS_NOT_RELATIVE
)

"""
    absolute_shelter_volume(X::YAXArray{T,4}, k_area::Vector{T}, inputs::DataFrameRow)::AbstractArray{T} where {T<:Real}
    absolute_shelter_volume(X::YAXArray{T,4}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
    absolute_shelter_volume(X::YAXArray{T,5}, k_area::Vector{T}, inputs::DataFrame)::AbstractArray{T} where {T<:Real}
    absolute_shelter_volume(X::YAXArray{T,5}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
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
    X::YAXArray{T,4},
    k_area::Vector{T},
    inputs::DataFrameRow
)::AbstractArray{T} where {T<:Real}
    _inputs::YAXArray = DataCube(
        Matrix(Vector(inputs)'); scenarios=1:1, factors=names(inputs)
    )
    return _absolute_shelter_volume(X, k_area, _inputs)
end
function _absolute_shelter_volume(
    X::YAXArray{T,4},
    k_area::Vector{T},
    inputs::YAXArray
)::AbstractArray{T} where {T<:Real}
    # Calculate colony mean area for each group and size class
    n_groups = size(X, 2)
    n_sizes = size(X, 3)
    _, _, cs_p::DataFrame = coral_spec()
    col_mask = inputs.factors .∈ Ref(cs_p.coral_id .* "_mean_colony_diameter_m")
    colony_mean_diams_cm::Array{Float64} =
        reshape(
            (collect(inputs[factors=col_mask]) .* 100.0),
            n_sizes, n_groups
        )'
    col_mean_area = colony_mean_area(colony_mean_diams_cm)

    # Colony planar area parameters (see Fig 2B in Aston et al., [1])
    # First column is `b`, second column is `a`
    # log(S) = b + a * log(x)
    pa_params::Array{Float64,3} = repeat(
        reshape(planar_area_params(), (n_groups, 1, 2)), 1, n_sizes, 1
    )

    ASV::YAXArray = ZeroDataCube(
        (:timesteps, :groups, :sizes, :locations), size(X), X.properties
    )
    ADRIAIndicators.absolute_shelter_volume!(
        X.data, colony_mean_diams_cm, pa_params, k_area, ASV.data
    )

    return dropdims(sum(ASV; dims=(2, 3)); dims=(2, 3))
end
function _absolute_shelter_volume(
    X::YAXArray{T,5},
    k_area::Vector{T},
    inputs::DataFrame
)::AbstractArray{T} where {T<:Real}
    _inputs::YAXArray = DataCube(
        Matrix(inputs); scenarios=1:size(X, 5), factors=names(inputs)
    )
    return _absolute_shelter_volume(X, k_area, _inputs)
end
function _absolute_shelter_volume(
    X::YAXArray{T,5},
    k_area::Vector{T},
    inputs::YAXArray
)::AbstractArray{T} where {T<:Real}
    n_sizes::Int64 = size(X, :sizes)
    n_groups::Int64 = size(X, :groups)
    n_scens::Int64 = size(X, :scenarios)
    pa_params::Array{Float64,3} = repeat(
        reshape(planar_area_params(), (n_groups, 1, 2)), 1, n_sizes, 1
    )

    _, _, cs_p::DataFrame = coral_spec()
    col_mask = inputs.factors .∈ Ref(cs_p.coral_id .* "_mean_colony_diameter_m")
    # Flattened Mean colony diameters with shape [scenarios ⋅ groups_sizes]
    flat_mean_diams::Matrix{Float64} = Matrix(inputs.data[:, col_mask])
    # Mean colony diameters with shape [groups ⋅ sizes ⋅ scenarios]
    colony_mean_diams_cm::Array{Float64,3} =
        permutedims(reshape(
                flat_mean_diams, (n_scens, n_sizes, n_groups)
            ), (3, 2, 1)) .* 100
    col_mean_area = colony_mean_area(colony_mean_diams_cm)

    ASV::YAXArray = ZeroDataCube(
        (:timesteps, :groups, :sizes, :locations, :scenarios), size(X), X.properties
    )
    for scen::Int64 in 1:n_scens
        ADRIAIndicators.absolute_shelter_volume!(
            view(X.data, :, :, :, :, scen),
            view(col_mean_area, :, :, scen),
            pa_params,
            k_area,
            view(ASV.data, :, :, :, :, scen)
        )
    end

    # Sum over groups and size classes to estimate total shelter volume per site
    return dropdims(sum(ASV; dims=(:groups, :sizes)); dims=(:groups, :sizes))
end
function _absolute_shelter_volume(rs::ResultSet)::AbstractArray
    return rs.outcomes[:absolute_shelter_volume]
end
absolute_shelter_volume = Metric(
    _absolute_shelter_volume,
    (:timesteps, :groups, :sizes, :locations, :scenarios),
    (:timesteps, :locations, :scenarios),
    "Volume",
    IS_NOT_RELATIVE,
    UNIT_VOLUME
)

"""
    relative_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::DataFrameRow)::AbstractArray{T} where {T<:Real}
    relative_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
    relative_shelter_volume(X::AbstractArray{T,5}, k_area::Vector{T}, inputs::DataFrame)::AbstractArray{T} where {T<:Real}
    relative_shelter_volume(X::AbstractArray{T,5}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
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
    X::YAXArray{T,4},
    k_area::Vector{T},
    inputs::DataFrameRow
)::AbstractArray{T} where {T<:Real}
    _inputs::YAXArray = DataCube(
        Matrix(Vector(inputs)'); scenarios=1:1, factors=names(inputs)
    )
    return _relative_shelter_volume(X, k_area, _inputs)
end
function _relative_shelter_volume(
    X::AbstractArray{T,4},
    k_area::Vector{T},
    inputs::YAXArray
)::AbstractArray{T} where {T<:Real}
    # Calculate colony mean area for each group and size class
    n_groups = size(X, 2)
    n_sizes = size(X, 3)
    _, _, cs_p::DataFrame = coral_spec()
    col_mask = inputs.factors .∈ Ref(cs_p.coral_id .* "_mean_colony_diameter_m")
    colony_mean_diams_cm::Array{Float64} =
        reshape(
            (collect(inputs[col_mask]) .* 100.0),
            n_sizes, n_groups
        )'
    col_mean_area = colony_mean_area(colony_mean_diams_cm)

    # Colony planar area parameters (see Fig 2B in Aston et al., [1])
    # First column is `b`, second column is `a`
    # log(S) = b + a * log(x)
    pa_params::Array{Float64,3} = repeat(
        reshape(planar_area_params(), (n_groups, 1, 2)), 1, n_sizes, 1
    )

    RSV::YAXArray = ZeroDataCube(
        (:timesteps, :groups, :sizes, :locations), size(X), X.properties
    )
    ADRIAIndicators.relative_shelter_volume!(
        X.data, colony_mean_diams_cm, pa_params, k_area, RSV.data
    )

    RSV = dropdims(sum(RSV; dims=(2, 3)); dims=(2, 3))
    clamp!(RSV, 0.0, 1.0)

    return RSV
end
function _relative_shelter_volume(
    X::YAXArray{T,5},
    k_area::Vector{T},
    inputs::DataFrame
)::AbstractArray{T} where {T<:Real}
    _inputs::YAXArray = DataCube(
        Matrix(inputs); scenarios=1:size(X, 5), factors=names(inputs)
    )
    return _relative_shelter_volume(X, k_area, _inputs)
end
function _relative_shelter_volume(
    X::AbstractArray{T,5},
    k_area::Vector{T},
    inputs::YAXArray
)::YAXArray where {T<:Real}
    n_sizes::Int64 = size(X, :sizes)
    n_groups::Int64 = size(X, :groups)
    n_scens::Int64 = size(X, :scenarios)
    pa_params::Array{Float64,3} = repeat(
        reshape(planar_area_params(), (n_groups, 1, 2)), 1, n_sizes, 1
    )

    _, _, cs_p::DataFrame = coral_spec()
    col_mask = inputs.factors .∈ Ref(cs_p.coral_id .* "_mean_colony_diameter_m")
    # Flattened Mean colony diameters with shape [scenarios ⋅ groups_sizes]
    flat_mean_diams::Matrix{Float64} = Matrix(inputs[:, col_mask])
    # Mean colony diameters with shape [groups ⋅ sizes ⋅ scenarios]
    colony_mean_diams_cm::Array{Float64,3} =
        permutedims(reshape(
                flat_mean_diams, (n_scens, n_sizes, n_groups)
            ), (3, 2, 1)) .* 100

    RSV::YAXArray = ZeroDataCube(
        (:timesteps, :groups, :sizes, :locations, :scenarios), size(X), X.properties
    )
    for scen::Int64 in 1:n_scens
        ADRIAIndicators.relative_shelter_volume!(
            view(X.data, :, :, :, :, scen),
            view(colony_mean_diams_cm, :, :, scen),
            pa_params,
            k_area,
            view(RSV.data, :, :, :, :, scen)
        )
    end

    # Sum over groups and size classes to estimate total shelter volume
    # proportional to the theoretical maximum (per site)
    RSV = dropdims(sum(RSV; dims=(:groups, :sizes)); dims=(:groups, :sizes))

    clamp!(RSV, 0.0, 1.0)
    return RSV
end
function _relative_shelter_volume(rs::ResultSet)::YAXArray
    return rs.outcomes[:relative_shelter_volume]
end
relative_shelter_volume = Metric(
    _relative_shelter_volume,
    (:timesteps, :groups, :sizes, :locations, :scenarios),
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
