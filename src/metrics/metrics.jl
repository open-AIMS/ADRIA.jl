module metrics

using Interpolations, Statistics, OnlineStats, NamedDims

using DataFrames
import ADRIA: coral_spec, ResultSet


abstract type Outcome end


struct Metric{F<:Function, T<:Tuple} <: Outcome
    func::F
    dims::T
end


function (f::Metric)(raw, args...; kwargs...)
    try
        return f.func(NamedDimsArray{(:timesteps, :species, :sites, :reps, :scenarios)[1:Base.ndims(raw)]}(raw), args...; kwargs...)
    catch
        return f.func(raw, args...; kwargs...)
    end
end


"""
    metric_name(m::Metric)::String

Get name of metric as a string.
"""
function metric_name(m::Metric)::String
    return replace(String(Symbol(m.func)), "_"=>"", count=1)
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
    call_metric(metric, raw, args...; timesteps=(:), species=(:), sites=(:), reps=(:), scens=(:))

Convenience method that slices the data in the specified manner.

# Arguments
- metric : Function, the metric function to apply to "raw" data.
- raw    : NamedDimsArray, raw data to pass into `metric`
- args   : Additional positional arguments to pass into `metric`
- dims   : dummy keyword argument, not used but defined to allow use with other methods
"""
function call_metric(metric::Function, raw::NamedDimsArray, args...; kwargs...)
    dims = haskey(kwargs, :dims) ? kwargs[:dims] : nothing
    if isnothing(dims)
        return metric(slice_results(raw; kwargs...), args...)
    else
        return metric(slice_results(raw; kwargs...), args...; dims=dims)
    end
end


"""
    slice_results(data::NamedDimsArray; timesteps=(:), species=(:), sites=(:), reps=(:), scenarios=(:))

Slice data as indicated.
Dimensions not found in target data are ignored.
"""
function slice_results(data::NamedDimsArray; timesteps=(:), species=(:), sites=(:), reps=(:), scenarios=(:))
    f_dims = (timesteps=timesteps, species=species, sites=sites, reps=reps, scenarios=scenarios)

    s_names = keys(f_dims)
    d_names = dimnames(data)
    common_dims = intersect(s_names, d_names)

    selected_slice = (; zip(common_dims, [getfield(f_dims, k) for k in common_dims])...)
    return data[selected_slice...]
end


"""
    _relative_cover(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    _relative_cover(rs::ResultSet)::AbstractArray{<:Real}

# Arguments
- X : Matrix of raw model results
"""
function _relative_cover(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    # sum over all species and size classes
    return dropdims(sum(X, dims=:species), dims=:species)
end
function _relative_cover(rs::ResultSet)
    return rs.outcomes[:relative_cover]
end


"""
    _total_absolute_cover(X::AbstractArray{<:Real}, site_area::Vector{<:Real})::AbstractArray{<:Real}
    _total_absolute_cover(rs::ResultSet)::AbstractArray{<:Real}

The Total Absolute Coral Cover.
Sum of proportional area taken up by all corals, multiplied by site area.

# Arguments
- X : Matrix of raw model results
- site_area : Vector of site areas, with sites following the same order as given indicated in X.
"""
function _total_absolute_cover(X::AbstractArray{<:Real}, site_area::Vector{<:Real})::AbstractArray{<:Real}
    return _relative_cover(X) .* site_area'
end
function _total_absolute_cover(rs::ResultSet)::AbstractArray{<:Real}
    return _relative_cover(rs) .* rs.site_area'
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
    rc::AbstractArray{<:Real} = _relative_cover(X)  # sum over all species and size classes

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
    juv_all = _juveniles(X)

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
    _juveniles(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    _juveniles(rs::ResultSet)::AbstractArray{<:Real}

Juvenile coral cover in absolute terms.
"""
function _juveniles(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    screen = (x, idx) -> findall(x .== idx)
    _, _, cs_p::DataFrame = coral_spec()

    # Cover of juvenile corals (< 5cm diameter)
    juv_groups::AbstractArray{<:Real} = X[:, screen(cs_p.class_id, 1), :, :, :] .+ X[:, screen(cs_p.class_id, 2), :, :, :]

    return dropdims(sum(juv_groups, dims=:species), dims=:species)
end
function _juveniles(rs::ResultSet)::AbstractArray{<:Real}
    return _juveniles(rs.raw)
end


"""
    _coral_evenness(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    _coral_evenness(rs::ResultSet)::AbstractArray{<:Real}

Calculates evenness across functional coral groups in ADRIA.
Inverse Simpsons diversity indicator.

# Notes
Number of taxa (distinct groups with enhanced lumped with unenhanced) is hardcoded in this function.
"""
function _coral_evenness(X::AbstractArray{<:Real})::AbstractArray{<:Real}
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
function _coral_evenness(rs::ResultSet)::AbstractArray{<:Real}
    return _coral_evenness(rs.raw)
end


"""
    _colony_Lcm2_to_m3_m2(inputs::DataFrame)::Tuple

Helper function to convert coral colony values from Litres/cm² to m³/m²

The metric applies log-log linear models developed by Urbina-Barreto et al., [1]
which uses colony diameter and planar area (2D metrics) to estimate 
shelter volume (a 3D metric).


# Arguments
- X : raw results
- inputs : DataFrame of scenarios

# References
1. Urbina-Barreto, I., Chiroleu, F., Pinel, R., Fréchon, L., Mahamadaly, V., 
     Elise, S., Kulbicki, M., Quod, J.-P., Dutrieux, E., Garnier, R., 
     Henrich Bruggemann, J., Penin, L., & Adjeroud, M. (2021). 
   Quantifying the shelter capacity of coral reefs using photogrammetric 
     3D modeling: From colonies to reefscapes. 
   Ecological Indicators, 121, 107151.
   https://doi.org/10.1016/j.ecolind.2020.107151
"""
function shelter_volume(X::AbstractArray{<:Real}, site_area::Vector{<:Real}, max_cover::Vector{<:Real}, inputs::DataFrame)::AbstractArray{<:Real}
    _, _, cs_p::DataFrame = coral_spec()
    n_corals::Int64 = length(unique(cs_p.taxa_id))

    # Extract assumed colony area (in cm^2) for each taxa/size class from scenario inputs
    colony_area_cm2::Array{Float64} = Array{Float64}(inputs[:, contains.(names(inputs), "colony_area_cm2")])'

    # Colony planar area parameters (see second column of Table 1 in Urbina-Barreto et al., [1])
    pa_params::Array{Float64} = Float64[
        -8.32 1.50;   # tabular from Urbina-Barretto 2021
        -8.32 1.50;   # tabular from Urbina-Barretto 2021
        -7.37 1.34;   # columnar from Urbina-Barretto 2021, assumed similar for Corymbose Acropora
        -7.37 1.34;   # columnar from Urbina-Barretto 2021, assumed similar for Corymbose Acropora
        -9.69 1.49;   # massives from Urbina-Barretto 2021, assumed similar for encrusting and small massives
        -9.69 1.49]   # massives from Urbina-Barretto 2021,  assumed similar for large massives

    # Expand planar area parameters defined above to cover all coral taxa/size classes
    pa_params = repeat(pa_params, n_corals, 1)

    nspecies::Int64 = size(X, :species)
    nscens::Int64 = size(X, :scenarios)

    # Estimate log colony volume (litres) based on relationship
    # established by Urbina-Barretto 2021, for each taxa/size class and scenario
    log_colony = pa_params[:, 1] .+ pa_params[:, 2] .* log10.(colony_area_cm2)
    max_log_colony = pa_params[:, 1] .+ pa_params[:, 2] .* log10.(maximum(colony_area_cm2, dims=1))

    colony_litres_per_cm2 = 10.0 .^ log_colony
    max_colony_litres_per_cm2 = 10.0 .^ max_log_colony

    # Convert from litres per cm2 to m^3 per m^2
    cm2_m3_per_m2::Float32 = (10^-3) * 10^4
    colony_vol_m3_per_m2::Array{Float32} = colony_litres_per_cm2 * cm2_m3_per_m2
    max_colony_vol_m3_per_m2::Array{Float32} = max_colony_litres_per_cm2 * cm2_m3_per_m2

    # calculate shelter volume of groups and size classes and multiply with covers
    sv::NamedDimsArray = NamedDimsArray{(:timesteps, :species, :sites, :reps, :scenarios)}(zeros(size(X)...))
    @inbounds Threads.@threads for sp::Int64 in 1:nspecies
        # Original:
        # sv[:, sp, :, :, :] = (colony_vol_m3_per_m2[sp] / max_colony_vol_m3_per_m2[sp]) .* X[:, sp, :, :, :]

        # Updated:
        # sv_m3 = [m^2 covered by species] * [m^3 volume per m^2]
        # max_sv_m3 = [theoretical max volume per m^2] .* [max possible absolute coral cover area for each site (in m^2)]
        # [SV as proportion of maximum possible SV] = sv_m3 / max_sv_m3
        sv_m3 = (X[species=sp] .* site_area') .* colony_vol_m3_per_m2[sp]
        max_sv_m3 = max_colony_vol_m3_per_m2[sp] .* (site_area .* max_cover)'

        # sv ∈ [0, 1], 0 = no shelter; 1 = maximum shelter
        sv[species=sp] .= sv_m3 ./ max_sv_m3
    end

    # sum over groups and size classes to estimate total proportional shelter volume per site
    return dropdims(sum(sv, dims=:species), dims=:species)
end
function shelter_volume(rs::ResultSet)::AbstractArray{<:Real}
    return shelter_volume(rs.raw, rs.site_area, rs.site_max_coral_cover ./ 100.0, rs.inputs)
end


"""
    _relative_shelter_volume(X::NamedDimsArray, site_area::Vector{<:Real}, inputs::DataFrame)
    _relative_shelter_volume(rs::ResultSet)

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
- X : raw results
- site_area : area in m^2 for each site
- max_cover : maximum possible coral cover for each site (in percentage of site_area)
- inputs : DataFrame of scenario inputs

# References
1. Urbina-Barreto, I., Chiroleu, F., Pinel, R., Fréchon, L., Mahamadaly, V.,
     Elise, S., Kulbicki, M., Quod, J.-P., Dutrieux, E., Garnier, R.,
     Henrich Bruggemann, J., Penin, L., & Adjeroud, M. (2021).
   Quantifying the shelter capacity of coral reefs using photogrammetric
     3D modeling: From colonies to reefscapes.
   Ecological Indicators, 121, 107151.
   https://doi.org/10.1016/j.ecolind.2020.107151
"""
function _relative_shelter_volume(X::NamedDimsArray, site_area::Vector{<:Real}, inputs::DataFrame)::AbstractArray{<:Real}
    nspecies::Int64 = size(X, :species)

    # Calculate shelter volume of groups and size classes and multiply with covers
    if nrow(inputs) > 1
        nscens::Int64 = size(X, :scenarios)
        rsv = NamedDimsArray{(:timesteps, :species, :sites, :reps, :scenarios)}(zeros(size(X)...))
        for scen::Int64 in 1:nscens
            colony_vol, max_colony_vol = _colony_Lcm2_to_m3m2(inputs[scen, :])
            rsv[scenarios=scen] .= _shelter_species_loop(X, nspecies, scen, colony_vol, max_colony_vol, site_area)
        end
    else
        # Collate for a single scenario
        colony_vol, max_colony_vol = _colony_Lcm2_to_m3m2(inputs)
        rsv = _shelter_species_loop(X, nspecies, 1, colony_vol, max_colony_vol, site_area)
        rsv = dropdims(rsv, dims=:scenarios)
    end

    # Sum over groups and size classes to estimate total shelter volume 
    # proportional to the theoretical maximum (per site)
    return dropdims(sum(rsv, dims=:species), dims=:species)
end
function _relative_shelter_volume(rs::ResultSet)::AbstractArray{<:Real}
    return rs.outcomes[:relative_shelter_volume]
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

Input dimensions: timesteps, species, sites, repeats, scenarios

# Returns
Dimensions: timesteps, sites, repeats, scenarios
"""
function _reef_condition_index(rc::AbstractArray{<:Real}, E::AbstractArray{<:Real}, SV::AbstractArray{<:Real}, juveniles::AbstractArray{<:Real})::AbstractArray{<:Real}
    # Compare outputs against reef condition criteria provided by experts

    # These are median values for 7 experts. TODO: draw from distributions
    #  Condition        TC       E       SV      Juv
    # {'VeryGood'}      0.45     0.45    0.45    0.35
    # {'Good'    }      0.35     0.35    0.35    0.25
    # {'Fair'    }      0.25     0.25    0.30    0.25
    # {'Poor'    }      0.15     0.25    0.30    0.25
    # {'VeryPoor'}      0.05     0.15    0.18    0.15

    rc .= min.(rc, 1.0)
    SV .= min.(SV, 1.0)
    juveniles .= min.(juveniles, 1.0)

    # Note that the scores for evenness and juveniles are slightly different
    lin_grid::Gridded{Linear{Throw{OnGrid}}} = Gridded(Linear())
    TC_func::GriddedInterpolation{Float64,1,Float64,Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.05, 0.15, 0.25, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0], lin_grid)
    # E_func::GriddedInterpolation{Float64, 1, Float64, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.15, 0.25, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.5, 0.7, 0.9, 1.0], lin_grid)
    SV_func::GriddedInterpolation{Float64,1,Float64,Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.18, 0.30, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.3, 0.5, 0.9, 1.0], lin_grid)
    juv_func::GriddedInterpolation{Float64,1,Float64,Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.15, 0.25, 0.35, 1.0],), Float64[0, 0.1, 0.5, 0.9, 1.0], lin_grid)

    rc_i::AbstractArray{<:Real} = TC_func.(rc)
    # E_i::T = E_func.(E);
    SV_i::AbstractArray{<:Real} = SV_func.(SV)
    juv_i::AbstractArray{<:Real} = juv_func.(juveniles)

    # Original
    # Y = (rc_i + E_i + SV_i + juv_i) ./ 4;

    # Weighted, giving evenness 10#  weight
    # Y = (rc_i*0.3) + (E_i*0.1) + (SV_i*0.3) + (juv_i*0.3);

    # Removing evenness completely
    # Y = mean([rc_i, SV_i, juv_i])
    # Y = (rc_i .+ SV_i .+ juv_i) ./ 3;

    return mean([rc_i, SV_i, juv_i])
end
function _reef_condition_index(rs::ResultSet)::AbstractArray{<:Real}
    rc::AbstractArray{<:Real} = _relative_cover(rs)

    # Divide across sites by the max possible proportional coral cover
    rc .= mapslices((s) -> s ./ (rs.site_max_coral_cover / 100.0), rc, dims=2)

    juv::AbstractArray{<:Real} = _juveniles(rs)
    # E::AbstractArray{<:Real} = _coral_evenness(rs)
    E = Array(Float32[])
    SV::AbstractArray{<:Real} = _relative_shelter_volume(rs)

    return _reef_condition_index(rc, E, SV, juv)
end


include("temporal.jl")
include("site_level.jl")
include("scenario.jl")


# Wrap base metric functions with dimension metadata
relative_cover = Metric(_relative_cover, (:timesteps, :sites, :reps, :scenarios))
total_absolute_cover = Metric(_total_absolute_cover, (:timesteps, :sites, :reps, :scenarios))
absolute_shelter_volume = Metric(_absolute_shelter_volume, (:timesteps, :sites, :reps, :scenarios))
relative_shelter_volume = Metric(_relative_shelter_volume, (:timesteps, :sites, :reps, :scenarios))
coral_evenness = Metric(_coral_evenness, (:timesteps, :sites, :reps, :scenarios))
juveniles = Metric(_juveniles, (:timesteps, :sites, :reps, :scenarios))
reef_condition_index = Metric(_reef_condition_index, (:timesteps, :sites, :reps, :scenarios))


"""
    @extend_metric(name, m, args)

Macro to extend a given metric with additional argument values that are made constant.

# Arguments
- name : arbitrary name for generated function
- m : metric function
- args : additional arguments whose values will be made constant

# Example

```julia
import ADRIA.metrics: total_absolute_cover
extended_metric = @extend_metric(example_func, total_absolute_cover, [site_area(domain)])

Y = extended_metric(raw_results)  # Equivalent to total_absolute_cover(raw_results, site_area(domain))
```
"""
macro extend_metric(name, m, args)
    eval(:(($name)(X) = ($m.func)(X, $args...)))
    return :(Metric(eval($name), $m.dims))
end


end
