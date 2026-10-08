using CSV, DataFrames, LinearAlgebra

"""Area-weighted reef aggregation contract for paired Lizard domain comparisons.

The model state is COTS ha^-1 per site.  Site polygon area (m²) is used only as
the within-reef averaging weight, so the result remains COTS ha^-1.  The
reef-to-site crosswalk is keyed by parent `UNIQUE_ID`; no row-order join is
permitted.  Observation-scale conversion to COTS/tow happens downstream.
"""
function reef_area_weights(
    loc_ids::AbstractVector{<:AbstractString},
    parent_ids::AbstractVector{<:AbstractString},
    areas_m2::AbstractVector{<:Real},
    names_by_parent::AbstractDict{String,String},
    target_reefs::AbstractVector{<:AbstractString},
)
    n = length(loc_ids)
    length(parent_ids) == n && length(areas_m2) == n || error("Site mapping lengths differ")
    length(unique(loc_ids)) == n || error("Site IDs are not unique")
    areas = Float64.(areas_m2)
    all(isfinite, areas) && all(areas .> 0.0) || error("Site areas must be finite and positive")
    all(haskey(names_by_parent, parent) for parent in parent_ids) ||
        error("Site parent missing from reef-name crosswalk")
    reef_names = [names_by_parent[parent] for parent in parent_ids]
    result = Dict{String,NamedTuple{(:indices, :weights, :area_m2),Tuple{Vector{Int},Vector{Float64},Float64}}}()
    for reef in target_reefs
        indices = findall(==(reef), reef_names)
        isempty(indices) && error("Target reef has no sites: $reef")
        total = sum(areas[indices])
        weights = areas[indices] ./ total
        isapprox(sum(weights), 1.0; atol=1e-12) || error("Weights do not sum to one: $reef")
        result[String(reef)] = (indices=indices, weights=weights, area_m2=total)
    end
    return result
end

function load_reef_area_weights(domain, domain_path::String, target_reefs)
    map_path = joinpath(domain_path, "site_to_reef.csv")
    map = CSV.read(map_path, DataFrame)
    (:UNIQUE_ID in propertynames(map) && :reef_name in propertynames(map)) ||
        error("Site-to-reef map lacks UNIQUE_ID or reef_name")
    names_by_parent = Dict{String,String}()
    for row in eachrow(map)
        parent = string(row.UNIQUE_ID)
        reef = first(split(String(row.reef_name), " ("))
        if haskey(names_by_parent, parent) && names_by_parent[parent] != reef
            error("Conflicting reef names for parent $parent")
        end
        names_by_parent[parent] = reef
    end
    loc_ids = String.(domain.loc_ids)
    parents = string.(domain.loc_data.UNIQUE_ID)
    Set(parents) == Set(keys(names_by_parent)) || error("Domain/map parent reef sets differ")
    if :site_id in propertynames(map)
        map_ids = String.(map.site_id)
        length(unique(map_ids)) == length(map_ids) || error("Site-to-reef map IDs are not unique")
        Set(map_ids) == Set(loc_ids) || error("Domain/map site ID sets differ")
        by_site = Dict(map_ids .=> string.(map.UNIQUE_ID))
        all(by_site[id] == parent for (id, parent) in zip(loc_ids, parents)) ||
            error("Domain/map site-to-parent assignment differs")
        if :area_m2 in propertynames(map)
            mapped_area = Dict(map_ids .=> Float64.(map.area_m2))
            all(isapprox(mapped_area[id], area; rtol=1e-9, atol=1e-6)
                for (id, area) in zip(loc_ids, domain.loc_data.area)) ||
                error("Domain/map polygon areas differ")
        end
    end
    return reef_area_weights(
        loc_ids, parents, Float64.(domain.loc_data.area), names_by_parent, target_reefs
    )
end

"""Aggregate a time-by-site matrix of site densities to one reef density series."""
function aggregate_reef_density(matrix::AbstractMatrix{<:Real}, mapping)
    maximum(mapping.indices) <= size(matrix, 2) || error("Site matrix is too narrow")
    return [dot(mapping.weights, @view matrix[t, mapping.indices]) for t in axes(matrix, 1)]
end
