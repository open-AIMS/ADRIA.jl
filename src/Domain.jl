abstract type Domain end

"""
    EnvLayer{S, TF}

Store environmental data layers used for scenario
"""
mutable struct EnvLayer{S<:AbstractString,TF}
    dpkg_path::S
    loc_data_fn::S
    const loc_id_col::S
    const cluster_id_col::S
    init_coral_cov_fn::S
    connectivity_fn::S
    DHW_fn::S
    wave_fn::S
    const timeframe::TF
end

"""
    load_domain(path::String)

Load ADRIA domain specification from data package.
No SSP/RCP data is preset.

# Arguments
- `path` : location of data package.
"""


function unique_loc_ids(d::Domain)::Vector{String}
    return d.loc_data[:, d.loc_id_col]
end

"""
    distance_matrix(coords::Vector{Tuple{Float64,Float64}})::Matrix{Float64}
    distance_matrix(dom::Domain)::Matrix{Float64}

Calculate the pairwise distance matrix for a set of location coordinates.

# Arguments
- `coords`: Vector of coordinate tuples (longitude, latitude)

# Returns
Matrix of pairwise distances between locations
"""
function distance_matrix(coords::Vector{Tuple{Float64,Float64}})::Matrix{Float64}
    n_locs = length(coords)
    dist_matrix = zeros(n_locs, n_locs)

    for i in 1:n_locs
        for j in (i + 1):n_locs
            dist = Distances.haversine(coords[i], coords[j])
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
        end
    end

    return dist_matrix
end

"""
    distance_matrix(dom::Domain)::Matrix{Float64}

Convenience method to extract distance matrix.

# Arguments
- `dom`: Domain holding geospatial location data.

# Returns
Matrix of pairwise distances between locations
"""
function distance_matrix(dom::Domain)::Matrix{Float64}
    return distance_matrix(centroids(dom))
end
function distance_matrix(loc_data::DataFrame)::Matrix{Float64}
    return distance_matrix(centroids(loc_data))
end

"""
    mean_distance(dist_matrix::Matrix{Float64})::Vector{Float64}

Calculate the mean pairwise distance between each location and all other locations to
give an indication of spatial spread.

# See also
- `nearest_neighbor_distances()`
"""
function mean_distance(dist_matrix::Matrix{Float64})::Vector{Float64}
    n_locs = size(dist_matrix, 1)
    mean_distances = vec(sum(dist_matrix; dims=2) ./ (n_locs - 1))

    return mean_distances
end

"""
    nearest_neighbor_distances(dist_matrix::Matrix{Float64}, n_neighbors::Int64)::Vector{Float64}

Calculate the mean distance to the n closest neighbors for each location.

# Arguments
- `dist_matrix`: Pairwise distance matrix between locations
- `n_neighbors`: Number of closest neighbors to consider

# Returns
Vector of mean distances (in meters) to nearest `n_neighbors` for each location.
"""
function nearest_neighbor_distances(
    dist_matrix::Matrix{Float64},
    n_neighbors::Int
)::Vector{Float64}
    n_locs = size(dist_matrix, 1)

    # Adjust n_neighbors if there aren't enough locations
    effective_n = min(n_neighbors, n_locs - 1)

    # Calculate mean distance to n closest neighbors for each location
    mean_distances = zeros(n_locs)
    for i in 1:n_locs
        # Get distances from this location to all others, excluding self
        other_indices = [j for j in 1:n_locs if j != i]
        distances = dist_matrix[i, other_indices]

        # Find the n closest neighbors
        if length(distances) <= effective_n
            # If we have few locations, use all available
            mean_distances[i] = mean(distances)
        else
            # Sort and take the effective_n closest
            sorted_indices = sortperm(distances)[1:effective_n]
            mean_distances[i] = mean(distances[sorted_indices])
        end
    end

    return mean_distances
end

"""
    nearest_neighbor_distances(dom::Domain, n_neighbors::Int64)::Vector{Float64}

Calculate the mean distance to the `n_neighbors` closest neighbors for each location.

# Arguments
- `dist_matrix`: Pairwise distance matrix between locations
- `n_neighbors`: Number of closest neighbors to consider

# Returns
Vector of mean distances (in m) to nearest neighbors for each location
"""
function nearest_neighbor_distances(dom::Domain, n_neighbors::Int64)::Vector{Float64}
    dist_matrix = distance_matrix(dom)
    return nearest_neighbor_distances(dist_matrix, n_neighbors)
end

"""
    param_table(d::ADRIADomain)::DataFrame

Get model fieldnames and their parameter values.
"""
function param_table(d::Domain)::DataFrame
    f_names::Vector{String} = collect(string.(d.model[:fieldname]))
    vals::Vector{<:Real} = collect(d.model[:val])
    p_df::DataFrame = DataFrame(OrderedDict(k => v for (k, v) in zip(f_names, vals)))

    return p_df
end

"""
    model_spec(d::Domain)::DataFrame
    model_spec(d::Domain, filepath::String)::Nothing
    model_spec(m::Model)::DataFrame
    model_spec(d::Domain, param_set::DataFrame)
    model_spec(d::Domain, param_names::Vector{String})
    model_spec(d::Domain, param_names::Vector{Symbol})

Get model specification as DataFrame with lower and upper bounds.
If a filepath is provided, writes the specification out to file with ADRIA metadata.
If a `param_set` or a vector of `param_set` names is provided, filters params not present in
the `param_set`.
"""
function model_spec(d::Domain)::DataFrame
    return model_spec(d.model)
end
function model_spec(d::Domain, filepath::String)::Nothing
    version = PkgVersion.Version(@__MODULE__)
    vers_id = "v$(version)"

    current_time = replace(string(now()), "T" => "_", ":" => "_", "." => "_")
    open(filepath, "w") do io
        write(io, "# Generated with ADRIA.jl $(vers_id) on $(current_time)\n")
    end

    ms = model_spec(d)
    ms[!, :] .= replace(Matrix(ms), nothing => "")
    CSV.write(filepath, ms; header=true, append=true)

    return nothing
end
function model_spec(m::Model)::DataFrame
    spec = DataFrame(m)
    # Model parameter distributions have at least one argument
    spec[!, :dist_params] = Vector{Tuple{Float64,Vararg{Float64}}}(spec.dist_params)

    DataFrames.hcat!(
        spec,
        DataFrame(
            :lower_bound => factor_lower_bounds.(eachrow(spec)),
            :upper_bound => factor_upper_bounds.(eachrow(spec))
        )
    )

    spec[!, :component] .= replace.(string.(spec[!, :component]), "ADRIA." => "")
    spec[!, :is_constant] .= spec[!, :lower_bound] .== spec[!, :upper_bound]

    # Reorder so name/description appears at end
    # makes viewing as CSV a little nicer given description can be very long
    select!(spec, Not([:name, :description]), [:name, :description])

    return spec
end
function model_spec(d::Domain, param_set::DataFrame)
    return model_spec(d, Symbol.(names(param_set)))
end
function model_spec(d::Domain, param_names::Vector{String})
    return model_spec(d, Symbol.(param_names))
end
function model_spec(d::Domain, param_names::Vector{Symbol})
    ms = model_spec(d)
    return ms[ms.fieldname .∈ Ref(param_names), :]
end

"""
    update_params!(d::Domain, params::Union{AbstractVector,DataFrameRow})::Nothing

Update given domain with new parameter values.
"""
function update_params!(d::Domain, params::Union{AbstractVector,DataFrameRow})::Nothing
    p_df::DataFrame = model_spec(d, names(params))[
        :, [:fieldname, :val, :ptype, :dist_params]
    ]

    try
        p_df[!, :val] .= collect(params[Not("RCP")])
    catch err
        if isa(err, ArgumentError) || isa(err, DimensionMismatch)
            if !occursin("RCP", "$err")
                error("Error occurred loading scenario samples. $err")
            end
            p_df[!, :val] .= collect(params)
        end
    end

    # Unused params need to be merged with `p_df` to keep the same number of factors as d.model
    ms_all = model_spec(d)
    p_df_complementar::DataFrame = ms_all[
        (ms_all.fieldname .∉ Ref(Symbol.(names(params)))),
        [:fieldname, :val, :ptype, :dist_params]
    ]

    # Update with new parameters
    update!(d.model, vcat(p_df, p_df_complementar))

    return nothing
end

"""
    component_params(dom::Domain, component)::DataFrame
    component_params(m::Model, component)::DataFrame
    component_params(spec::DataFrame, component)::DataFrame
    component_params(m::Model, components::Vector{T})::DataFrame
    component_params(spec::DataFrame, components::Vector{T})::DataFrame

Extract parameters for a specific model component.
"""
function component_params(dom::Domain, component)::DataFrame
    return component_params(model_spec(dom), component)
end
function component_params(m::Model, component)::DataFrame
    return component_params(model_spec(m), component)
end
function component_params(spec::DataFrame, component)::DataFrame
    return spec[spec.component .== replace.(string(component), "ADRIA." => ""), :]
end
function component_params(m::Model, components::Vector{T})::DataFrame where {T}
    return component_params(model_spec(m), components)
end
function component_params(spec::DataFrame, components::Vector{T})::DataFrame where {T}
    return spec[spec.component .∈ [replace.(string.(components), "ADRIA." => "")], :]
end

"""
    _convert_abs_to_k(coral_cover::Union{YAXArray,Matrix{Float64}}, spatial::DataFrame)::Union{YAXArray,Matrix{Float64}}

Convert coral cover data from being relative to absolute location area to relative to
\$k\$ area.
"""
function _convert_abs_to_k(
    coral_cover::Union{YAXArray,Matrix{Float64}}, spatial::DataFrame
)::Union{YAXArray,Matrix{Float64}}
    # Initial coral cover is provided as values relative to location area.
    # Convert coral covers to be relative to k area, ignoring locations with 0 carrying
    # capacity (k area = 0.0).
    absolute_k_area = (spatial.k .* spatial.area)'  # max possible coral area in m^2
    valid_locs::BitVector = absolute_k_area' .> 0.0
    coral_cover[:, valid_locs] .= (
        (coral_cover[:, valid_locs] .* spatial.area[valid_locs]') ./
        absolute_k_area[valid_locs]'
    )

    # Ensure initial coral cover values are <= maximum carrying capacity
    @assert all(sum(coral_cover; dims=1) .<= 1.0)

    return coral_cover
end

"""
    loc_area(domain::Domain)::Vector{Float64}

Get location area for the given domain.
"""
function loc_area(domain::Domain)::Vector{Float64}
    return domain.loc_data.area
end

@deprecate site_k_area(dom) loc_k_area(dom)

"""
    loc_k_area(domain::Domain)::Vector{Float64}

Get maximum coral cover area for the given domain in absolute area.
"""
function loc_k_area(domain::Domain)::Vector{Float64}
    return location_k(domain) .* loc_area(domain)
end

"""
    n_locations(domain::Domain)::Int64

Returns the number of locations (sites/reefs/clusters) represented within the domain.
"""
function n_locations(domain::Domain)::Int64
    return size(domain.loc_data, 1)
end

"""
    relative_leftover_space(loc_cover::AbstractArray)::AbstractArray

Get proportion of leftover space, given site_k and proportional cover on each site, summed
over species.

# Arguments
- `loc_cover` : Proportion of coral cover relative to `k` (maximum carrying capacity).

# Returns
Leftover space ∈ [0, 1]
"""
function relative_leftover_space(
    loc_cover::AbstractArray
)::AbstractArray
    return max.(1.0 .- loc_cover, 0.0)
end

"""
    loc_coral_cover(C_cover_t::Array{Float64,3})::Vector{Float64}

Sum coral cover across all functional groups and size classes of a single timestep for each location.
"""
function loc_coral_cover(C_cover_t::AbstractArray{Float64,3})::Vector{Float64}
    return dropdims(sum(C_cover_t; dims=(1, 2)); dims=(1, 2))
end

"""
    loc_recruits_cover(recruits::Matrix{Float64})::Vector{Float64}

Absolute cover of recruits on each location.
"""
function loc_recruits_cover(recruits::Matrix{Float64})::Vector{Float64}
    return vec(sum(recruits; dims=1))
end

"""
    location_k(domain::Domain)::Vector{Float64}

Get maximum coral habitable area as a proportion of a location's area (\$k ∈ [0, 1]\$).
"""
function location_k(domain::Domain)
    return domain.loc_data.k
end

"""Extract the time steps represented in the data package."""
function timesteps(domain::Domain)
    return domain.env_layer_md.timeframe
end

"""
    update!(dom::Domain, spec::DataFrame)::Nothing

Update a Domain model with new values specified in spec.
Assumes all `val` and `bounds` are to be updated.

# Arguments
- `dom` : Domain
- `spec` : updated model specification
"""
function update!(dom::Domain, spec::DataFrame)::Nothing
    dom.model[:val] = spec.val
    dom.model[:dist_params] = spec.dist_params

    return nothing
end

# Dummy interface to allow precompilation
function switch_RCPs!() end

"""
    set_seed_target_locations!(domain::Domain, location_ids::Vector{String})

Set the locations eligible for seeding interventions.

# Arguments
- `domain`: Domain to modify
- `location_ids`: Vector of location IDs to target for seeding

# Example
```julia
dom = ADRIA.load_domain("path/to/domain")
# Only seed in marine park zones
ADRIA.set_seed_target_locations!(dom, ["reef_01", "reef_05", "reef_12"])
```
"""
function set_seed_target_locations!(domain::Domain, location_ids::Vector{String})
    # Validate that all locations exist in domain
    invalid_locs = setdiff(location_ids, domain.loc_ids)
    if !isempty(invalid_locs)
        error("Invalid location IDs: $(join(invalid_locs, ", "))")
    end

    domain.seed_target_locations = location_ids
    return nothing
end

"""
    set_fog_target_locations!(domain::Domain, location_ids::Vector{String})

Set the locations eligible for fogging interventions.

# Arguments
- `domain`: Domain to modify
- `location_ids`: Vector of location IDs to target for fogging

# Example
```julia
dom = ADRIA.load_domain("path/to/domain")
# Only fog high-value tourism reefs
ADRIA.set_fog_target_locations!(dom, ["reef_03", "reef_07"])
```
"""
function set_fog_target_locations!(domain::Domain, location_ids::Vector{String})
    # Validate that all locations exist in domain
    invalid_locs = setdiff(location_ids, domain.loc_ids)
    if !isempty(invalid_locs)
        error("Invalid location IDs: $(join(invalid_locs, ", "))")
    end

    domain.fog_target_locations = location_ids
    return nothing
end
