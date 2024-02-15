abstract type Domain end

"""
    EnvLayer{S, TF}

Store environmental data layers used for scenario
"""
mutable struct EnvLayer{S<:AbstractString,TF}
    dpkg_path::S
    site_data_fn::S
    const site_id_col::S
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
function load_domain(path::String)::Domain
    return load_domain(path, "")
end

function unique_sites(d::Domain)::Vector{String}
    return d.site_data[:, d.cluster_id_col]
end

"""
    param_table(d::ADRIADomain)::DataFrame

Get model fieldnames and their parameter values.
"""
function param_table(d::Domain)::DataFrame
    f_names::Vector{String} = collect(string.(d.model[:fieldname]))
    vals::Vector{<:Real} = collect(d.model[:val])
    p_df::DataFrame = DataFrame(OrderedDict(k => v for (k, v) in zip(f_names, vals)))

    p_df[!, :RCP] .= d.RCP  # Add entry to indicate which RCP scenario was used

    return p_df
end

"""
    model_spec(d::Domain)::DataFrame
    model_spec(d::Domain, filepath::String)::Nothing
    model_spec(m::Model)::DataFrame

Get model specification as DataFrame with lower and upper bounds.
If a filepath is provided, writes the specification out to file with ADRIA metadata.
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

    CSV.write(filepath, model_spec(d); header=true, append=true)

    return nothing
end
function model_spec(m::Model)::DataFrame
    spec = DataFrame(m)
    dist_params = spec[!, :dist_params]

    DataFrames.hcat!(
        spec, DataFrame(
            :lower_bound => first.(dist_params),
            :upper_bound => getindex.(dist_params, 2)
        )
    )

    spec[!, :component] .= replace.(string.(spec[!, :component]), "ADRIA." => "")
    spec[!, :is_constant] .= spec[!, :lower_bound] .== spec[!, :upper_bound]

    # Reorder so name/description appears at end
    # makes viewing as CSV a little nicer given description can be very long
    select!(spec, Not([:name, :description]), [:name, :description])

    return spec
end

"""
    update_params!(d::Domain, params::Union{AbstractVector,DataFrameRow})::Nothing

Update given domain with new parameter values.
"""
function update_params!(d::Domain, params::Union{AbstractVector,DataFrameRow})::Nothing
    p_df::DataFrame = DataFrame(d.model)[:, [:fieldname, :val, :ptype, :dist_params]]

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

    # Update with new parameters
    update!(d.model, p_df)

    return nothing
end

"""
    component_params(m::Model, component)::DataFrame
    component_params(spec::DataFrame, component)::DataFrame
    component_params(m::Model, components::Vector)::DataFrame
    component_params(spec::DataFrame, components::Vector)::DataFrame

Extract parameters for a specific model component.
"""
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
    _convert_abs_to_k(coral_cover::Union{NamedDimsArray,Matrix{Float64}}, site_data::DataFrame)::Union{NamedDimsArray,Matrix{Float64}}

Convert coral cover data from being relative to absolute location area to relative to
\$k\$ area.
"""
function _convert_abs_to_k(
    coral_cover::Union{NamedDimsArray,Matrix{Float64}}, site_data::DataFrame
)::Union{NamedDimsArray,Matrix{Float64}}
    # Initial coral cover is provided as values relative to location area.
    # Convert coral covers to be relative to k area, ignoring locations with 0 carrying
    # capacity (k area = 0.0).
    absolute_k_area = (site_data.k .* site_data.area)'  # max possible coral area in m^2
    valid_locs::BitVector = absolute_k_area' .> 0.0
    coral_cover[:, valid_locs] .= (
        (coral_cover[:, valid_locs] .* site_data.area[valid_locs]') ./
        absolute_k_area[valid_locs]'
    )

    # Ensure initial coral cover values are <= maximum carrying capacity
    @assert all(sum(coral_cover; dims=1) .<= 1.0)

    return coral_cover
end

"""
    site_area(domain::Domain)::Vector{Float64}

Get site area for the given domain.
"""
function site_area(domain::Domain)::Vector{Float64}
    return domain.site_data.area
end

"""
    site_k_area(domain::Domain)::Vector{Float64}

Get maximum coral cover area for the given domain in absolute area.
"""
function site_k_area(domain::Domain)::Vector{Float64}
    return location_k(domain) .* site_area(domain)
end

"""
    n_locations(domain::Domain)::Int64

Returns the number of locations (sites/reefs/clusters) represented within the domain.
"""
function n_locations(domain::Domain)::Int64
    return size(domain.site_data, 1)
end

"""
    relative_leftover_space(loc_coral_cover::AbstractArray)::AbstractArray

Get proportion of leftover space, given site_k and proportional cover on each site, summed
over species.

# Arguments
- `loc_coral_cover` : Proportion of coral cover relative to `k` (maximum carrying capacity).

# Returns
Leftover space ∈ [0, 1]
"""
function relative_leftover_space(
    loc_coral_cover::AbstractArray,
)::AbstractArray
    return max.(1.0 .- loc_coral_cover, 0.0)
end

"""
    location_k(domain::Domain)::Vector{Float64}

Get maximum coral habitable area as a proportion of a location's area (\$k ∈ [0, 1]\$).
"""
function location_k(domain::Domain)
    return domain.site_data.k
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
