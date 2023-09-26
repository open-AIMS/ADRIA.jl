using NCDatasets


abstract type Domain end


"""
    EnvLayer{S, TF}

Store environmental data layers used for scenario
"""
mutable struct EnvLayer{S<:AbstractString,TF}
    dpkg_path::S
    site_data_fn::S
    const site_id_col::S
    const unique_site_id_col::S
    init_coral_cov_fn::S
    connectivity_fn::S
    DHW_fn::S
    wave_fn::S
    const timeframe::TF
end


"""
    site_distance(site_data::DataFrame)::Matrix

Calculate matrix of unique distances between sites.

# Returns
tuple, matrix of distance between sites, median site distance for domain
"""

function site_distances(site_data::DataFrame)::Tuple{Matrix{Float64},Float64}
    site_centroids = centroids(site_data)
    longitudes = first.(site_centroids)
    latitudes = last.(site_centroids)

    n_sites = size(site_data, 1)
    dist = fill(NaN, n_sites, n_sites)
    for jj in axes(dist, 2)
        for ii in axes(dist, 1)
            if ii == jj
                continue
            end

            @inbounds dist[ii, jj] = haversine((longitudes[ii], latitudes[ii]), (longitudes[jj], latitudes[jj]))
        end
    end

    median_site_dist = median(dist[.!isnan.(dist)])
    return dist, median_site_dist
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
    return d.site_data[:, d.unique_site_id_col]
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

    open(filepath, "w") do io
        write(io, "# Generated with ADRIA.jl $(vers_id) on $(replace(string(now()), "T"=>"_", ":"=>"_", "."=>"_"))\n")
    end

    CSV.write(filepath, model_spec(d); header=true, append=true)

    return nothing
end
function model_spec(m::Model)::DataFrame
    spec = DataFrame(m)
    bnds = spec[!, :bounds]

    DataFrames.hcat!(spec, DataFrame(
        :lower_bound => first.(bnds),
        :upper_bound => getindex.(bnds, 2)
    ))

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
Maps sampled continuous values to discrete values for categorical variables.
"""
function update_params!(d::Domain, params::Union{AbstractVector,DataFrameRow})::Nothing
    p_df::DataFrame = DataFrame(d.model)[:, [:fieldname, :val, :ptype, :bounds]]

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

    to_floor = (p_df.ptype .== "integer") .| (p_df.ptype .== "categorical")
    if any(to_floor)
        p_df[to_floor, :val] .= map_to_discrete.(p_df[to_floor, :val], Int64.(getindex.(p_df[to_floor, :bounds], 2)))
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
    return spec[spec.component.==replace.(string(component), "ADRIA." => ""), :]
end
function component_params(m::Model, components::Vector{T})::DataFrame where {T}
    return component_params(model_spec(m), components)
end
function component_params(spec::DataFrame, components::Vector{T})::DataFrame where {T}
    return spec[spec.component.∈[replace.(string.(components), "ADRIA." => "")], :]
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
    return site_k(domain) .* site_area(domain)
end

"""
    n_locations(domain::Domain)::Int64

Returns the number of locations (sites/reefs/clusters) represented within the domain.
"""
function n_locations(domain::Domain)::Int64
    return size(domain.site_data, 1)
end

"""
    relative_leftover_space(domain::Domain)::Vector{Float64}
    relative_leftover_space(site_k::Matrix{Float64}, site_coral_cover::Matrix{Float64})::Matrix{Float64}

Get proportion of leftover space, given site_k and proportional cover on each site, summed over species.
"""
function relative_leftover_space(domain::Domain, site_coral_cover::Matrix{Float64})::Matrix{Float64}
    return relative_leftover_space(site_k(domain)', site_coral_cover)
end
function relative_leftover_space(site_k::AbstractArray{Float64,2}, site_coral_cover::Matrix{Float64})::Matrix{Float64}
    return max.(site_k .- site_coral_cover, 0.0)
end

"""
    site_k(domain::Domain)::Vector{Float64}

Get maximum coral cover area as a proportion of site area.
"""
function site_k(domain::Domain)::Vector{Float64}
    return domain.site_data.k ./ 100.0
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
    dom.model[:bounds] = spec.bounds

    return nothing
end
