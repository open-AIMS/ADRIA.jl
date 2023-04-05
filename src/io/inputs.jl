using JSON


"""
    _check_compat(dpkg_details::Dict)

Checks for version compatibility.

# Arguments
- `dpkg_details` : Datapackage spec
"""
function _check_compat(dpkg_details::Dict{String,Any})::Nothing
    if haskey(dpkg_details, "version") || haskey(dpkg_details, "dpkg_version")
        dpkg_version::String = dpkg_details["version"]
        if dpkg_version ∉ COMPAT_DPKG
            error("Incompatible Domain data package. Detected $(dpkg_version), but only support one of $(COMPAT_DPKG)")
        end
    else
        error("Incompatible Domain data package.")
    end

    return nothing
end

"""
    _load_dpkg(dpkg_path::String)

Load and parse datapackage.

# Arguments
- `dpkg_path` : path to datapackage
"""
function _load_dpkg(dpkg_path::String)::Dict{String,Any}
    local dpkg_md::Dict{String,Any}
    open(joinpath(dpkg_path, "datapackage.json"), "r") do fp
        dpkg_md = JSON.parse(read(fp, String))
    end
    _check_compat(dpkg_md)

    return dpkg_md
end

"""
    load_scenarios(domain::Domain, filepath::String)::DataFrame

Load and pre-process scenario values.
Parameters intended to be of Integer type or casted as such.
"""
function load_scenarios(domain::Domain, filepath::String)::DataFrame
    df = CSV.read(filepath, DataFrame, comment="#")

    if columnindex(df, :RCP) > 0
        df = df[!, Not("RCP")]
    end
    _process_inputs!(domain, df)

    return df
end

"""
    _process_inputs!(d::Domain, df::DataFrame)::Nothing
    _process_inputs!(spec::DataFrame, df::DataFrame)::Nothing
    _process_inputs!(bnds::AbstractArray, p_types::AbstractArray, df::DataFrame)::Nothing
    
Map sampled values in `df` back to discrete bounds for parameters
indicated to be of integer type in the Domain spec.

# Arguments
- `d` : Domain
- `bnds` : Tuple containing list of parameter bounds (lower,upper).
- `spec` : Model specification defining parameter bounds, types and distributions.
- `df` : parameter sets defining scenarios
"""
function _process_inputs!(d::Domain, df::DataFrame)::Nothing
    return _process_inputs!(d.model[:bounds], d.model[:ptype], df)
end
function _process_inputs!(spec::DataFrame, df::DataFrame)::Nothing
    return _process_inputs!(Tuple(spec[:, :bounds]), Tuple(spec[:, :ptype]), df)
end
function _process_inputs!(bnds::Tuple, p_types::Tuple, df::DataFrame)::Nothing
    for (i, dt) in enumerate(p_types)
        if dt == "integer" && (bnds[i][1] < bnds[i][2])
            @inbounds df[!, i] .= map_to_discrete.(df[!, i], Int64(bnds[i][2]))
        end
    end
    return nothing
end

function load_mat_data(data_fn::String, attr::String, site_data::DataFrame)::NamedDimsArray{Float32}
    data = matread(data_fn)
    local loaded::NamedDimsArray{Float32}
    local site_order::Vector{String}

    # Attach site names to each dimension
    try
        site_order = Vector{String}(vec(data["reef_siteid"]))
        loaded = NamedDimsArray(data[attr]; Source=site_order, Receiving=site_order)
    catch err
        if isa(err, KeyError)
            @warn "Provided file $(data_fn) did not have reef_siteid! There may be a mismatch in sites."
            if size(loaded, 2) != nrow(site_data)
                @warn "Mismatch in number of sites ($(data_fn)).\nTruncating so that data size matches!"

                # Subset down to number of sites
                tmp = selectdim(data[attr], 2, 1:nrow(site_data))
                loaded = NamedDimsArray(tmp; Source=site_order, Receiving=site_order)
            end
        else
            rethrow(err)
        end
    end

    return loaded
end


"""
    load_nc_data(data_fn::String, attr::String, site_data::DataFrame)::NamedDimsArray

Load cluster-level data for a given attribute in a netCDF as a NamedDimsArray.
"""
function load_nc_data(data_fn::String, attr::String, site_data::DataFrame)::NamedDimsArray
    local loaded::NamedDimsArray

    ds = Dataset(data_fn, "r")
    dim_names = keys(dimsize(ds[attr]))
    data = ds[attr][:, :]

    if "reef_siteid" in keys(ds)
        sites = _char_to_string(ds["reef_siteid"][:])
    else
        sites = 1:size(data, 2)
    end
    close(ds)

    # Note: cannot trust indicated dimension metadata
    # because this can be incorrect!
    # instead, match by number of sites
    dim_keys = Union{UnitRange{Int64},Vector{String}}[1:n for n in size(data)]
    i = first(findall(size(data) .== length(sites)))
    dim_keys[i] = sites

    try
        loaded = NamedDimsArray(data; zip(dim_names, dim_keys)...)
    catch err
        if isa(err, KeyError)
            n_sites = size(data, 2)
            @warn "Provided file $(data_fn) did not have the expected dimensions (one of: timesteps, reef_siteid, scenarios)."
            if n_sites != nrow(site_data)
                error("Mismatch in number of sites ($(data_fn)). Expected $(nrow(site_data)), got $(n_sites)")
            end
        else
            rethrow(err)
        end
    end

    return loaded
end

"""
    _char_to_string(vals)::Vector{String}

Convert character array entries in netCDFs to string.

Some packages used to write out netCDFs do not yet support
string values, and instead reverts to writing out character arrays.
"""
function _char_to_string(vals::AbstractVecOrMat)::Vector{String}
    if vals isa Matrix
        vals = map(x -> join(skipmissing(x)), eachcol(vals))
    end

    # R's ncdf4 package does not yet support string values
    # so strip the null terminator from the joined character array.
    vals = replace.(vals, "\0" => "")

    return vals
end


"""
    load_covers(data_fn::String, attr::String, site_data::DataFrame)::NamedDimsArray

Load initial coral cover data from netCDF.
"""
function load_covers(data_fn::String, attr::String, site_data::DataFrame)::NamedDimsArray
    data::NamedDimsArray = load_nc_data(data_fn, attr, site_data)
    data = NamedDims.rename(data, :covers => :species, :reef_siteid => :sites)

    # Reorder sites to match site_data
    data = data[sites=Key(site_data[:, :reef_siteid])]

    return data
end


"""
    load_env_data(data_fn::String, attr::String, site_data::DataFrame)::NamedDimsArray

Load environmental data layers (DHW, Wave) from netCDF.
"""
function load_env_data(data_fn::String, attr::String, site_data::DataFrame)::NamedDimsArray
    data::NamedDimsArray = load_nc_data(data_fn, attr, site_data)

    # Re-attach correct dimension names
    data = NamedDims.rename(data, (:timesteps, :sites, :scenarios))

    # Reorder sites to match site_data
    data = data[sites=Key(site_data[:, :reef_siteid])]

    return data
end
