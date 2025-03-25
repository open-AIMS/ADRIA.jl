using Distributions
using JSON
using NetCDF
using YAXArrays

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
            error("""Incompatible Domain data package. Detected $(dpkg_version),
            but only support one of $(COMPAT_DPKG)""")
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
    df = CSV.read(filepath, DataFrame; comment="#")

    if columnindex(df, :RCP) > 0
        df = df[!, Not("RCP")]
    end

    return df
end

"""
    load_nc_data(data_fn::String, attr::String; dim_names::Vector{Symbol}=Symbol[], dim_names_replace::Vector{Pair{Symbol,Symbol}}=Pair{Symbol,Symbol}[])::YAXArray

Load cluster-level data for a given attribute in a netCDF.
"""
function load_nc_data(
    data_fn::String,
    attr::String;
    dim_names::Vector{Symbol}=Symbol[],
    dim_names_replace::Vector{Pair{Symbol,Symbol}}=Pair{Symbol,Symbol}[]
)::YAXArray
    data = try
        sort_axis(Cube(data_fn), :locations)
    catch
        fallback_nc_data(data_fn, attr; dim_names, dim_names_replace)
    end

    return data
end
function fallback_nc_data(
    data_fn::String,
    attr::String;
    dim_names::Vector{Symbol}=Symbol[],
    dim_names_replace::Vector{Pair{Symbol,Symbol}}=Pair{Symbol,Symbol}[]
)::YAXArray
    NetCDF.open(data_fn; mode=NC_NOWRITE) do nc_file
        data::Array{<:AbstractFloat} = NetCDF.readvar(nc_file, attr)

        if isempty(dim_names)
            dim_names = [Symbol(dim.name) for dim in nc_file.vars[attr].dim]
        end

        if !isempty(dim_names_replace)
            replace!(dim_names, dim_names_replace...)
        end

        dim_labels = _nc_dim_labels(data_fn, data, nc_file)
        return sort_axis(DataCube(data; zip(dim_names, dim_labels)...), :sites)
    end
end

"""
Return vector of labels for each dimension.

**Important:** Cannot trust indicated dimension metadata to get site labels,
because this can be incorrect. Instead, match by number of sites.
"""
function _nc_dim_labels(
    data_fn::String, data::Array{<:Real}, nc_file::NetCDF.NcFile
)::Vector{Union{UnitRange{Int64},Vector{String}}}
    local locs_idx::Int64

    sites = "reef_siteid" in keys(nc_file.vars) ? _site_labels(nc_file) : 1:size(data, 2)

    try
        # This will be an issue if the number of elements for two or more dimensions have
        # the same number of elements, but so far that hasn't happened...
        locs_idx = first(findall(size(data) .== length(sites)))
    catch err
        error(
            "Error loading $data_fn : could not determine number of locations." *
            "Detected size: $(size(data)) | Known number of locations: $(length(sites))"
        )
    end

    dim_labels = Union{UnitRange{Int64},Vector{String}}[1:n for n in size(data)]
    dim_labels[locs_idx] = sites

    time_idx = findfirst([k == "timesteps" for k in keys(nc_file.dim)])
    if !isnothing(time_idx)
        time_vals =
            "timesteps" in keys(
                nc_file.vars
            ) ? NetCDF.readvar(nc_file["timesteps"]) :
            1:Int64(nc_file.dim["timesteps"].dimlen)
        dim_labels[1] = minimum(time_vals):maximum(time_vals)
    end

    return dim_labels
end

"""
Some packages used to write out netCDFs do not yet support string values, and instead
reverts to writing out character arrays.
"""
function _site_labels(nc_file::NetCDF.NcFile)::Vector{String}
    loc_ids = NetCDF.readvar(nc_file, "reef_siteid")
    # Converts character array entries in netCDFs to string if needed
    return loc_ids isa Matrix ? nc_char2string(loc_ids) : loc_ids
end

"""
    load_env_data(data_fn::String, attr::String)::YAXArray
    load_env_data(timeframe, sites, n_scenarios)::YAXArray

Load environmental data layers (DHW, Wave) from netCDF.
"""
function load_env_data(data_fn::String, attr::String, timeframe::Vector{Int64})::YAXArray
    _dim_names::Vector{Symbol} = [:timesteps, :sites, :scenarios]
    env_data::YAXArray = load_nc_data(data_fn, attr; dim_names=_dim_names)
    env_data_tf = collect(env_data.timesteps)

    if all(timeframe .∈ Ref(env_data_tf))
        return env_data[timesteps=At(timeframe)]
    elseif all(env_data_tf .== 1:length(env_data_tf))
        # if the data file has timesteps starting at 1 return the environmental data
        # starting at one to the length of the input timeframe
        return env_data[timesteps=1:length(timeframe)]
    end

    msg = "Timeframe in environmental netcdf does not contain given timeframe."
    throw(ArgumentError(msg))
end
function load_env_data(
    timeframe::Vector{Int64}, sites::Vector{String}, n_scenarios::Int64
)::YAXArray
    return ZeroDataCube(;
        T=Float32, timesteps=timeframe, sites=sites, scenarios=1:n_scenarios
    )
end

"""
    load_cyclone_mortality(data_fn::String, timeframe)::YAXArray
    load_cyclone_mortality(timeframe::Vector{Int64}, loc_data::DataFrame)::YAXArray

Load cyclone mortality datacube from NetCDF file. The returned cyclone_mortality datacube is
ordered by :locations
"""
function load_cyclone_mortality(data_fn::String, timeframe::Vector{Int64})::YAXArray
    cyclone_cube::YAXArray = Cube(data_fn)
    cyclone_cube = sort_axis(cyclone_cube, :locations)
    cyclone_tf = collect(cyclone_cube.timesteps)

    # Return cyclones with the correct timeframe specified in the data package.
    if all(timeframe .∈ Ref(cyclone_tf))
        return cyclone_cube[timesteps=At(timeframe)]
    elseif all(cyclone_tf .== 1:length(cyclone_tf))
        # if the data file has timesteps starting at 1 return the cyclone data
        # starting at one to the length of the input timeframe
        return cyclone_cube[timesteps=1:length(timeframe)]
    end

    throw(ArgumentError("Timeframe in cyclone netcdf does not contain given timeframe."))
end
function load_cyclone_mortality(
    timeframe::Vector{Int64}, loc_data::DataFrame, location_id_col::String
)::YAXArray
    return ZeroDataCube(;
        timesteps=1:length(timeframe),
        locations=loc_data[:, location_id_col],
        species=ADRIA.coral_spec().taxa_names,
        scenarios=[1]
    )
end
