using JSON
using NetCDF

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
    dim_names_replace::Vector{Pair{Symbol,Symbol}}=Pair{Symbol,Symbol}[],
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
    local sites_idx::Int64

    sites = "reef_siteid" in keys(nc_file.vars) ? _site_labels(nc_file) : 1:size(data, 2)

    try
        # This will be an issue if number of elements is equal the number of sites, two or
        # more dimensions but so far that hasn't happened...
        sites_idx = first(findall(size(data) .== length(sites)))
    catch err
        error(
            "Error loading $data_fn : could not determine number of locations." *
            "Detected size: $(size(data)) | Known number of locations: $(length(sites))",
        )
    end

    dim_labels = Union{UnitRange{Int64},Vector{String}}[1:n for n in size(data)]
    dim_labels[sites_idx] = sites

    return dim_labels
end

"""
Some packages used to write out netCDFs do not yet support string values, and instead
reverts to writing out character arrays.
"""
function _site_labels(nc_file::NetCDF.NcFile)::Vector{String}
    site_ids = NetCDF.readvar(nc_file, "reef_siteid")
    # Converts character array entries in netCDFs to string if needed
    return site_ids isa Matrix ? nc_char2string(site_ids) : site_ids
end

"""
    load_cover(data_fn::String, site_data::DataFrame)::YAXArray
    load_cover(n_species::Int64, n_sites::Int64)::YAXArray

Load initial coral cover data from netCDF.
"""
function load_cover(data_fn::String, site_data::DataFrame)::YAXArray
    _dim_names_replace = [:covers => :species, :reef_siteid => :sites]
    data = load_nc_data(data_fn, "covers"; dim_names_replace=_dim_names_replace)

    return _convert_abs_to_k(data, site_data)
end
function load_cover(n_species::Int64, n_sites::Int64)::YAXArray
    @warn "Using random initial coral cover"

    return DataCube(rand(Float32, n_species, n_sites); species=1:(n_species), sites=1:n_sites)
end

"""
    load_env_data(data_fn::String, attr::String, site_data::DataFrame)::YAXArray
    load_env_data(timeframe, sites)::YAXArray

Load environmental data layers (DHW, Wave) from netCDF.
"""
function load_env_data(data_fn::String, attr::String)::YAXArray
    _dim_names::Vector{Symbol} = [:timesteps, :sites, :scenarios]
    return load_nc_data(data_fn, attr; dim_names=_dim_names)
end
function load_env_data(timeframe::Vector{Int64}, sites::Vector{String})::YAXArray
    return ZeroDataCube(Float32; timesteps=timeframe, sites=sites, scenarios=1:50)
end

"""
    load_cyclone_mortality(data_fn::String)::YAXArray
    load_cyclone_mortality(timeframe::Vector{Int64}, site_data::DataFrame)::YAXArray

Load cyclone mortality datacube from NetCDF file. The returned cyclone_mortality datacube is
ordered by :locations
"""
function load_cyclone_mortality(data_fn::String)::YAXArray
    cyclone_cube::YAXArray = Cube(data_fn)
    return sort_axis(cyclone_cube, :locations)
end
function load_cyclone_mortality(timeframe::Vector{Int64}, site_data::DataFrame)::YAXArray
    return ZeroDataCube(;
        timesteps=1:length(timeframe),
        locations=sort(site_data.reef_siteid),
        species=ADRIA.coral_spec().taxa_names,
        scenarios=[1]
    )
end

"""
    DataCube(data::AbstractArray; kwargs...)::YAXArray

Constructor for YAXArray.
"""
function DataCube(data::AbstractArray; kwargs...)::YAXArray
    return YAXArray(Tuple(Dim{name}(val) for (name, val) in kwargs), data)
end

"""
    ZeroDataCube(T=Float64; kwargs...)::YAXArray

Constructor for YAXArray with all entries equal zero.
"""
function ZeroDataCube(T=Float64; kwargs...)::YAXArray
    return DataCube(zeros(T, [length(val) for (name, val) in kwargs]...); kwargs...)
end

function axes_names(cube::YAXArray)
    return name.(cube.axes)
end

"""
    sort_axis(cube::YAXArray, axis_name::Symbol)

Sorts axis labels of a given YAXArray datacube.
"""
function sort_axis(cube::YAXArray, axis_name::Symbol)::YAXArray
    axis_labels = collect(lookup(cube, axis_name))
    labels_sort_idx::Vector{Int64} = sortperm(axis_labels)
    ordered_labels = axis_labels[labels_sort_idx]

    # Get selector to order data
    _axes_names = axes_names(cube)
    n_axes = length(_axes_names)
    selector::Vector{Union{Colon,Vector{Int64}}} = fill(:, n_axes)
    axis_idx::Int64 = findfirst(x -> x == axis_name, _axes_names)
    selector[axis_idx] = labels_sort_idx

    return cube[selector...]
end

"""
    copy(cube::YAXArray)::YAXArray

Copy a YAXArray data cube.
"""
function copy_datacube(cube::YAXArray)::YAXArray
    new_axlist = Tuple(ax for ax in deepcopy(cube.axes))
    return YAXArray(new_axlist, copy(cube.data))
end

function yaxarray2nameddimsarray(cube::YAXArray)::NamedDimsArray
    axis_names = axes_names(cube)
    axis_labels = collect.(lookup(cube, axis_names))
    return NamedDimsArray(cube.data; NamedTuple{axis_names}(axis_labels)...)
end
