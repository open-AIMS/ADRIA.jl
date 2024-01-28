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

function load_mat_data(
    data_fn::String, attr::String, site_data::DataFrame
)::NamedDimsArray{Float32}
    data = matread(data_fn)
    local loaded::NamedDimsArray{Float32}
    local site_order::Vector{String}

    # Attach site names to each dimension
    try
        site_order = Vector{String}(vec(data["reef_siteid"]))
        loaded = NamedDimsArray(data[attr]; Source=site_order, Sink=site_order)
    catch err
        if isa(err, KeyError)
            @warn "Provided file $(data_fn) did not have reef_siteid! There may be a mismatch in sites."
            if size(loaded, 2) != nrow(site_data)
                @warn "Mismatch in number of sites ($(data_fn)).\nTruncating so that data size matches!"

                # Subset down to number of sites
                tmp = selectdim(data[attr], 2, 1:nrow(site_data))
                loaded = NamedDimsArray(tmp; Source=site_order, Sink=site_order)
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
    local loaded_data::NamedDimsArray

    NetCDF.open(data_fn; mode=NC_NOWRITE) do nc_file
        data::Array{<:AbstractFloat} = NetCDF.readvar(nc_file, attr)
        dim_names::Vector{Symbol} = Symbol[
            Symbol(dim.name) for dim in nc_file.vars[attr].dim
        ]
        dim_labels::Vector{Union{UnitRange{Int64},Vector{String}}} = _nc_dim_labels(
            data_fn, data, nc_file
        )

        try
            loaded_data = NamedDimsArray(data; zip(dim_names, dim_labels)...)
        catch err
            if isa(err, KeyError)
                n_sites = size(data, 2)
                @warn "Provided file $(data_fn) did not have the expected dimensions " *
                    "(one of: timesteps, reef_siteid, scenarios)."
                if n_sites != nrow(site_data)
                    error(
                        "Mismatch in number of sites ($(data_fn)). " *
                        "Expected $(nrow(site_data)), got $(n_sites)",
                    )
                end
            else
                rethrow(err)
            end
        end
    end

    return loaded_data
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
        # This will be an issue if two or more dimensions have the same number of elements
        # as the number of sites, but so far it hasn't happened...
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
    load_covers(data_fn::String, attr::String, site_data::DataFrame)::NamedDimsArray

Load initial coral cover data from netCDF.
"""
function load_covers(data_fn::String, attr::String, site_data::DataFrame)::NamedDimsArray
    data::NamedDimsArray = load_nc_data(data_fn, attr, site_data)
    data = NamedDims.rename(data, :covers => :species, :reef_siteid => :sites)

    # Reorder sites to match site_data
    data = data[sites=Key(site_data[:, :reef_siteid])]
    data = _convert_abs_to_k(data, site_data)

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

"""
    load_cyclone_mortality(data_fn::String)::NamedDimsArray
    load_cyclone_mortality(timeframe::Vector{Int64}, site_data::DataFrame)::NamedDimsArray

Load cyclone mortality datacube from NetCDF file. The returned cyclone_mortality datacube is
ordered by :locations
"""
function load_cyclone_mortality(data_fn::String)::NamedDimsArray
    # Read file as YAXArray and convert to NamedDimsArray
    # as we intend to move to use only YAXArrays soon
    cyclone_cube::YAXArray = Cube(data_fn)

    # Get locations sort indexes
    locations = collect(cyclone_cube.locations)
    location_sort_idx::Vector{Int64} = sortperm(locations)

    # Order locations labels
    ordered_locations::Vector{String} = locations[location_sort_idx]
    YAXArrays.DD.set(cyclone_cube, :locations => ordered_locations)

    # Order data according to locations
    ordered_data = cyclone_cube.data[:, location_sort_idx, :, :]

    ordered_cyclone_cube::YAXArray = YAXArray(cyclone_cube.axes, ordered_data)
    return _yaxarray2nameddimsarray(ordered_cyclone_cube)
end
function load_cyclone_mortality(
    timeframe::Vector{Int64}, site_data::DataFrame
)::NamedDimsArray
    locations::Vector{String} = sort(site_data.reef_siteid)
    species::Vector{String} = ADRIA.coral_spec().taxa_names
    scenarios::Vector{Int64} = [1]

    n_timesteps::Int64 = length(timeframe)
    n_locations::Int64 = length(locations)
    n_species::Int64 = length(species)
    n_scenarios::Int64 = length(scenarios)

    axlist::Tuple = (
        Dim{:timesteps}(1:n_timesteps),
        Dim{:locations}(locations),
        Dim{:species}(species),
        Dim{:scenarios}(scenarios),
    )

    data::Array{Float64,4} = zeros(n_timesteps, n_locations, n_species, n_scenarios)

    return _yaxarray2nameddimsarray(YAXArray(axlist, data))
end

function _yaxarray2nameddimsarray(yarray::YAXArray)::NamedDimsArray
    data = yarray.data

    dim_names::NTuple{4,Symbol} = name.(yarray.axes)
    dim_labels::Vector = lookup.([yarray], dim_names)

    return NamedDimsArray(data; zip(dim_names, dim_labels)...)
end
