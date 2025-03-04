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
    load_env_data(timeframe, sites)::YAXArray

Load environmental data layers (DHW, Wave) from netCDF.
"""
function load_env_data(data_fn::String, attr::String)::YAXArray
    _dim_names::Vector{Symbol} = [:timesteps, :sites, :scenarios]
    return load_nc_data(data_fn, attr; dim_names=_dim_names)
end
function load_env_data(timeframe::Vector{Int64}, sites::Vector{String})::YAXArray
    return ZeroDataCube(; T=Float32, timesteps=timeframe, sites=sites, scenarios=1:50)
end

"""
    load_cyclone_mortality(data_fn::String)::YAXArray
    load_cyclone_mortality(timeframe::Vector{Int64}, loc_data::DataFrame)::YAXArray

Load cyclone mortality datacube from NetCDF file. The returned cyclone_mortality datacube is
ordered by :locations
"""
function load_cyclone_mortality(data_fn::String)::YAXArray
    cyclone_cube::YAXArray = Cube(data_fn)
    return sort_axis(cyclone_cube, :locations)
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

function load_cyclone_mortality_from_category(
    data_fn::String, 
    spatial_data::DataFrame, 
    timeframe::Vector{Int64}
)::YAXArray
    cyclone_categories::YAXArray = Cube(data_fn)[
        locations=At(spatial_data.UNIQUE_ID), 
        timesteps=At(timeframe)
    ]

    return _cyclone_mortality_from_category(
        cyclone_categories, collect(cyclone_categories.locations), spatial_data
    )
end

"""Convert cyclone categories to mortalities."""
function _cyclone_mortality_from_category(
    cyclone::YAXArray,
    location_ids::Vector{String},
    spatial_data::DataFrame
)::YAXArray
    # Convert categories into indices for cyclone_mr
    cyclone_cats = cyclone .+ 1
    species = functional_group_names()
    cyclone_mortality::YAXArray{Float64} = ZeroDataCube(;
        T=Float64,
        timesteps=collect(cyclone_cats.timesteps),
        locations=location_ids,
        species=species,
        scenarios=1:length(cyclone_cats.scenario)
    )

    # Mortality rate was generated using rrap_dg
    cyclone_mr::NamedTuple = _cyclone_mortalities()

    # To filter massives/branchings
    massives::BitVector = contains.(String.(species), ["massives"])
    branchings::BitVector = .!massives

    # Set massives mortality rates
    mr_massives::Vector{Float64} = cyclone_mr[:massives]
    cm_scens_massives = mr_massives[cyclone_cats].data
    for m in species[massives]
        cyclone_mortality[species=At(m)] .= cm_scens_massives
    end

    # Set mortality rates for branching corals at <= 5m depth
    below_5::BitVector = spatial_data.depth_med .<= -5
    if sum(below_5) > 0
        mr_bd5::Vector{Float64} = cyclone_mr[:branching_deeper_than_5]
        cm_scens_bd5::Array{Float64} = mr_bd5[cyclone_cats[location=below_5]].data
        for b in species[branchings]
            cyclone_mortality[locations=below_5, species=At(b)] .= cm_scens_bd5
        end
    end

    # Set mortality rates for branching corals at > 5m depth
    above_5::BitVector = spatial_data.depth_med .> -5
    if sum(above_5) > 0
        mr_bs5::Vector{Float64} = cyclone_mr[:branching_shallower_than_5]
        cm_scens_bs5::Array{Float64} = mr_bs5[cyclone_cats[location=above_5]].data
        for b in species[branchings]
            cyclone_mortality[locations=above_5, species=At(b)] .= cm_scens_bs5
        end
    end

    return cyclone_mortality
end

"""
    _cyclone_mortalities()

Mortality rates due to cyclones for each category (0 to 5) and for coral groups and depths.

# Notes
- Cyclone categories are represented through indices 1 to 6 of the arrays, there can be
distinct mortality rates for each coral functional group and depth.
- For Acropora (branching) the mortality depends on the depth (if it is deeper or shallower
than 5 meters).
- For massives the mortality does not depend on the depth.

# See also
https://github.com/open-AIMS/rrap-dg/blob/main/rrap_dg/cyclones/mortality_regression.jl
"""
function _cyclone_mortalities()::NamedTuple
    return (
        branching_deeper_than_5=[0.0, 0.0, 0.00993649, 0.497092, 0.989037, 0.99994],
        branching_shallower_than_5=[0.0, 0.0, 0.104957, 0.979102, 0.999976, 1.0],
        massives=[0.0, 0.0121482, 0.0155069, 0.0197484, 0.0246357, 0.0302982]
    )
end
