using ADRIA: GDF, ResultSet
using ArchGDAL: centroid
using DataFrames,
    DimensionalData,
    NetCDF,
    YAXArrays

struct RMEResultSet{T1, T2, D, G, D1} <: ResultSet
    name::String
    RCP::String

    site_ids::T1
    site_area::Vector{Float64}
    site_max_coral_cover::Vector{Float64}
    site_centroids::T2
    env_layer_md::EnvLayer
    connectivity_data::D
    site_data::G
    scenario_groups
    inputs::DataFrame

    sim_constants::D1
    model_spec::DataFrame

    outcomes::Dict{Symbol, YAXArray}
end

"""
    load_results(::Type{RMEResultSet}, data_dir::String)::RMEResultSet
    load_results(::Type{RMEResultSet}, data_dir::String, result_dir::String)::RMEResultSet

ReefMod Engine result interface.

# Arguments
- `data_dir`: Path to directory containing geospatial and connectivity data
- `result_dir`: Path to directory containing results in a netcdf file. Defaults to subdirectory of `data_dir` named "results"

# Returns
RMEResultSet struct compatible with most ADRIA analysis functionality.
"""
function load_results(
    ::Type{RMEResultSet}, data_dir::String
)::RMEResultSet
    return load_results(RMEResultSet, data_dir, joinpath(data_dir, "results"))
end
function load_results(
    ::Type{RMEResultSet}, data_dir::String, result_dir::String
)::RMEResultSet
    !isdir(data_dir) ? error("Expected a directory but received $(data_dir)") : nothing

    result_path = joinpath(result_dir, "results.nc")
    raw_set = open_dataset(result_path)

    # Name and RCP haven't been added to outputs yet.
    name::String = "ReefModEngine Result Set"
    netcdf_rcp::String = "Unknown"

    geodata_path = joinpath(data_dir, "region", "reefmod_gbr.gpkg")
    geodata = GDF.read(geodata_path)

    # Correct site ids in the same manner
    _standardize_cluster_ids!(geodata)

    reef_id_col = "UNIQUE_ID"
    cluster_id_col = "UNIQUE_ID"
    site_ids = geodata[:, reef_id_col]

    # Load accompanying ID list
    id_list_fn = _find_file(joinpath(data_dir, "id"), Regex("id_list.*.csv"))
    id_list = CSV.read(
        id_list_fn,
        DataFrame;
        header=false,
        comment="#",
    )

    _manual_id_corrections!(geodata, id_list)

    # Load reef area and convert from km^2 to m^2
    geodata[:, :area] = id_list[:, 2] .* 1e6

    # Calculate `k` area (1.0 - "ungrazable" area)
    geodata[:, :k] = 1 .- id_list[:, 3]

    # Connectivity is loaded using the same method as the Domain
    connectivity = load_connectivity(RMEDomain, data_dir, site_ids)

    location_max_coral_cover = 1 .- geodata.k ./ 100
    location_centroids = [centroid(multipoly) for multipoly ∈ geodata.geom]

    timeframe = 2008:2101
    if !haskey(raw_set.axes, :timesteps)
        @warn "Unable to find timestep axes. Defaulting to $(timeframe[1]):$(timeframe[end])"
    else
        timeframe = raw_set.timesteps[1]:raw_set.timesteps[end]
    end

    input_path = joinpath(result_dir, "scenarios.csv")
    inputs::DataFrame = _get_inputs(input_path)

    # Counterfactual scenario if outplant area and enrichment area are both 0
    scenario_groups::Dict{Symbol, BitVector} =
        _construct_scenario_groups(inputs, length(raw_set.scenarios))

    env_layer_md::EnvLayer = EnvLayer(
        data_dir,
        geodata_path,
        reef_id_col,
        cluster_id_col,
        "",
        joinpath(data_dir, "con_bin"),
        "",
        "",
        timeframe
    )

    outcomes::Dict{Symbol, YAXArray} = Dict();
    for (key, cube) in raw_set.cubes
        outcomes[key] = _reformat_cube(RMEResultSet, cube)
    end

    # Calculate relative cover
    n_locations::Int = length(raw_set.locations)
    outcomes[:relative_cover] =
        (outcomes[:total_cover] ./ 100) ./ reshape(geodata.k, (1, n_locations, 1))
    outcomes[:relative_taxa_cover] =
        (outcomes[:total_taxa_cover] ./ 100) ./ reshape(geodata.k, (1, n_locations, 1))

    mod_spec::DataFrame = _create_model_spec(inputs)

    return RMEResultSet(
        name,
        netcdf_rcp,
        site_ids,
        geodata.area,
        location_max_coral_cover,
        location_centroids,
        env_layer_md,
        connectivity,
        geodata,
        scenario_groups,
        inputs,
        SimConstants(),
        mod_spec,
        outcomes
    )
end

"""
    _get_inputs(filepath::String)

Read in scenario dataframe from given location. Warn and default to empty dataframe if not
found.
"""
function _get_inputs(filepath::String)::DataFrame
    if !isfile(filepath)
        @warn "Unable to find scenario spec at $(filepath). Skipping."
        return DataFrame()
    end
    return inputs = CSV.read(filepath, DataFrame, header=true)
end

"""
    _construct_scenario_groups(inputs::DataFrame, n_scenarios::Int)::Dict{Symbol, BitVector}

Using the inputs dataframe to determine if a run was counterfactual, construct the scenario
groupings for compatibility with ADRIA scenario plotting.
"""
function _construct_scenario_groups(
    inputs::DataFrame,
    n_scenarios::Int
)::Dict{Symbol, BitVector}
    if size(inputs) == (0, 0)
        return Dict(:counterfactual => BitVector([true for _ in 1:n_scenarios]))
    end

    counterfactual_scens::BitVector = BitVector([
        p_a == 0.0 && e_a == 0 for (p_a, e_a)
            in zip(inputs.outplant_area_pct, inputs.enrichment_area_pct)
    ])

    # Intervened if not counterfactual
    intervened_scens::BitVector = BitVector([
        p_a != 0 || e_a != 0 for (p_a, e_a)
            in zip(inputs.outplant_area_pct, inputs.enrichment_area_pct)
    ])

    scenario_groups::Dict{Symbol, BitVector} = Dict()

    # If the runs do not contain counterfactual or intervened runs exclude the key
    if any(counterfactual_scens)
        scenario_groups[:counterfactual] = counterfactual_scens
    end

    if any(intervened_scens)
        scenario_groups[:intervened] = intervened_scens
    end

    return scenario_groups
end

"""
    _reformat_cube(::Type{RMEResultSet}, cube::YAXArray)::YAXArray

Reformat an outcome cube to match dimension names and ordering of ADRIA cubes.
"""
function _reformat_cube(::Type{RMEResultSet}, cube::YAXArray)::YAXArray
    # Eventually ADRIA will use locations instead of sites
    if :locations in name.(cube.axes)
        cube = renameaxis!(cube, :locations => :sites)
    end

    return cube
end

"""
    _create_model_spec(scenario_spec::DataFrame)::DataFrame

Recreate a partial model specification to allow more analysis utilities.

Note: The model specification is incomplete as not all scenario data is extracted from
ReefModEngine. This partial specification will need to be updated as updates are made to
ReefModEngine.
"""
function _create_model_spec(scenario_spec::DataFrame)::DataFrame
    # Construct default model spec
    fieldname::Vector{Symbol} = [
        :dhw_tolerance,
        :outplant_count_per_m2,
        :outplant_area_pct,
        :n_outplant_locs,
        :enrichment_count_per_m2,
        :enrichment_area_pct,
        :n_enrichment_locs
    ]
    descriptions::Vector{String} = [
        "Extra DHW tolerance of thermally resistant outplanted coral",
        "Density of coral outplants over the entire deployment area",
        "Percentage of deployment area corals were deployed at",
        "Number of locations corals were outplanted at",
        "Density of coral enrichment over the entire deployment area",
        "Percentage of deployment area enrichment took place",
        "Number of locations enrichment took place at"
    ]
    field_names::Vector{String} = [
        "Restoration DHW Tolerance Outplants",
        "Outplant Density",
        "Outplant Area Percentage",
        "Outplant Location Count",
        "Enrichment Density",
        "Enrichment Area Percentage",
        "Enrichment Location Count"
    ]

    default_df = DataFrame(
        component="Intervention",
        fieldname=fieldname,
        description=descriptions,
        name=field_names
    )

    # Field names are the same as column names
    fieldnames::Vector{Symbol} = Symbol.(names(scenario_spec))
    inds::Vector{Int} = [findfirst(x -> x == fname, default_df.fieldname) for fname in fieldnames]

    return default_df[inds, :]
end

function Base.show(io::IO, mime::MIME"text/plain", rs::RMEResultSet)
    rcps = rs.RCP
    scens = length(rs.outcomes[:total_cover].scenarios)
    locations = length(rs.outcomes[:total_cover].sites)
    tf = rs.env_layer_md.timeframe

    println("""
    Name: $(rs.name)

    Results stored at: $(rs.env_layer_md.dpkg_path)

    RCP(s) represented: $(rcps)
    Scenarios run: $(scens)
    Number of locations: $(locations)
    Timesteps: $(tf)
    """)
end
