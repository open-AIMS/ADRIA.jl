using ADRIA: SimConstants, Domain, GDF
using AxisKeys
using DataFrames
using DimensionalData
using Distributions
using MAT
using ModelParameters
using NamedDims
using NetCDF
using Statistics
using YAXArrays

using Infiltrator


mutable struct ReefModDomain <: Domain
    const name::String
    RCP::String
    env_layer_md
    scenario_invoke_time::String  # time latest set of scenarios were run
    const conn
    const in_conn
    const out_conn
    const strong_pred
    const site_data
    const site_id_col
    const cluster_id_col
    init_coral_cover
    const coral_growth::CoralGrowth
    const site_ids
    dhw_scens

    # `wave_scens` holds empty wave data to maintain compatibility with
    # ADRIA's dMCDA methods
    wave_scens

    # `cyclone_mortality_scens` holds dummy cyclone mortality data to maintain compatibility
    # with ADRIA's dMCDA methods
    cyclone_mortality_scens

    model
    sim_constants::SimConstants
    netcdf_dir::String
    loc_ids_path::String
end

"""
    load_domain(
        ::Type{ReefModDomain},
        netcdf_dir::String,
        geo_path::String,
        RCP::String
    )::ReefModDomain

Load ReefMod Matlab Dataset stored in netcdf file format. 
Uses a path ReefMod Engine data to fill missing required data

# Arguments
- `ReefModDomain` : DataType
- `netcdf_dir` : path to netcdf ReefMod Matlab Dataset Directory
- `geo_path` : path to geo data used in RMEDomain

# Returns
RMEDomain
"""
function load_domain(
    ::Type{ReefModDomain},
    netcdf_dir::String,
    rfmod_engine_path::String,
    RCP::String,
    timeframe =  (2022, 2100)
)::ReefModDomain
    netcdf_file = _find_netcdf(netcdf_dir, RCP)
    dom_dataset::Dataset = open_dataset(netcdf_file)
    data_files = joinpath(rfmod_engine_path, "data_files")
    loc_ids_path = joinpath(data_files, "dhw")
    loc_ids = _get_loc_ids(loc_ids_path)

    # force YAXArrays to load data into NamedDimsArray
    dhws = Cube(dom_dataset[["record_applied_DHWs"]])[timestep = At(timeframe[1] : timeframe[2])].data[:, :, :]
    dhw_scens = NamedDimsArray(
        dhws,
        timesteps=timeframe[1]:timeframe[2],
        locs=loc_ids,
        scenarios=1:size(dhws)[3]
    )
    
    site_data_path = joinpath(data_files, "region", "reefmod_gbr.gpkg")
    lat_col_id = "Y_COORD"
    lon_col_id = "X_COORD"
    site_data = GDF.read(site_data_path)
    site_id_col = "LOC_NAME_S"
    cluster_id_col = "LOC_NAME_S"
    site_ids = site_data[:, cluster_id_col]
    
    site_data[:, lat_col_id] = Cube(dom_dataset[["lat"]]).data
    site_data[:, lon_col_id] = Cube(dom_dataset[["lon"]]).data
    site_data[:, :area] = Cube(dom_dataset[["reef_area"]]).data
    site_data[:, :k] = dropdims(mean(Cube(dom_dataset[["nongrazable"]]), dims=2), dims=2).data

    # Convert non grazable area to a proportion of area and calculate grazable proportion
    site_data[:, :k] = (1 .- site_data[:, :k] ./ 100.0)
    # km^2 to m^2
    site_data[:, :area] *= 1e6 

    init_coral_cover = load_initial_cover(ReefModDomain, dom_dataset, site_data, loc_ids, timeframe[1])
    conn_data = load_connectivity(RMEDomain, data_files, loc_ids)
    in_conn, out_conn, strong_pred = ADRIA.connectivity_strength(
        conn_data, vec(site_data.area .* site_data.k), similar(conn_data)
    )

    # Set all site depths to 6m below sea level
    # (ReefMod does not account for depth)
    # Ensure column is of float type
    site_data[:, :depth_med] .= 6.0
    site_data[!, :depth_med] = convert.(Float64, site_data[!, :depth_med])

    # Add GBRMPA zone type info as well
    gbr_zt_path = joinpath(data_files, "region", "gbrmpa_zone_type.csv")
    gbr_zone_types = CSV.read(gbr_zt_path, DataFrame; types=String)
    missing_rows = ismissing.(gbr_zone_types[:, "GBRMPA Zone Types"])
    gbr_zone_types[missing_rows, "GBRMPA Zone Types"] .= ""
    zones = gbr_zone_types[:, "GBRMPA Zone Types"]
    zones = replace.(zones, "Zone" => "", " " => "")
    site_data[:, :zone_type] .= zones

    timeframe = (2022, 2100)

    # timesteps, location, scenario
    wave_scens = zeros(length(timeframe[1]:timeframe[2]), nrow(site_data), 1)
    wave_scens = NamedDimsArray(
        wave_scens;
        timesteps=timeframe[1]:timeframe[2],
        locs=loc_ids,
        scenarios=[1]
    )
    
    # current ReefMod mat data only contains cyclone classifications not mortality
    # timesteps, location, species, scenario
    cyc_scens = zeros(length(timeframe[1]:timeframe[2]), nrow(site_data), 6, 1)
    cyc_scens = NamedDimsArray(
        cyc_scens;
        timesteps=timeframe[1]:timeframe[2],
        locs=loc_ids,
        species=1:6,
        scenarios=[1]
    )

    env_md = EnvLayer(
        netcdf_dir,
        site_data_path,
        site_id_col,
        cluster_id_col,
        "",
        "",
        "",
        "",
        timeframe[1]:timeframe[2],
    )

    model::Model = Model((
        EnvironmentalLayer(dhw_scens, wave_scens, cyc_scens),
        Intervention(),
        CriteriaWeights(),
        Coral()
    ))

    return ReefModDomain(
        "ReefMod",
        RCP,
        env_md,
        "",
        conn_data,
        in_conn,
        out_conn,
        strong_pred,
        site_data,
        site_id_col,
        cluster_id_col,
        init_coral_cover,
        CoralGrowth(nrow(site_data)),
        site_ids,
        dhw_scens,
        wave_scens,
        cyc_scens,
        model,
        SimConstants(),
        netcdf_dir,
        loc_ids_path
    )
end

function _get_loc_ids(file_path::String)::Vector{String}
    possible_match = _get_relevant_files(file_path, "SSP")
    first_file = CSV.read(possible_match[1], DataFrame)
    loc_ids = String.(first_file[:, 1])
    return loc_ids
end

"""
    load_initial_cover(
        ::Type{ReefModDomain}, 
        dom_data::Dataset, 
        site_data::DataFrame, 
        loc_ids::Vector{String},
        init_yr::Int=2022
    )::NamedDimsArray
"""
function load_initial_cover(
    ::Type{ReefModDomain}, 
    dom_data::Dataset, 
    site_data::DataFrame, 
    loc_ids::Vector{String},
    init_yr::Int=2022
)::NamedDimsArray
    if !haskey(dom_data.cubes, :coral_cover_per_taxa)
        @error "coral_cover_per_taxa variable not found in YAXArrays dataset"
    end
    init_cc_per_taxa::YAXArray = Cube(dom_data[["coral_cover_per_taxa"]])[time = At(init_yr)]

    # the following class weight calculations are taken from ReefModEngine Domain calculation
    
    # Use ReefMod distribution for coral size class population (shape parameters have units log(cm^2))
    # as suggested by YM (pers comm. 2023-08-08 12:55pm AEST). Distribution is used to split ReefMod initial
    # species covers into ADRIA's 6 size classes by weighting with the cdf.
    reef_mod_area_dist = LogNormal(log(700), log(4))
    bin_edges_area = colony_mean_area(Float64[0, 2, 5, 10, 20, 40, 80])

    # Find integral density between bounds of each size class areas to create weights for each size class.
    cdf_integral = cdf.(reef_mod_area_dist, bin_edges_area)
    size_class_weights = (cdf_integral[2:end] .- cdf_integral[1:(end - 1)])
    size_class_weights = size_class_weights ./ sum(size_class_weights)

    # Take the mean over repeats, as suggested by YM (pers comm. 2023-02-27 12:40pm AEDT).
    # Convert from percent to relative values.
    # YAXArray ordering is [time ⋅ location ⋅ scenario]
    icc_data = ((dropdims(mean(init_cc_per_taxa; dims=:scenario); dims=:scenario)) ./ 100.0).data

    # Repeat species over each size class and reshape to give ADRIA compatible size (36 * n_locs).
    # Multiply by size class weights to give initial cover distribution over each size class.
    icc_data = Matrix(hcat(reduce.(vcat, eachrow(icc_data .* [size_class_weights]))...))

    # Convert values relative to absolute area to values relative to k area
    icc_data = _convert_abs_to_k(icc_data, site_data)

    n_species = length(init_cc_per_taxa[location=1, group=:, scenario=1])

    return NamedDimsArray(icc_data; species=1:(n_species * 6), locs=loc_ids)
end

function _find_netcdf(directory::String, scenario::String)::String
    pos_files = filter(isfile, readdir(directory, join=true))
    pos_files = filter(x -> occursin(".nc", x), pos_files)
    pos_files = filter(x -> occursin(scenario, x), pos_files)
    if length(pos_files) == 0
        ArgumentError("Unable to find NetCDF file relating to scenario: $(scenario)")
    elseif length(pos_files) > 1
        @info "Find multiple NetCDF files relating to scenario, using first"
    end
    return pos_files[1]
end

function switch_RCPs!(d::ReefModDomain, RCP::String)::ReefModDomain
    new_scen_fn = _find_netcdf(d.netcdf_dir, RCP)
    new_scen_dataset = open_dataset(new_scen_fn)
    dhws = Cube(new_scen_dataset[["record_applied_DHWs"]])[timestep = At(d.env_layer_md.timeframe)].data[:, :, :]
    scens = 1:size(dhws)[3]
    loc_ids = _get_loc_ids(d.loc_ids_path)
    d.dhw_scens = NamedDimsArray(
        dhws,
        timesteps=d.env_layer_md.timeframe,
        locs=loc_ids,
        scenarios=scens
    )
    return d
end
