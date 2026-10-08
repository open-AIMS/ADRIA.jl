using CSV, DataFrames, ModelParameters
using ADRIA: SimConstants, Domain, GDF, EnvLayer, CoralGrowth, YAXArray, ZeroDataCube, DataCube, Model, EnvironmentalLayer, Intervention, Coral, GrowthAcceleration, CotsParams, COTSConnectivityForcing, adjust_sampling_bounds, load_env_data, load_initial_cover, location_connectivity, load_cots_connectivity_forcing, distance_matrix, nearest_neighbor_distances
using ADRIA.decision:
    DecisionThresholds,
    DecisionWeights,
    DepthThresholds,
    SeedCriteriaWeights,
    FogCriteriaWeights,
    MCCriteriaWeights

mutable struct LizardDomain <: Domain
    const name::String
    RCP::String
    env_layer_md::EnvLayer
    scenario_invoke_time::String
    const conn::YAXArray
    const cots_conn::YAXArray
    const cots_forcing::Union{COTSConnectivityForcing,Nothing}
    loc_data::DataFrame
    const loc_id_col::String
    const cluster_id_col::String
    init_coral_cover::YAXArray
    const coral_growth::CoralGrowth
    const loc_ids::Vector{String}
    const removed_locs::Vector{String}
    dhw_scens::YAXArray
    wave_scens::YAXArray
    cyclone_mortality_scens::Union{Matrix{<:Real},YAXArray}

    seed_target_locations::Vector{@NamedTuple{weight::Float64, target_locs::Vector{String}}}
    fog_target_locations::Vector{String}
    mc_target_locations::Vector{@NamedTuple{weight::Float64, target_locs::Vector{String}}}
    shade_target_locations::Vector{String}

    model::Model
    sim_constants::SimConstants
    cots_init_density::Union{Vector{Float64}, Nothing}
end

function load_domain(::Type{LizardDomain}, path::String, rcp::String)::LizardDomain
    isdir(path) || error("Path does not exist or is not a directory.")

    domain_name = basename(path)
    if isempty(domain_name)
        domain_name = basename(dirname(path))
    end

    dpkg_details = ADRIA._load_dpkg(path)

    # Hardcoded/assumed columns for Lizard domain
    location_id_col = "site_id"
    cluster_id_col = "UNIQUE_ID"

    spatial_path = joinpath(path, "spatial")
    gpkg_path = joinpath(spatial_path, "lizard_cluster.gpkg")
    loc_data = GDF.read(gpkg_path)

    ADRIA._standardise_location_columns!(loc_data; cluster_id_col=cluster_id_col, k_area_col="k", area_col="area")
    
    # Missing location_id_col? Let's assume site_id exists. If not, generate it.
    if !(location_id_col in names(loc_data))
        loc_data[!, location_id_col] = string.(1:nrow(loc_data))
    else
        loc_data[!, location_id_col] = string.(loc_data[!, location_id_col])
    end

    # Make loc_ids unique to match how they were saved in the connectivity matrix
    loc_ids = String.(DataFrames.make_unique(Symbol.(loc_data[!, location_id_col]); makeunique=true))
    loc_data[!, location_id_col] = loc_ids

    if !("CB_CALIB_GROUPS" in names(loc_data))
        loc_data[!, "CB_CALIB_GROUPS"] .= 1
    end

    dhw_path = joinpath(path, "DHWs", "dhw_RCP$(rcp).nc")
    dhw_raw = YAXArrays.Cube(dhw_path)
    # Rebuild YAXArray with correct labels
    dhw = DataCube(
        Float64.(Array(dhw_raw));
        timesteps=collect(dhw_raw.time),
        sites=loc_ids,
        scenarios=1:size(dhw_raw, 3)
    )
    timeframe = collect(dhw.timesteps)

    env_layer_md = EnvLayer(
        path,
        gpkg_path,
        location_id_col,
        cluster_id_col,
        joinpath(path, "initial_cover.nc"),
        joinpath(path, "connectivity"),
        joinpath(path, "DHWs", "dhw_RCP$(rcp).nc"),
        "",
        timeframe
    )

    connectivity = location_connectivity(joinpath(path, "connectivity"), loc_ids)
    cots_conn_root = joinpath(path, "cots_connectivity")
    cots_connectivity_dataset = get(ENV, "ADRIA_COTS_CONNECTIVITY_DATASET", "legacy")
    occursin(r"^[A-Za-z0-9_.-]+$", cots_connectivity_dataset) || error(
        "ADRIA_COTS_CONNECTIVITY_DATASET contains unsupported path characters"
    )
    cots_conn_path = if cots_connectivity_dataset == "legacy"
        cots_conn_root
    else
        joinpath(cots_conn_root, "datasets", cots_connectivity_dataset)
    end
    cots_connectivity_mode = lowercase(get(ENV, "ADRIA_COTS_CONNECTIVITY_MODE", "auto"))
    cots_connectivity_mode in ["auto", "cots", "coral"] || error(
        "ADRIA_COTS_CONNECTIVITY_MODE must be auto, cots, or coral; got '$cots_connectivity_mode'"
    )
    has_cots_connectivity = isdir(cots_conn_path) &&
        any(endswith(".csv"), readdir(cots_conn_path))
    if cots_connectivity_dataset != "legacy" &&
       cots_connectivity_mode != "coral" && !has_cots_connectivity
        error("Selected COTS dataset is missing its connectivity CSV: $cots_conn_path")
    end
    cots_connectivity = if cots_connectivity_mode != "coral" && has_cots_connectivity
        loaded = location_connectivity(cots_conn_path, loc_ids)
        isempty(loaded.truncated) || error(
            "COTS connectivity does not contain the complete Lizard site set: $(loaded.truncated)"
        )
        loaded
    elseif cots_connectivity_mode == "cots"
        error("COTS connectivity was required but no CSV was found in $cots_conn_path")
    else
        if cots_connectivity_mode == "auto"
            @warn "No COTS-specific connectivity found; using coral connectivity for COTS dispersal" path cots_connectivity_dataset
        end
        connectivity
    end
    cots_temporal_mode = lowercase(get(ENV, "COTS_CONNECTIVITY_TEMPORAL_MODE", "mean"))
    cots_temporal_mode in ["mean", "cycle", "sample"] || error(
        "COTS_CONNECTIVITY_TEMPORAL_MODE must be mean, cycle, or sample; got '$cots_temporal_mode'"
    )
    cots_forcing_dir = joinpath(cots_conn_path, "annual_forcing")
    cots_forcing = if isfile(joinpath(cots_forcing_dir, "provenance.toml"))
        load_cots_connectivity_forcing(cots_forcing_dir, loc_ids)
    else
        nothing
    end
    cots_temporal_mode == "mean" || !isnothing(cots_forcing) || error(
        "Annual COTS forcing is required for temporal mode '$cots_temporal_mode'; run build_lizard_cots_connectivity.jl"
    )
    
    # Drop truncated locations from loc_data to ensure matrix size alignment
    if length(connectivity.truncated) > 0
        invalid_ids = Set(connectivity.truncated)
        loc_data = loc_data[.!in.(loc_ids, Ref(invalid_ids)), :]
        loc_ids = loc_data[!, location_id_col]
    end
    
    # Distance matrix
    dist_matrix = distance_matrix(loc_data)
    loc_data.mean_to_neighbor .= nearest_neighbor_distances(dist_matrix, 10)

    coral_growth = CoralGrowth(nrow(loc_data))

    # initial_cover was already split into 36 bins in the domain builder
    icc_raw = YAXArrays.Cube(joinpath(path, "initial_cover.nc"))
    
    # We need species labels for the 36 bins
    n_sizes = coral_growth.n_sizes
    n_groups = coral_growth.n_groups
    
    # Generate labels using dummy cube
    dummy_cube = DataCube(
        zeros(n_groups, 1); 
        species=String.(collect(ADRIA.functional_group_names())), 
        locations=["1"]
    )
    cover_labels = ADRIA._cover_labels(dummy_cube, n_sizes, trues(n_groups))
    
    init_coral_cover = DataCube(
        Float64.(Array(icc_raw));
        species=cover_labels.species,
        locations=loc_ids
    )

    # dhw already loaded above
    
    # Wave and cyclone data are zeroes for now since we don't have them
    functional_groups = ADRIA.functional_group_names()
    wave_scens = ZeroDataCube(; T=Float64, timesteps=timeframe, locs=loc_ids, scenarios=[1])
    cyclone_mortality = ZeroDataCube(; T=Float64, timesteps=timeframe, locs=loc_ids, species=functional_groups, scenarios=[1])

    criteria_weights = [
        SeedCriteriaWeights(),
        FogCriteriaWeights(),
        MCCriteriaWeights(),
        DepthThresholds()
    ]

    cots_init = if "cots_density" in names(loc_data)
        [isnan(x) ? 0.0 : Float64(x) for x in coalesce.(loc_data.cots_density, 0.0)]
    else
        nothing
    end

    model = Model((
        EnvironmentalLayer(dhw, wave_scens, cyclone_mortality),
        Intervention(),
        criteria_weights...,
        Coral(),
        GrowthAcceleration(),
        CotsParams()
    ))

    return LizardDomain(
        domain_name,
        rcp,
        env_layer_md,
        "",
        connectivity.conn,
        cots_connectivity.conn,
        cots_forcing,
        loc_data,
        location_id_col,
        cluster_id_col,
        init_coral_cover,
        coral_growth,
        connectivity.loc_ids,
        connectivity.truncated,
        dhw,
        wave_scens,
        cyclone_mortality,
        [(weight=1.0, target_locs=loc_ids)],
        loc_ids,
        [(weight=1.0, target_locs=loc_ids)],
        loc_ids,
        model,
        SimConstants(),
        cots_init
    )
end

function ADRIA.switch_RCPs!(d::LizardDomain, RCP::String)::LizardDomain
    return d
end

