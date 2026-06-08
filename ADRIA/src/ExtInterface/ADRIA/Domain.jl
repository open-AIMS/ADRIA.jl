using ADRIA.decision:
    DecisionThresholds,
    DecisionWeights,
    DepthThresholds

using ADRIA.decision:
    SeedCriteriaWeights,
    FogCriteriaWeights,
    MCCriteriaWeights

"""
    ADRIADomain{Σ,M,I,D,X,Y,Z}

Core ADRIA domain. Represents study area.
"""
mutable struct ADRIADomain <: Domain
    const name::String  # human-readable name
    RCP::String  # RCP scenario represented
    env_layer_md::EnvLayer  # Layers used
    scenario_invoke_time::String  # time latest set of scenarios were run
    const conn::YAXArray  # connectivity data
    loc_data::DataFrame  # table of location data (depth, carrying capacity, etc)
    const loc_id_col::String  # column to use as location ids, also used by the connectivity dataset (indicates order of `conn`)
    const cluster_id_col::String  # column of unique cluster ids
    init_coral_cover::YAXArray  # initial coral cover dataset
    const coral_growth::CoralGrowth  # coral
    const loc_ids::Vector{String}  # Location IDs that are represented (i.e., subset of loc_data[:, location_id_col], after missing locations are filtered)
    const removed_locs::Vector{String}  # indices of locations that were removed. Used to align loc_data, DHW, connectivity, etc.
    dhw_scens::YAXArray  # DHW scenarios
    wave_scens::YAXArray  # wave scenarios
    cyclone_mortality_scens::Union{Matrix{<:Real},YAXArray}  # Cyclone mortality scenarios

    # Strategy target locations
    # Each element of these vector is a pair weight and list of location ids.
    # The weights represents the share of the intervention capacity that will be distributed
    # to each list of locations
    seed_target_locations::Vector{@NamedTuple{weight::Float64,target_locs::Vector{String}}}  # locations eligible for seeding
    fog_target_locations::Vector{String}   # locations eligible for fogging
    mc_target_locations::Vector{@NamedTuple{weight::Float64,target_locs::Vector{String}}}  # locations eligible for moving corals
    shade_target_locations::Vector{String}    # locations eligible for shading

    # Parameters
    model::Model  # core model
    sim_constants::SimConstants
end

"""
    _coral_calib_overrides(nc_ds)::Dict{String,Float64}

Read calibrated coral parameter values from a NetCDF dataset (produced by CoralBlox
calibration) and return a Dict mapping Coral struct field names to their calibrated values.
Covers `linear_extension`, `mb_rate`, `dist_mean`, `linear_extension_scale`, and
`mb_rate_scale`. Any field not present in the dataset falls back to the ADRIA default.
"""
function _coral_calib_overrides(nc_ds)::Dict{String,Float64}
    overrides = Dict{String,Float64}()
    fg_names = string.(functional_group_names())

    for (param_name, varname) in (
        ("linear_extension", "linear_extension"),
        ("mb_rate", "mb_rate"),
        ("dist_mean", "dist_mean")
    )
        data = Array(nc_ds[varname])  # (n_groups, n_sizes)
        for (fg_idx, fg) in enumerate(fg_names), sc = 1:size(data, 2)
            overrides["$(fg)_$(fg_idx)_$(sc)_$(param_name)"] = data[fg_idx, sc]
        end
    end

    le_scale_da = nc_ds["linear_extension_scale"]
    mb_scale_da = nc_ds["mb_rate_scale"]
    bg_ids = collect(DimensionalData.lookup(le_scale_da, :cb_calib_group))
    le_scale = Array(le_scale_da)
    mb_scale = Array(mb_scale_da)

    for (fg_idx, fg) in enumerate(fg_names), (bg_idx, bg) in enumerate(bg_ids)
        overrides["linear_extension_scale_cb_group_$(bg)_$(fg)"] = le_scale[fg_idx, bg_idx]
        overrides["mb_rate_scale_cb_group_$(bg)_$(fg)"] = mb_scale[fg_idx, bg_idx]
    end

    return overrides
end

"""
    _growth_accel_calib_overrides(nc_ds)::Dict{String,Float64}

Read calibrated growth acceleration parameter values from a NetCDF dataset and return a
Dict mapping GrowthAcceleration struct field names to their calibrated values.
"""
function _growth_accel_calib_overrides(nc_ds)::Dict{String,Float64}
    overrides = Dict{String,Float64}()
    ga_da = nc_ds["growth_acceleration"]  # (cb_calib_group=12, accel_param=3)
    bg_ids = collect(DimensionalData.lookup(ga_da, :cb_calib_group))
    ap_vals = collect(DimensionalData.lookup(ga_da, :accel_param))
    ga = Array(ga_da)

    for (bg_idx, bg) in enumerate(bg_ids), (ap_idx, ap) in enumerate(ap_vals)
        overrides["growth_acceleration_cb_group_$(bg)_$(ap)"] = ga[bg_idx, ap_idx]
    end

    return overrides
end

"""
Barrier function to create Domain struct without specifying Intervention/Criteria/Coral/SimConstant parameters.
"""
function Domain(
    name::String,
    rcp::String,
    env_layers::EnvLayer,
    TP_base::YAXArray{T},
    location_data::DataFrame,
    location_id_col::String,
    cluster_id_col::String,
    init_coral_cover::YAXArray,
    coral_growth::CoralGrowth,
    location_ids::Vector{String},
    removed_locations::Vector{String},
    DHW::YAXArray,
    wave::YAXArray,
    cyclone_mortality::YAXArray;
    calib_params_fn::String=""
)::ADRIADomain where {T<:Union{Float32,Float64}}
    sim_constants::SimConstants = SimConstants()

    if has_mcb_scenarios(DHW)
        albedos = collect(DHW.albedo)
        durations = collect(DHW.mcb_durations)
    else
        albedos = [0.0, 0.0]
        durations = [0.0, 0.0]
    end

    intervention_params = (
        mcb_albedo=Factor(
            albedos[1];
            ptype="ordered categorical",
            dist=CategoricalDistribution,
            dist_params=(Tuple(albedos)),
            name="MCB Albedo",
            description="Albedo level to use from 5D DHW dataset."
        ),
        mcb_duration=Factor(
            durations[1];
            ptype="ordered categorical",
            dist=CategoricalDistribution,
            dist_params=(Tuple(durations)),
            name="MCB Duration",
            description="Duration level (yearly days) to use from 5D DHW dataset."
        )
    )

    local coral_instance::Coral
    local growth_accel_instance::GrowthAcceleration

    if !isempty(calib_params_fn) && isfile(calib_params_fn)
        @info "Loading calibrated coral parameters from $(calib_params_fn)"

        nc_ds = open_dataset(calib_params_fn)
        coral_overrides, growth_accel_overrides = _coral_calib_overrides(nc_ds),
        _growth_accel_calib_overrides(nc_ds)

        coral_instance = create_coral_instance(; overrides=coral_overrides)
        growth_accel_instance = create_growth_acceleration_instance(;
            overrides=growth_accel_overrides
        )
    else
        coral_instance = Coral()
        growth_accel_instance = GrowthAcceleration()
    end

    model::Model = Model((
        EnvironmentalLayer(DHW, wave, cyclone_mortality),
        Intervention(; intervention_params...),
        SeedCriteriaWeights(),
        FogCriteriaWeights(),
        MCCriteriaWeights(),
        DepthThresholds(),
        coral_instance,
        growth_accel_instance,
        CotsParams()
    ))

    return ADRIADomain(
        name,
        rcp,
        env_layers,
        "",
        TP_base,
        location_data,
        location_id_col,
        cluster_id_col,
        init_coral_cover,
        coral_growth,
        location_ids,
        removed_locations,
        DHW,
        wave,
        cyclone_mortality,
        [(weight=1.0, target_locs=location_ids)],
        location_ids,
        [(weight=1.0, target_locs=location_ids)],
        location_ids,
        model,
        sim_constants
    )
end

"""
    Domain(name::String, rcp::String, timeframe::Vector, location_data_fn::String, location_id_col::String, cluster_id_col::String, init_coral_fn::String, conn_path::String, dhw_fn::String, wave_fn::String, cyclone_mortality_fn::String)::Domain

Convenience constructor for Domain.

# Arguments
- `name` : Name of domain
- `dpkg_path` : location of data package
- `rcp` : RCP scenario represented
- `timeframe` : Time steps represented
- `location_data_fn` : File name of spatial data used
- `location_id_col` : Column holding name of reef the location is associated with (non-unique)
- `cluster_id_col` : Column holding unique cluster names/ids
- `init_coral_fn` : Name of file holding initial coral cover values
- `conn_path` : Path to directory holding connectivity data
- `dhw_fn` : Filename of DHW data cube in use
- `wave_fn` : Filename of wave data cube
- `cyclone_mortality_fn` : Filename of cyclone mortality data cube
"""
function Domain(
    name::String,
    dpkg_path::String,
    rcp::String,
    timeframe::Vector,
    location_data_fn::String,
    location_id_col::String,
    cluster_id_col::String,
    k_area_col::String,
    area_col::String,
    init_coral_fn::String,
    conn_path::String,
    dhw_fn::String,
    wave_fn::String,
    cyclone_mortality_fn::String;
    calib_params_fn::String=""
)::ADRIADomain
    location_data = load_location_data(location_data_fn)

    _validate_Domain(
        location_data, location_id_col, cluster_id_col, k_area_col, area_col, init_coral_fn,
        dhw_fn
    )

    _standardise_location_columns!(
        location_data;
        cluster_id_col=cluster_id_col, k_area_col=k_area_col, area_col=area_col
    )

    env_layer_md::EnvLayer = EnvLayer(
        dpkg_path,
        location_data_fn,
        location_id_col,
        cluster_id_col,
        init_coral_fn,
        conn_path,
        dhw_fn,
        wave_fn,
        timeframe
    )

    # Sort data to maintain consistent order
    sort!(location_data, Symbol[Symbol(location_id_col)])

    # Ensure location_id_col values are strings
    if !all(isa.(location_data[!, location_id_col], String))
        location_data[!, location_id_col] = string.(location_data[!, location_id_col])
    end

    u_sids::Vector{String} = location_data[!, location_id_col]

    # If location id column is missing then derive it from the Unique IDs
    if !in(location_id_col, names(location_data))
        # Extract ID component from unique ID (second part after split on "_"; limit=2)
        location_data[!, location_id_col] = [
            last(split(id, "_"; limit=2)) for id in u_sids
        ]
    end

    location_data.row_id = 1:nrow(location_data)

    conn_ids::Vector{String} = u_sids
    connectivity::NamedTuple = location_connectivity(conn_path, u_sids)

    # Filter out missing entries
    location_data = location_data[
        coalesce.(in.(conn_ids, [connectivity.loc_ids]), false), (:)
    ]

    # Make `k` non-dimensional (if it is a percentage)
    if any(x -> x > 1.0, location_data.k)
        @info "Values in 'k' are > 1.0, assuming they are percentages and dividing by 100."
        location_data.k .= location_data.k / 100.0
    end

    dist_matrix = distance_matrix(location_data)
    location_data.mean_to_neighbor .= nearest_neighbor_distances(dist_matrix, 10)

    n_locs::Int64 = nrow(location_data)
    coral_growth::CoralGrowth = CoralGrowth(n_locs)

    # Load initial coral cover relative to k area
    init_coral_cover = load_initial_cover(init_coral_fn)

    dhw::YAXArray{Float64} = load_env_data(dhw_fn, "dhw")
    waves = load_wave_data(wave_fn, timeframe, conn_ids)
    cyclone_mortality = load_cyclone_data(cyclone_mortality_fn, timeframe, u_sids)

    # If environmental data is over different timeframes align them.
    if !all(dhw.timesteps .== 1:length(dhw.timesteps))
        dhw = dhw[timesteps = At(timeframe)]
    end
    if !all(waves.timesteps .== 1:length(waves.timesteps))
        waves = waves[timesteps = At(timeframe)]
    end
    if !all(cyclone_mortality.timesteps .== 1:length(cyclone_mortality.timesteps))
        cyclone_mortality = cyclone_mortality[timesteps = At(timeframe)]
    end

    msg::String =
        "Provided time frame must match timesteps in DHW and wave data\n" *
        "Got: $(length(timeframe)) | $(size(dhw, 1)) | $(size(waves, 1))"

    @assert length(timeframe) == size(dhw, 1) == size(waves, 1) msg

    return Domain(
        name,
        rcp,
        env_layer_md,
        connectivity.conn,
        location_data,
        location_id_col,
        cluster_id_col,
        init_coral_cover,
        coral_growth,
        connectivity.loc_ids,
        connectivity.truncated,
        dhw,
        waves,
        cyclone_mortality;
        calib_params_fn=calib_params_fn
    )
end

function _validate_Domain(
    location_data::DataFrame,
    location_id_col::String,
    cluster_id_col::String,
    k_area_col::String,
    area_col::String,
    init_coral_fn::String,
    dhw_fn::String
)::Nothing
    # Ensure other required columns are present
    required_cols = [location_id_col, cluster_id_col, k_area_col, area_col]
    existence_mask = [col in names(location_data) for col in required_cols]
    if !all(existence_mask)
        error("Columns '$(required_cols[existence_mask])' not found in location data.")
    end

    if !ispath(init_coral_fn)
        error(
            "Initial coral cover data file not found at $(init_coral_fn). Initial coral cover data is required."
        )
    end

    if !ispath(dhw_fn)
        error("DHW data file not found at $(dhw_fn). DHW data is always required.")
    end

    # Ensure other required columns are present
    required_cols = [location_id_col, cluster_id_col, k_area_col, area_col]
    existence_mask = [col in names(location_data) for col in required_cols]
    if !all(existence_mask)
        error("Columns '$(required_cols[existence_mask])' not found in location data.")
    end

    return nothing
end

"""
    _standardise_location_columns!(location_data::DataFrame; cluster_id_col::Union{Nothing, String}=nothing, k_area_col::Union{Nothing, String}=nothing, area_col::Union{Nothing, String}=nothing,)::Nothing

Rename the columns of a given locationd dataframe to intenral ADRIA standards.
ADRIA Domains assume specific column names for specific features in the location data.

|  Feature | Internal Standard |
| ------------- |:-------------:|
| Habitable Space | "k" |
| Cluster ID | "cluster_id" |
| Location Area | "area" |
"""
function _standardise_location_columns!(
    location_data::DataFrame;
    cluster_id_col::Union{Nothing,String}=nothing,
    k_area_col::Union{Nothing,String}=nothing,
    area_col::Union{Nothing,String}=nothing
)::Nothing
    if !isnothing(cluster_id_col)
        location_data[!, :cluster_id] = location_data[!, Symbol(cluster_id_col)]
    end
    if !isnothing(k_area_col)
        location_data[!, :k] = location_data[!, Symbol(k_area_col)]
    end
    if !isnothing(area_col)
        location_data[!, :area] = location_data[!, Symbol(area_col)]
    end

    return nothing
end

function _get_spatial_column_names(dpkg_details::Dict{String,Any})::Dict{String,String}
    # Default column names
    col_names = Dict(
        "location_id_col" => "reef_siteid",
        "cluster_id_col" => "cluster_id",
        "k_col" => "k",
        "area_col" => "area"
    )

    if haskey(dpkg_details, "resources")
        for resource in dpkg_details["resources"]
            if resource["name"] == "spatial_data"
                for (key, val) in col_names
                    col_names[key] = get(resource, key, val)
                end
                break # Found the spatial data resource, no need to check others
            end
        end
    end

    return col_names
end

"""
    load_domain(ADRIADomain, path::String, rcp::String)::ADRIADomain
    load_domain(path::String, rcp::String)
    load_domain(path::String, rcp::Int64)

# Arguments
- `path` : location of data package
- `rcp` : RCP scenario to run. If none provided, no data path is set.
"""
function load_domain(
    ::Type{ADRIADomain}, path::String, rcp::String; calib_params_fn::String=""
)::ADRIADomain
    isdir(path) ? true : error("Path does not exist or is not a directory.")

    domain_name::String = basename(path)
    if length(domain_name) == 0
        domain_name = basename(dirname(path))
    end

    dpkg_details::Dict{String,Any} = _load_dpkg(path)

    # Extract spatial column names from datapackage.json
    spatial_col_names = _get_spatial_column_names(dpkg_details)
    location_id_col::String = spatial_col_names["location_id_col"]
    cluster_id_col::String = spatial_col_names["cluster_id_col"]
    k_area_col::String = spatial_col_names["k_col"]
    area_col::String = spatial_col_names["area_col"]

    # Handle compatibility
    # Extract the time frame represented in this data package
    md_timeframe::Tuple{Int64,Int64} = Tuple(
        dpkg_details["simulation_metadata"]["timeframe"]
    )

    if length(md_timeframe) == 2
        @assert md_timeframe[1] < md_timeframe[2] "Start date/year specified in data package must be < end date/year"
        # If only two elements, assume a range is specified.
        # Collate the time steps as a full list if necessary
        timeframe::Vector{Int64} = collect(md_timeframe[1]:md_timeframe[2])
    else
        # Otherwise assume entry specifies yearly timesteps
        timeframe = parse.(Int64, md_timeframe)
    end

    conn_path::String = joinpath(path, "connectivity/")
    spatial_path::String = joinpath(path, "spatial")

    gpkg_path::String = joinpath(spatial_path, "$(domain_name).gpkg")
    init_coral_cov::String = joinpath(spatial_path, "coral_cover.nc")

    dhw_fn::String = joinpath(path, "DHWs", "dhwRCP$(rcp).nc")
    wave_fn::String = !isempty(rcp) ? joinpath(path, "waves", "wave_RCP$(rcp).nc") : ""
    cyclone_mortality_fn::String = joinpath(path, "cyclones", "cyclone_mortality.nc")

    return Domain(
        domain_name,
        path,
        rcp,
        timeframe,
        gpkg_path,
        location_id_col,
        cluster_id_col,
        k_area_col,
        area_col,
        init_coral_cov,
        conn_path,
        dhw_fn,
        wave_fn,
        cyclone_mortality_fn;
        calib_params_fn=calib_params_fn
    )
end
function load_domain(path::String, rcp::String; calib_params_fn::String="")::ADRIADomain
    return load_domain(ADRIADomain, path, rcp; calib_params_fn=calib_params_fn)
end
function load_domain(path::String, rcp::Int64; calib_params_fn::String="")::ADRIADomain
    return load_domain(ADRIADomain, path, string(rcp); calib_params_fn=calib_params_fn)
end

"""Get the path to the DHW data associated with the domain."""
function get_DHW_data(d::ADRIADomain, RCP::String)::String
    return joinpath(d.env_layer_md.dpkg_path, "DHWs", "dhwRCP$(RCP).nc")
end

"""Get the path to the wave data associated with the domain."""
function get_wave_data(d::ADRIADomain, RCP::String)::String
    return joinpath(d.env_layer_md.dpkg_path, "waves", "wave_RCP$(RCP).nc")
end

"""
    switch_RCPs!(d::Domain, RCP::String)::Domain

Switch environmental datasets to represent the given RCP.

# Examples
```julia
dom = ADRIA.load_domain(<path to domain>, "45")

@info dom.RCP
# "45"

ADRIA.switch_RCPs!(dom, "85")
@info dom.RCP
# "85"
```
"""
function switch_RCPs!(d::ADRIADomain, RCP::String)::ADRIADomain
    d.env_layer_md.DHW_fn = get_DHW_data(d, RCP)
    d.env_layer_md.wave_fn = get_wave_data(d, RCP)
    d.RCP = RCP

    d.dhw_scens = load_env_data(d.env_layer_md.DHW_fn, "dhw")

    # Constraint dhw timeframe to timeframe in datapackage
    if !all(d.dhw_scens.timesteps .== 1:length(d.dhw_scens.timesteps))
        d.dhw_scens = d.dhw_scens[timesteps = At(
            d.env_layer_md.timeframe
        )]
    end

    return d
end

function Base.show(io::IO, mime::MIME"text/plain", d::ADRIADomain)::Nothing
    df = model_spec(d)
    println("""
        Domain: $(d.name)

        Number of locations: $(n_locations(d))
        Location data file: $(d.env_layer_md.loc_data_fn)
        Connectivity file: $(d.env_layer_md.connectivity_fn)
        DHW file: $(d.env_layer_md.DHW_fn)
        Wave file: $(d.env_layer_md.wave_fn)
        Cyclone mortality file:
        Timeframe: $(d.env_layer_md.timeframe[1]) - $(d.env_layer_md.timeframe[end])

        Model Specification:
        - Parameters: $(nrow(df))
        - Number of constants: $(nrow(df[df.is_constant .== true, :]))
        """)

    # println("\nEcosystem model specification:")
    # show(io, mime, model_spec(d)[:, [:component, :fieldname, :val, :full_bounds, :dists, :is_constant]])
    return nothing
end
