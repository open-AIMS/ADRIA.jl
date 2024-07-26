using ADRIA.decision:
    DecisionThresholds,
    DecisionWeights,
    DepthThresholds

using ADRIA.decision:
    SeedCriteriaWeights,
    FogCriteriaWeights

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
    site_data::DataFrame  # table of site data (depth, carrying capacity, etc)
    const site_id_col::String  # column to use as site ids, also used by the connectivity dataset (indicates order of `conn`)
    const cluster_id_col::String  # column of unique site ids
    init_coral_cover::YAXArray  # initial coral cover dataset
    const coral_growth::CoralGrowth  # coral
    const site_ids::Vector{String}  # Site IDs that are represented (i.e., subset of site_data[:, site_id_col], after missing sites are filtered)
    const removed_sites::Vector{String}  # indices of sites that were removed. Used to align site_data, DHW, connectivity, etc.
    dhw_scens::YAXArray  # DHW scenarios
    wave_scens::YAXArray  # wave scenarios
    cyclone_mortality_scens::Union{Matrix{<:Real},YAXArray}  # Cyclone mortality scenarios

    # Parameters
    model::Model  # core model
    sim_constants::SimConstants
end

"""
Barrier function to create Domain struct without specifying Intervention/Criteria/Coral/SimConstant parameters.
"""
function Domain(
    name::String,
    rcp::String,
    env_layers::EnvLayer,
    TP_base::YAXArray{T},
    site_data::DataFrame,
    site_id_col::String,
    cluster_id_col::String,
    init_coral_cover::YAXArray,
    coral_growth::CoralGrowth,
    site_ids::Vector{String},
    removed_sites::Vector{String},
    DHW::YAXArray,
    wave::YAXArray,
    cyclone_mortality::YAXArray,
)::ADRIADomain where {T<:Union{Float32,Float64}}
    sim_constants::SimConstants = SimConstants()
    criteria_weights::Vector{Union{DecisionWeights,DecisionThresholds}} = [
        SeedCriteriaWeights(),
        FogCriteriaWeights(),
        DepthThresholds()
    ]

    model::Model = Model((
        EnvironmentalLayer(DHW, wave, cyclone_mortality),
        Intervention(),
        criteria_weights...,
        Coral(),
    ))
    return ADRIADomain(
        name,
        rcp,
        env_layers,
        "",
        TP_base,
        site_data,
        site_id_col,
        cluster_id_col,
        init_coral_cover,
        coral_growth,
        site_ids,
        removed_sites,
        DHW,
        wave,
        cyclone_mortality,
        model,
        sim_constants,
    )
end

"""
    Domain(name::String, rcp::String, timeframe::Vector, site_data_fn::String, site_id_col::String, cluster_id_col::String, init_coral_fn::String, conn_path::String, dhw_fn::String, wave_fn::String, cyclone_mortality_fn::String)::Domain

Convenience constructor for Domain.

# Arguments
- `name` : Name of domain
- `dpkg_path` : location of data package
- `rcp` : RCP scenario represented
- `timeframe` : Time steps represented
- `site_data_fn` : File name of spatial data used
- `site_id_col` : Column holding name of reef the site is associated with (non-unique)
- `cluster_id_col` : Column holding unique site names/ids
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
    site_data_fn::String,
    site_id_col::String,
    cluster_id_col::String,
    init_coral_fn::String,
    conn_path::String,
    dhw_fn::String,
    wave_fn::String,
    cyclone_mortality_fn::String,
)::ADRIADomain
    local site_data::DataFrame
    try
        site_data = GDF.read(site_data_fn)
    catch err
        if !isfile(site_data_fn)
            error("Provided site data path is not valid or missing: $(site_data_fn).")
        else
            rethrow(err)
        end
    end

    if cluster_id_col ∉ names(site_data)
        @warn "Cluster ID column $(cluster_id_col) not found. Defaulting to UNIQUE_ID."
        cluster_id_col = "UNIQUE_ID"
    end
    if typeof(site_data[:, cluster_id_col][1]) != String
        site_data[!, cluster_id_col] .= string.(Int64.(site_data[:, cluster_id_col]))
    end

    env_layer_md::EnvLayer = EnvLayer(
        dpkg_path,
        site_data_fn,
        site_id_col,
        cluster_id_col,
        init_coral_fn,
        conn_path,
        dhw_fn,
        wave_fn,
        timeframe,
    )

    # Sort data to maintain consistent order
    sort!(site_data, Symbol[Symbol(site_id_col)])
    u_sids::Vector{String} = string.(collect(site_data[!, site_id_col]))
    # If site id column is missing then derive it from the Unique IDs
    if !in(site_id_col, names(site_data))
        site_data[!, site_id_col] .= String[d[2] for d in split.(u_sids, "_"; limit=2)]
    end

    site_data.row_id = 1:nrow(site_data)

    conn_ids::Vector{String} = u_sids
    connectivity::NamedTuple = site_connectivity(conn_path, u_sids)

    # Filter out missing entries
    site_data = site_data[coalesce.(in.(conn_ids, [connectivity.site_ids]), false), :]
    site_data.k .= site_data.k / 100.0  # Make `k` non-dimensional (provided as a percent)

    n_locs::Int64 = nrow(site_data)
    n_groups::Int64, n_sizes::Int64 = size(linear_extensions())
    coral_growth::CoralGrowth = CoralGrowth(n_locs, n_groups, n_sizes)
    n_group_and_size = coral_growth.n_group_and_size

    # Load initial coral cover relative to k area
    cover_params = ispath(init_coral_fn) ? (init_coral_fn,) : (n_group_and_size, n_locs)
    init_coral_cover = load_cover(cover_params...)

    dhw_params = ispath(dhw_fn) ? (dhw_fn, "dhw") : (timeframe, conn_ids)
    dhw = load_env_data(dhw_params...)

    waves_params = ispath(wave_fn) ? (wave_fn, "Ub") : (timeframe, conn_ids)
    waves = load_env_data(waves_params...)

    cyc_params = ispath(cyclone_mortality_fn) ? (cyclone_mortality_fn,) : (timeframe, site_data)
    cyclone_mortality = load_cyclone_mortality(cyc_params...)

    # Add compatability with non-migrated datasets but always default current coral spec
    if size(cyclone_mortality, 3) == 6
        n_groups = coral_growth.n_groups
        @warn """
        Cyclone mortality uses 6 functional groups. ADRIA uses $(n_groups).
        Skipping first functional group.
        """
        cyclone_mortality = cyclone_mortality[:, :, 2:end, :]
    end

    msg::String = "Provided time frame must match timesteps in DHW and wave data"
    msg = msg * "\n Got: $(length(timeframe)) | $(size(dhw, 1)) | $(size(waves, 1))"

    @assert length(timeframe) == size(dhw, 1) == size(waves, 1) msg

    return Domain(
        name,
        rcp,
        env_layer_md,
        connectivity.conn,
        site_data,
        site_id_col,
        cluster_id_col,
        init_coral_cover,
        coral_growth,
        connectivity.site_ids,
        connectivity.truncated,
        dhw,
        waves,
        cyclone_mortality,
    )
end

"""
    load_domain(ADRIADomain, path::String, rcp::String)::ADRIADomain
    load_domain(path::String, rcp::String)
    load_domain(path::String, rcp::Int64)

# Arguments
- `path` : location of data package
- `rcp` : RCP scenario to run. If none provided, no data path is set.
"""
function load_domain(::Type{ADRIADomain}, path::String, rcp::String)::ADRIADomain
    isdir(path) ? true : error("Path does not exist or is not a directory.")

    domain_name::String = basename(path)
    if length(domain_name) == 0
        domain_name = basename(dirname(path))
    end

    dpkg_details::Dict{String,Any} = _load_dpkg(path)

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

    site_id_col::String = "reef_siteid"
    cluster_id_col::String = "cluster_id"

    conn_path::String = joinpath(path, "connectivity/")
    spatial_path::String = joinpath(path, "spatial")

    gpkg_path::String = joinpath(spatial_path, "$(domain_name).gpkg")
    init_coral_cov::String = joinpath(spatial_path, "coral_cover.nc")

    dhw_fn::String = !isempty(rcp) ? joinpath(path, "DHWs", "dhwRCP$(rcp).nc") : ""
    wave_fn::String = !isempty(rcp) ? joinpath(path, "waves", "wave_RCP$(rcp).nc") : ""
    cyclone_mortality_fn::String = joinpath(path, "cyclones", "cyclone_mortality.nc")

    return Domain(
        domain_name,
        path,
        rcp,
        timeframe,
        gpkg_path,
        site_id_col,
        cluster_id_col,
        init_coral_cov,
        conn_path,
        dhw_fn,
        wave_fn,
        cyclone_mortality_fn,
    )
end
function load_domain(path::String, rcp::String)::ADRIADomain
    return load_domain(ADRIADomain, path, rcp)
end
function load_domain(path::String, rcp::Int64)::ADRIADomain
    return load_domain(ADRIADomain, path, string(rcp))
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
"""
function switch_RCPs!(d::ADRIADomain, RCP::String)::ADRIADomain
    @set! d.env_layer_md.DHW_fn = get_DHW_data(d, RCP)
    @set! d.env_layer_md.wave_fn = get_wave_data(d, RCP)
    @set! d.RCP = RCP

    @set! d.dhw_scens = load_env_data(d.env_layer_md.DHW_fn, "dhw")
    # @set! d.wave_scens = load_env_data(d.env_layer_md.wave_fn, "Ub")

    return d
end

function Base.show(io::IO, mime::MIME"text/plain", d::ADRIADomain)::Nothing
    df = model_spec(d)
    println("""
        Domain: $(d.name)

        Number of sites: $(n_locations(d))
        Site data file: $(d.env_layer_md.site_data_fn)
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
