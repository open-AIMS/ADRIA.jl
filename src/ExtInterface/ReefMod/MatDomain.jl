using ADRIA: SimConstants, Domain, GDF

using 
    Distributions
    Statistics

using
    AxisKeys
    DataFrames
    MAT
    ModelParameters
    NamedDims
    NetCDF
    YAXArrays

import ArchGDAL: createpoint

import YAXArrays.DD: At

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
end

"""
    load_domain(
        ::Type{ReefModDomain},
        fn_path::String,
        RCP::String
    )::ReefModDomain

Load ReefMod Matlab Dataset stored in netcdf file format. 
Uses a path ReefMod Engine data to fill missing required data

# Arguments
- `ReefModDomain` : DataType
- `netcdf_dir` : path to netcdf ReefMod Matlab Dataset Directory
- `RCP` : Representative Concentration Pathway scenario ID

# Returns
ReefModDomain
"""
function load_domain(
    ::Type{ReefModDomain},
    netcdf_dir::String,
    RCP::String;
    timeframe::Tuple=(2022, 2100)
)::ReefModDomain
    netcdf_file = _find_netcdf(fn_path, RCP)
    dom_dataset::Dataset = open_dataset(netcdf_file)

    site_ids = dom_dataset.properties["reef_ids"]

    # Force YAXArrays to load data into NamedDimsArray
    dhws = Cube(dom_dataset[["record_applied_DHWs"]])[timestep = At(timeframe[1] : timeframe[2])].data[:, :, :]
    dhw_scens = NamedDimsArray(
        dhws,
        timesteps=timeframe[1]:timeframe[2],
        locs=site_ids,
        scenarios=1:size(dhws)[3]
    )

    n_locations::Int = dom_dataset.properties["n_locations"]

    latitudes = Cube(dom_dataset[["lat"]]).data
    longitudes = Cube(dom_dataset[["lon"]]).data

    site_data = DataFrame((
        geom = [createpoint(x, y) for (x, y) in zip(longitudes, latitudes)],
        LOC_NAME_S = site_ids,
        UNIQUE_ID = site_ids,
        X_COORD = longitudes,
        Y_COORD = latitudes,
        depth_med = [Float64(6.0) for _ in 1:n_locations],
        area = Cube(dom_dataset[["reef_area"]]).data .* 1e6, # Convert to m^2 from km^2
        k = 1 .- dropdims(mean(Cube(dom_dataset[["nongrazable"]]), dims=2), dims=2).data ./ 100.0
    ))
    site_id_col = "LOC_NAME_S"
    cluster_id_col = "LOC_NAME_S"
    
    # Initial coral cover is loaded from the first year of reefmod 'coral_cover_per_taxa' data
    init_coral_cover = load_initial_cover(ReefModDomain, dom_dataset, site_data, site_ids, timeframe[1])

    # Connectivity data is retireved from a subdirectory because it's not contained in matfiles
    conn_data = load_connectivity(RMEDomain, netcdf_dir, site_ids)
    in_conn, out_conn, strong_pred = ADRIA.connectivity_strength(
        conn_data, vec(site_data.area .* site_data.k), similar(conn_data)
    )

    # GBRMPA zone types are not contained in matfiles
    site_data[:, :zone_type] .= ["" for _ in 1:n_locations]

    # timesteps, location, scenario
    wave_scens = zeros(length(timeframe[1]:timeframe[2]), n_locations, 1)
    wave_scens = NamedDimsArray(
        wave_scens;
        timesteps=timeframe[1]:timeframe[2],
        locs=site_ids,
        scenarios=[1]
    )
    
    # Current ReefMod mat data only contains cyclone classifications not mortality
    # timesteps, location, species, scenario
    cyc_scens = zeros(length(timeframe[1]:timeframe[2]), n_locations, 6, 1)
    cyc_scens = NamedDimsArray(
        cyc_scens;
        timesteps=timeframe[1]:timeframe[2],
        locs=site_ids,
        species=1:6,
        scenarios=[1]
    )

    env_md = EnvLayer(
        netcdf_dir,
        "",
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
    )
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
        @error "coral_cover_per_taxa variable not found in ReefMod data"
    end
    init_cc_per_taxa::YAXArray = Cube(dom_data[["coral_cover_per_taxa"]])[time = At(init_yr)]

    # The following class weight calculations are taken from ReefModEngine Domain calculation
    
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

"""
    _find_file(dir::String)::String
"""
function _find_file(dir::String, ident::Union{Regex, String})::String
    pos_files = filter(isfile, readdir(dir, join=true))
    pos_files = filter(x -> occursin(ident, x), pos_files)
    if length(pos_files) == 0
        ArgumentError("Unable to find file in $(dir)")
    elseif length(pos_files) > 1
        @info "Find multiple files matching identifier, using first"
    end
    return pos_files[1]
end

function _find_netcdf(dir::String, scenario::String)::String
    pos_files = filter(isfile, readdir(dir, join=true))
    pos_files = filter(x -> occursin(".nc", x), pos_files)
    pos_files = filter(x -> occursin(scenario, x), pos_files)
    if length(pos_files) == 0
        ArgumentError("Unable to find NetCDF file relating to scenario: $(scenario)")
    elseif length(pos_files) > 1
        @info "Find multiple NetCDF files relating to scenario, using first"
    end
    return pos_files[1]
end

"""
    switch_RCPs!(d::ReefModDomain, RCP::String)::ReefModDomain

Load different RCP into domain.
"""
function switch_RCPs!(d::ReefModDomain, RCP::String)::ReefModDomain
    new_scen_fn = _find_netcdf(d.env_layer_md.dpkg_path, RCP)
    new_scen_dataset = open_dataset(new_scen_fn)

    dhws = Cube(new_scen_dataset[["record_applied_DHWs"]])[timestep = At(d.env_layer_md.timeframe)].data[:, :, :]
    scens = 1:size(dhws)[3]
    loc_ids = d.site_ids

    d.dhw_scens = NamedDimsArray(
        dhws,
        timesteps=d.env_layer_md.timeframe,
        locs=loc_ids,
        scenarios=scens
    )
    return d
end

function Base.show(io::IO, mime::MIME"test/plain", d::ReefModDomain)::Nothing
    println("""
        ReefMod Domain: $(d.name)

        Number of site: $(n_locations(d))
    """)
    return nothing
end
