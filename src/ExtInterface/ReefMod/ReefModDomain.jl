using ADRIA: SimConstants, Domain, GDF, DataCube

using Distributions, Statistics

using DataFrames,
    MAT,
    ModelParameters,
    NetCDF,
    YAXArrays

import ArchGDAL: createpoint

import YAXArrays.DD: At


mutable struct ReefModDomain <: AbstractReefModDomain
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
- `fn_path` : path to netcdf ReefMod Matlab Dataset Directory
- `RCP` : Representative Concentration Pathway scenario ID

# Returns
ReefModDomain
"""
function load_domain(
    ::Type{ReefModDomain},
    fn_path::String,
    RCP::String;
    timeframe::Tuple=(2022, 2100)
)::ReefModDomain
    netcdf_file = _find_netcdf(fn_path, RCP)
    dom_dataset::Dataset = open_dataset(netcdf_file)

    # Read location data
    geodata_dir = joinpath(fn_path, "region")
    geodata_fn = _find_file(geodata_dir, ".gpkg")
    spatial_data = GDF.read(geodata_fn)

    _standardize_cluster_ids!(spatial_data)

    reef_id_col = "UNIQUE_ID"
    cluster_id_col = "UNIQUE_ID"
    site_ids = spatial_data[:, reef_id_col]

    # Load accompanying ID list
    # TODO: Create canonical geopackage file that aligns all IDs.
    #       Doing this removes the need for the manual correction below and removes the
    #       dependency on this file.
    id_list_fn = _find_file(joinpath(fn_path, "id"), Regex("id_list.*.csv"))
    id_list = CSV.read(
        id_list_fn,
        DataFrame;
        header=false,
        comment="#",
    )

    _manual_id_corrections!(spatial_data, id_list)

    # Load reef area and convert from km^2 to m^2
    spatial_data[:, :area] = id_list[:, 2] .* 1e6

    # Calculate `k` area (1.0 - "ungrazable" area)
    spatial_data[:, :k] = 1 .- id_list[:, 3]

    # Load DHWs
    dhws = Cube(
        dom_dataset[["record_applied_DHWs"]]
    )[timestep=At(timeframe[1]:timeframe[2])]

    # Redfine dimensions as ReefMod Matfiles do not contain reef ids.
    # Forcibly load data as disk arrays are not fully support.
    dhw_scens = DataCube(
        dhws.data[:, :, :];
        timesteps=timeframe[1]:timeframe[2],
        locs=site_ids,
        scenarios=1:size(dhws)[3]
    )

    # Initial coral cover is loaded from the first year of reefmod 'coral_cover_per_taxa' data
    init_coral_cover = load_initial_cover(
        ReefModDomain, dom_dataset, spatial_data, site_ids, timeframe[1]
    )

    # Connectivity data is retireved from a subdirectory because it's not contained in matfiles
    conn_data = load_connectivity(RMEDomain, fn_path, site_ids)
    in_conn, out_conn, strong_pred = ADRIA.connectivity_strength(
        conn_data, vec(spatial_data.area .* spatial_data.k), similar(conn_data)
    )

    spatial_data[:, :depth_med] .= 6.0
    spatial_data[!, :depth_med] = convert.(Float64, spatial_data[!, :depth_med])
    # GBRMPA zone types are not contained in matfiles
    spatial_data[:, :zone_type] .= ["" for _ in 1:nrow(spatial_data)]

    # timesteps, location, scenario
    wave_scens = ZeroDataCube(;
        T=Float64,
        timesteps=timeframe[1]:timeframe[2],
        locs=site_ids,
        scenarios=[1]
    )

    env_md = EnvLayer(
        fn_path,
        geodata_fn,
        reef_id_col,
        cluster_id_col,
        "",
        "",
        "",
        "",
        timeframe[1]:timeframe[2],
    )

    criteria_weights::Vector{Union{DecisionWeights,DecisionThresholds}} = [
        SeedCriteriaWeights(),
        FogCriteriaWeights(),
        DepthThresholds()
    ]

    cyclone_mortality_scens::YAXArray{Float64,4} = _cyclone_mortality_scens(
        dom_dataset,
        spatial_data,
        site_ids,
        timeframe
    )

    model::Model = Model((
        EnvironmentalLayer(dhw_scens, wave_scens, cyclone_mortality_scens),
        Intervention(),
        criteria_weights...,
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
        spatial_data,
        reef_id_col,
        cluster_id_col,
        init_coral_cover,
        CoralGrowth(nrow(spatial_data)),
        site_ids,
        dhw_scens,
        wave_scens,
        cyclone_mortality_scens,
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
    )::YAXArray
"""
function load_initial_cover(
    ::Type{ReefModDomain},
    dom_data::Dataset,
    site_data::DataFrame,
    loc_ids::Vector{String},
    init_yr::Int=2022
)::YAXArray
    if !haskey(dom_data.cubes, :coral_cover_per_taxa)
        @error "coral_cover_per_taxa variable not found in ReefMod data"
    end
    init_cc_per_taxa::YAXArray = Cube(dom_data[["coral_cover_per_taxa"]])[time=At(init_yr)]

    # The following class weight calculations are taken from ReefModEngine Domain calculation

    # Use ReefMod distribution for coral size class population (shape parameters have units log(cm^2))
    # as suggested by YM (pers comm. 2023-08-08 12:55pm AEST). Distribution is used to split ReefMod initial
    # species covers into ADRIA's 6 size classes by weighting with the cdf.
    reef_mod_area_dist = LogNormal(log(700), log(4))
    bin_edges_area = colony_mean_area(Float64[0, 2, 5, 10, 20, 40, 80])

    # Find integral density between bounds of each size class areas to create weights for each size class.
    cdf_integral = cdf.(reef_mod_area_dist, bin_edges_area)
    size_class_weights = (cdf_integral[2:end] .- cdf_integral[1:(end-1)])
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

    return DataCube(icc_data; species=1:(n_species*6), locs=loc_ids)
end

"""
    _find_file(dir::String)::String
"""
function _find_file(dir::String, ident::Union{Regex,String})::String
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

    dhws = Cube(
        new_scen_dataset[["record_applied_DHWs"]]
    )[timestep=At(d.env_layer_md.timeframe)].data[:, :, :]

    scens = 1:size(dhws)[3]
    loc_ids = d.site_ids

    d.dhw_scens = DataCube(
        dhws,
        timesteps=d.env_layer_md.timeframe,
        locs=loc_ids,
        scenarios=scens
    )
    return d
end

function Base.show(io::IO, mime::MIME"text/plain", d::ReefModDomain)::Nothing
    println("""
        ReefMod Domain: $(d.name)

        Number of locations: $(n_locations(d))
    """)

    return nothing
end

"""
    _cyclone_mortality_scens(dom_dataset, spatial_data, site_ids, timeframe)::YAXArray{Float64}

Cyclone scenarios (from 0 to 5) are converted to mortality rates. More details of how this
happens van be found in https://github.com/open-AIMS/rrap-dg.
"""
function _cyclone_mortality_scens(dom_dataset, spatial_data, site_ids, timeframe)::YAXArray{Float64}
    # Add 1 to every scenarios so they represent indexes in cyclone_mr vectors
    cyclone_scens::YAXArray = Cube(
        dom_dataset[["record_applied_cyclone"]]
    )[timestep=At(timeframe[1]:timeframe[2])] .+ 1

    species::Vector{Symbol} = functional_group_names()
    cyclone_mortality_scens::YAXArray{Float64} = ZeroDataCube(;
        T=Float64,
        timesteps=timeframe[1]:timeframe[2],
        locations=site_ids,
        species=species,
        scenarios=1:length(cyclone_scens.scenario)
    )

    # Mortality rate was generated using rrap_dg
    cyclone_mr::NamedTuple = _cyclone_mr()

    # To filter massives/branchings
    massives::BitVector = contains.(String.(species), ["massives"])
    branchings::BitVector = .!massives

    # Set massives mortality rates
    mr_massives::Vector{Float64} = cyclone_mr[:massives]
    cm_scens_massives = mr_massives[cyclone_scens].data
    for m in species[massives]
        cyclone_mortality_scens[species=At(m)] .= cm_scens_massives
    end

    # Set branchings deeper than 5 mortality rates
    mask_d5::BitVector = spatial_data.Y_COORD .<= -5
    if sum(mask_d5) > 0
        mr_bd5::Vector{Float64} = cyclone_mr[:branching_deeper_than_5]
        cm_scens_bd5::Array{Float64} = mr_bd5[cyclone_scens[location=mask_d5]].data
        for b in species[branchings]
            cyclone_mortality_scens[location=(mask_d5), species=At(b)] .= cm_scens_bd5
        end
    end

    # Set branchings shallower than 5 mortality rates
    mask_s5::BitVector = spatial_data.Y_COORD .> -5
    if sum(mask_s5) > 0
        mr_bs5::Vector{Float64} = cyclone_mr[:branching_shallower_than_5]
        cm_scens_bs5::Array{Float64} = mr_bs5[cyclone_scens[location=(mask_s5)]].data
        for b in species[branchings]
            cyclone_mortality_scens[locations=(mask_s5), species=At(b)] .= cm_scens_bs5
        end
    end

    return cyclone_mortality_scens
end

"""
    _cyclone_mr()

For each cyclone category 0 to 5, represented as indexes 1 to 6 of the arrays, there can be
distinct mortality rates for each coral functional group and location.
For Acropora (branching) the mortality depends on the depth (if it is deeper or shalower
than 5 meters);
For massives the mortality does not depend on the depth.
"""
function _cyclone_mr()::NamedTuple
    return (branching_deeper_than_5=[0.0, 0.0, 0.104957, 0.979102, 0.999976, 1.0],
        branching_shallower_than_5=[0.0, 0.0, 0.00993649, 0.497092, 0.989037, 0.99994],
        massives=[0.0, 0.0121482, 0.0155069, 0.0197484, 0.0246357, 0.0302982])
end
