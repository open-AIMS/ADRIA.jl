using ADRIA: SimConstants, Domain, GDF, DataCube

using Distributions, Statistics

using
    DataFrames,
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
    const loc_data
    const loc_id_col
    const cluster_id_col
    init_coral_cover
    const coral_growth::CoralGrowth
    const loc_ids
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
    isdir(fn_path) ? true : error("Path does not exist or is not a directory.")
    netcdf_file = _find_netcdf(fn_path, RCP)
    dom_dataset::Dataset = open_dataset(netcdf_file)

    # Read location data
    geodata_dir = joinpath(fn_path, "region")
    geodata_fn = _find_file(geodata_dir, ".gpkg")
    spatial_data = GDF.read(geodata_fn)

    _standardize_cluster_ids!(spatial_data)

    reef_id_col = "UNIQUE_ID"
    cluster_id_col = "UNIQUE_ID"
    location_ids = spatial_data[:, reef_id_col]

    # Load accompanying ID list
    # TODO: Create canonical geopackage file that aligns all IDs.
    #       Doing this removes the need for the manual correction below and removes the
    #       dependency on this file.
    id_list_fn = _find_file(joinpath(fn_path, "id"), Regex("id_list.*.csv"))
    id_list = CSV.read(
        id_list_fn,
        DataFrame;
        header=false,
        comment="#"
    )

    _manual_id_corrections!(spatial_data, id_list)

    # Load reef area and convert from km^2 to m^2
    spatial_data[:, :area] = id_list[:, 2] .* 1e6

    # Calculate `k` area (1.0 - "ungrazable" area)
    spatial_data[:, :k] = 1 .- id_list[:, 3]

    # Load DHWs
    # Redefine dimensions as ReefMod Matfiles do not contain reef ids.
    # Forcibly load data as disk arrays are not fully support.
    cube_axes = caxes(dom_dataset.record_applied_DHWs)
    dhw_scens = DataCube(
        read(dom_dataset.record_applied_DHWs);
        timestep=Int64.(collect(cube_axes[1])),
        location=location_ids,
        scenario=Int64.(collect(cube_axes[3]))
    )[timestep=At(timeframe[1]:timeframe[2])]

    dhw_scens = correct_axis_names!(dhw_scens)

    # Initial coral cover is loaded from the first year of reefmod 'coral_cover_per_taxa' data
    init_coral_cover = load_initial_cover(
        ReefModDomain, dom_dataset, spatial_data, location_ids, timeframe[1]
    )

    # Connectivity data is retireved from a subdirectory because it's not contained in matfiles
    conn_data = load_connectivity(RMEDomain, fn_path, location_ids)

    spatial_data[:, :depth_med] .= 7.0
    spatial_data[!, :depth_med] = convert.(Float64, spatial_data[!, :depth_med])
    # GBRMPA zone types are not contained in matfiles
    spatial_data[:, :zone_type] .= ["" for _ in 1:nrow(spatial_data)]

    dist_matrix = distance_matrix(spatial_data)
    spatial_data.mean_to_neighbor .= nearest_neighbor_distances(dist_matrix, 10)

    # timesteps, location, scenario
    wave_scens = ZeroDataCube(;
        T=Float64,
        timesteps=timeframe[1]:timeframe[2],
        locs=location_ids,
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
        timeframe[1]:timeframe[2]
    )

    criteria_weights::Vector{Union{DecisionWeights,DecisionThresholds}} = [
        SeedCriteriaWeights(),
        FogCriteriaWeights(),
        DepthThresholds()
    ]

    cyclone_mortality_scens::YAXArray{Float64,4} = _cyclone_mortality_scens(
        dom_dataset,
        spatial_data,
        location_ids,
        timeframe
    )

    model::Model = Model((
        EnvironmentalLayer(dhw_scens, wave_scens, cyclone_mortality_scens),
        Intervention(),
        criteria_weights...,
        Coral(),
        GrowthAcceleration()
    ))

    return ReefModDomain(
        "ReefMod",
        RCP,
        env_md,
        "",
        conn_data,
        spatial_data,
        reef_id_col,
        cluster_id_col,
        init_coral_cover,
        CoralGrowth(nrow(spatial_data)),
        location_ids,
        dhw_scens,
        wave_scens,
        cyclone_mortality_scens,
        model,
        SimConstants()
    )
end

"""
    load_initial_cover(
        ::Type{ReefModDomain},
        dom_data::Dataset,
        location_data::DataFrame,
        loc_ids::Vector{String},
        init_yr::Int=2022
    )::YAXArray
"""
function load_initial_cover(
    ::Type{ReefModDomain},
    dom_data::Dataset,
    location_data::DataFrame,
    loc_ids::Vector{String},
    init_yr::Int=2022
)::YAXArray
    if !haskey(dom_data.cubes, :coral_cover_per_taxa)
        @error "coral_cover_per_taxa variable not found in ReefMod data"
    end

    c_axes = caxes(dom_data.coral_cover_per_taxa)
    init_cc_per_taxa::YAXArray = DataCube(
        read(dom_data.coral_cover_per_taxa);
        timestep=Int64.(collect(c_axes[1])),
        location=Int64.(collect(c_axes[2])),
        group=Int64.(collect(c_axes[3])),
        scenario=Int64.(collect(c_axes[4]))
    )[timestep=At(init_yr)]

    init_cc_per_taxa = init_cc_per_taxa[group=At(2:6)]
    # The following class weight calculations are taken from ReefModEngine Domain calculation

    # Use ReefMod distribution for coral size class population (shape parameters have units log(cm^2))
    # as suggested by YM (pers comm. 2023-08-08 12:55pm AEST). Distribution is used to split ReefMod initial
    # species covers into ADRIA's 6 size classes by weighting with the cdf.
    reef_mod_area_dist = LogNormal(log(700), log(4))
    bin_edges_area = colony_mean_area(bin_edges(; unit=:cm))

    # Find integral density between bounds of each size class areas to create weights for each size class.
    cdf_integral = cdf.(reef_mod_area_dist, bin_edges_area)
    size_class_weights = (cdf_integral[:, 2:end] .- cdf_integral[:, 1:(end - 1)])
    size_class_weights = size_class_weights ./ sum(size_class_weights; dims=2)

    # Take the mean over repeats, as suggested by YM (pers comm. 2023-02-27 12:40pm AEDT).
    # Convert from percent to relative values.
    # YAXArray ordering is [time ⋅ location ⋅ scenario]
    icc_data =
        (
            (dropdims(mean(init_cc_per_taxa; dims=:scenario); dims=:scenario)) ./ 100.0
        ).data
    # Repeat species over each size class and reshape to give ADRIA compatible size
    # [(n_groups × n_sizes) ⋅ n_locs]
    # Multiply by size class weights to give initial cover distribution over each size class.
    n_sizes = size(size_class_weights, 2)
    n_groups = size(size_class_weights, 1)
    n_locs = size(icc_data, 1)
    icc_data = icc_data .* reshape(size_class_weights, (1, n_groups, n_sizes))
    # Reshape and Permute icc_data from
    # [n_locs ⋅ n_groups ⋅ n_sizes] to [(n_groups × n_sizes) ⋅ n_locs]
    icc_data = reshape(
        permutedims(icc_data, (3, 2, 1)), (n_sizes * n_groups, n_locs)
    )

    # Convert values relative to absolute area to values relative to k area
    icc_data = icc_data ./ location_data.k'
    icc_data_sum = dropdims(sum(icc_data; dims=1); dims=1)
    if any(icc_data_sum .> 1.0)
        msg = "Initial cover exceeds habitable area, "
        msg *= "there is most likely an issue with the data package. Constraining to 1.0"
        @warn msg
        icc_data[:, icc_data_sum .> 1.0] ./= icc_data_sum[icc_data_sum .> 1.0]'
    end

    return DataCube(icc_data; species=1:(n_groups * n_sizes), locs=loc_ids)
end

"""
    _find_file(dir::String)::String
"""
function _find_file(dir::String, ident::Union{Regex,String})::String
    pos_files = filter(isfile, readdir(dir; join=true))
    pos_files = filter(x -> occursin(ident, x), pos_files)
    if length(pos_files) == 0
        ArgumentError("Unable to find file in $(dir)")
    elseif length(pos_files) > 1
        @info "Find multiple files matching identifier, using first"
    end
    return pos_files[1]
end

function _find_netcdf(dir::String, scenario::String)::String
    pos_files = filter(isfile, readdir(dir; join=true))
    pos_files = filter(x -> occursin(".nc", x), pos_files)
    pos_files = filter(x -> occursin(scenario, x), pos_files)
    if length(pos_files) == 0
        ArgumentError("Unable to find NetCDF file relating to scenario: $(scenario)")
    elseif length(pos_files) > 1
        @info "Find multiple NetCDF files relating to scenario, using first"
    end
    return pos_files[1]
end

function correct_axis_names!(reefmod_cube::YAXArray)::YAXArray
    reefmod_cube = renameaxis!(reefmod_cube, :timestep => :timesteps)
    reefmod_cube = renameaxis!(reefmod_cube, :location => :locs)
    reefmod_cube = renameaxis!(reefmod_cube, :scenario => :scenarios)
    return reefmod_cube
end

"""
    switch_RCPs!(d::ReefModDomain, RCP::String)::ReefModDomain

Load different RCP into domain.
"""
function switch_RCPs!(d::ReefModDomain, RCP::String)::ReefModDomain
    new_scen_fn = _find_netcdf(d.env_layer_md.dpkg_path, RCP)
    new_scen_dataset = open_dataset(new_scen_fn)

    cube_axes = caxes(new_scen_dataset.record_applied_DHWs)
    dhws = DataCube(
        read(new_scen_dataset.record_applied_DHWs);
        timestep=Int64.(collect(cube_axes[1])),
        location=d.loc_ids,
        scenario=Int64.(collect(cube_axes[3]))
    )[timestep=At(d.env_layer_md.timeframe)]

    dhws = correct_axis_names!(dhws)

    d.dhw_scens = dhws

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
    _cyclone_mortality_scens(
        dom_dataset::Dataset,
        spatial_data::DataFrame,
        location_ids::Vector{String},
        timeframe::Tuple{Int64,Int64}
    )::YAXArray{Float64}

Cyclone scenarios (from 0 to 5) are converted to mortality rates.

# See also
1. https://github.com/open-AIMS/rrap-dg/blob/main/rrap_dg/cyclones/datacube_generator.jl
"""
function _cyclone_mortality_scens(
    dom_dataset::Dataset,
    spatial_data::DataFrame,
    location_ids::Vector{String},
    timeframe::Tuple{Int64,Int64}
)::YAXArray{Float64}
    # Add 1 to every scenarios so they represent indexes in cyclone_mr vectors
    cyclone_data::YAXArray = Cube(dom_dataset[["record_applied_cyclone"]])
    c_axes = caxes(cyclone_data)
    cyclone_scens::YAXArray =
        DataCube(
            read(cyclone_data);
            timestep=Int64.(collect(c_axes[1])),
            location=Int64.(collect(c_axes[2])),
            scenario=Int64.(collect(c_axes[3]))
        )[timestep=At(timeframe[1]:timeframe[2])] .+ 1

    species::Vector{Symbol} = functional_group_names()
    cyclone_mortality_scens::YAXArray{Float64} = ZeroDataCube(;
        T=Float64,
        timesteps=timeframe[1]:timeframe[2],
        locations=location_ids,
        species=species,
        scenarios=1:length(cyclone_scens.scenario)
    )

    # Mortality rate was generated using rrap_dg
    cyclone_mr::NamedTuple = _cyclone_mortalities()

    # To filter massives/branchings
    massives::BitVector = contains.(String.(species), ["massives"])
    branchings::BitVector = .!massives

    # Set massives mortality rates
    mr_massives::Vector{Float64} = cyclone_mr[:massives]
    cm_scens_massives = mr_massives[cyclone_scens].data
    for m in species[massives]
        cyclone_mortality_scens[species=At(m)] .= cm_scens_massives
    end

    # Set mortality rates for branching corals at <= 5m depth
    below_5::BitVector = spatial_data.depth_med .<= -5
    if sum(below_5) > 0
        mr_bd5::Vector{Float64} = cyclone_mr[:branching_deeper_than_5]
        cm_scens_bd5::Array{Float64} = mr_bd5[cyclone_scens[location=below_5]].data
        for b in species[branchings]
            cyclone_mortality_scens[locations=below_5, species=At(b)] .= cm_scens_bd5
        end
    end

    # Set mortality rates for branching corals at > 5m depth
    above_5::BitVector = spatial_data.depth_med .> -5
    if sum(above_5) > 0
        mr_bs5::Vector{Float64} = cyclone_mr[:branching_shallower_than_5]
        cm_scens_bs5::Array{Float64} = mr_bs5[cyclone_scens[location=above_5]].data
        for b in species[branchings]
            cyclone_mortality_scens[locations=above_5, species=At(b)] .= cm_scens_bs5
        end
    end

    return cyclone_mortality_scens
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
