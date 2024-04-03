using ADRIA: GDF, ResultSet

using YAXArrays
using NetCDF


struct ReefModResultSet{T1, T2, D, G, D1} <: ResultSet
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

    sim_constants::D1

    # raw::AbstractArray
    outcomes::Dict{Symbol, YAXArray}
    # seed_log::B  # Values stored in m^2
    # fog_log::C   # Reduction in bleaching mortality (0.0 - 1.0)
    # shade_log::C # Reduction in bleaching mortality (0.0 - 1.0)
    # coral_dhw_tol_log::D3
end

"""
    load_results(::Type{ReefModResultSet}, result_loc::String, RCP::String)::ReefModResultSet

Reefmod result interface.
"""
function load_results(::Type{ReefModResultSet}, result_loc::String, RCP::String)::ReefModResultSet
    !isdir(result_loc) ? error("Expected a directory but received $(result_loc)") : nothing

    result_path = _get_netcdf(result_loc, RCP)
    raw_set = open_dataset(result_path)
    
    name::String = raw_set.properties["reef model"]
    netcdf_rcp::String = raw_set.properties["climate scenario (RCP)"]
    if netcdf_rcp != RCP
        @warn "RCPs given and file found do not match. Using RCP of file found."
    end
    
    geodata_path = joinpath(result_loc, "region", "reefmod_gbr.gpkg")
    geodata = GDF.read(geodata_path)
    
    # Correct site ids in the same manner 
    _standardize_cluster_ids!(geodata)

    reef_id_col = "UNIQUE_ID"
    cluster_id_col = "UNIQUE_ID"
    site_ids = geodata[:, reef_id_col]

    # Load accompanying ID list
    id_list_fn = _find_file(joinpath(result_loc, "id"), Regex("id_list.*.csv"))
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
    connectivity = load_connectivity(RMEDomain, result_loc, site_ids)

    location_max_coral_cover = 1 .- geodata.k ./ 100
    location_centroids = [centroid(multipoly) for multipoly ∈ geodata.geom]

    timeframe = 2008:2101
    if !haskey(raw_set.axes, :timestep)
        @warn "Unable to find timestep axes. Defaulting to $(timeframe[1]):$(timeframe[end])"
    else
        timeframe = raw_set.timestep[1]:raw_set.timestep[end]
    end

    scenario_groups = Dict(:counterfactual => BitVector(true for _ in raw_set.scenario))

    env_layer_md::EnvLayer = EnvLayer(
        result_loc,
        geodata_path,
        reef_id_col,
        cluster_id_col,
        "",
        joinpath(result_loc, "con_bin"),
        "",
        "",
        timeframe
    )
    
    outcome_keys = [
        :coral_larval_supply,
        :coral_cover_per_taxa,
        :nb_coral_offspring,
        :nb_coral_adol,
        :COTS_larval_supply,
        :nb_coral_adult,
        :nb_coral_juv,
        :nb_coral_recruit,
        :reef_shelter_volume_relative
    ]
    
    outcomes::Dict{Symbol, YAXArray} = Dict()
    for outcome in outcome_keys
        if !haskey(raw_set.cubes, outcome)
            @warn "Unable to find $(string(outcome)). Skipping."
            continue
        end
        if outcome == :reef_shelter_volume_relative
            outcomes[:relative_shelter_volume] = reformat_cube(ReefModResultSet, raw_set.cubes[outcome])
            continue
        end
        outcomes[outcome] = reformat_cube(ReefModResultSet, raw_set.cubes[outcome])
    end
    outcomes[:relative_cover] = _reefmod_relative_cover(raw_set)

    return ReefModResultSet(
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
        SimConstants(),
        outcomes
    )
end

function _reefmod_relative_cover(dataset::Dataset)::YAXArray
    relative_cover = reformat_cube(ReefModResultSet, dataset.coral_cover_per_taxa)
    return dropdims(sum(relative_cover, dims=:taxa), dims=:taxa) ./ 100
end

"""
    _get_netcdf(result_loc::String, RCP::String)::String
"""
function _get_netcdf(result_loc::String, RCP::String)::String
    possible_files = filter(x->occursin(RCP, x), readdir(result_loc))
    possible_files = filter(x->occursin(".nc", x), possible_files)
    if length(possible_files) == 0
        throw(ArgumentError("Unable to find ReefMod NetCDF file"))
    elseif length(possible_files) > 1
        @warn "Found multiple matching netcdf files. Using first."
    end
    return possible_files[1]
end 

"""
   reformat_cube(cube::YAXArray)::YAXArray

Temporary function. Reefmod NetCDF files should follow expected dimension naming and ordering.
"""
function reformat_cube(::Type{ReefModResultSet}, cube::YAXArray)::YAXArray
    dim_names = name.(cube.axes)
    reefmod_names = [:timestep, :location, :group, :scenario]
    adria_names = [:timesteps, :sites, :taxa, :scenarios]
    final_ordering::Vector{Int} = Vector{Int}(undef, length(dim_names))
    current_index = 1
    # rename expected dimensions and update the permutation vector
    for (cscape_nm, adria_nm) in zip(reefmod_names, adria_names)
        if cscape_nm ∉ dim_names
            continue
        end
        cube = renameaxis!(cube, cscape_nm=>adria_nm)
        final_ordering[current_index] = findfirst(x->x==cscape_nm, dim_names)
        current_index += 1
    end
     # append the indices of non-adria axis to end of permutation vector
    for dim_name in dim_names
        if dim_name in reefmod_names
            continue
        end
        final_ordering[current_index] = findfirst(x->x==dim_name, dim_names)
        current_index += 1;
    end
    if haskey(cube.properties, "units")
        if cube.properties["units"] == "percent"
            cube = cube ./ 100
            cube.properties["units"] == "proportion [0, 1]"
        end
    end
    return permutedims(cube, final_ordering)
end

function Base.show(io::IO, mime::MIME"text/plain", rs::ReefModResultSet)
    rcps = rs.RCP
    scens = length(rs.outcomes[:coral_cover_per_taxa].scenarios)
    sites = length(rs.outcomes[:coral_cover_per_taxa].sites)
    tf = rs.env_layer_md.timeframe

    println("""
    Name: $(rs.name)

    Results stored at: $(rs.env_layer_md.dpkg_path)

    RCP(s) represented: $(rcps)
    Scenarios run: $(scens)
    Number of sites: $(sites)
    Timesteps: $(tf)
    """)
end
