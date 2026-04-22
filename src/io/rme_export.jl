module RMEExport

using DataFrames, NetCDF, CSV, JSON, YAXArrays, DimensionalData
import GeoDataFrames as GDF
import ArchGDAL as AG
import GeoFormatTypes as GFT

using ADRIA: ResultSet, timesteps
using ADRIA.metrics: total_absolute_cover, relative_cover, relative_loc_taxa_cover,
    ltmp_cover, relative_shelter_volume, relative_juveniles
using ADRIA: to_coral_spec, Domain, switch_RCPs!

export export_to_rme

"""
    extract_GCM_from_RME(filepath::String)::String

Extract GCM name from the scenarios dimension in the DHW dimensions.
"""
function extract_GCM_from_RME(filepath::String)::String
    filename = basename(filepath)
    gcm = split(filename, "_SSP")[1]

    return gcm
end

function extract_GCM_from_results(
    dom::Domain, rs::ResultSet, scen_idx::Int64
)::String
    RCP::String = split(string(rs.inputs.RCP[scen_idx]), ".")[1]
    if dom.RCP != RCP
        switch_RCPs!(dom, RCP)
    end

    dhw_scen_idx::Int64 = Int64(rs.inputs.dhw_scenario[scen_idx])

    if eltype(dom.dhw_scens.scenarios) == String
        scen_name = dom.dhw_scens.scenarios[dhw_scen_idx]
        # RME Domain lists dimensions as filepaths to netcdf files.
        if ispath(scen_name)
            return extract_GCM_from_RME(scen_name)
        else
            return scen_name
        end
    end

    return "ADRIA_GCM"
end

"""
    export_to_rme(dom::Domain, rs::ResultSet, out_dir::String; map_counterfactuals::Bool=false, pdp_geopackage_path::Union{String, Nothing}=nothing, scen_spec::Union{DataFrame, Nothing}=nothing)

Exports an ADRIA `ResultSet` into the format expected by `cost-eco-model-linker`,
emulating the output of ReefModEngine.jl.

# Arguments
- `dom` : Domain used to run scenarios
- `rs` : ResultSet to export
- `out_dir` : Directory to save the exported files
- `map_counterfactuals` : Whether to map guided scenarios to their counterfactuals. Defaults to false.
- `pdp_geopackage_path` : Path to GeoPackage containing PDP_region column. Defaults to nothing.
- `scen_spec` : Scenario specification to use. If nothing, uses `rs.inputs`. Defaults to nothing.
"""
function export_to_rme(
    dom::Domain,
    rs::ResultSet,
    out_dir::String;
    map_counterfactuals::Bool=false,
    pdp_geopackage_path::Union{String,Nothing}=nothing,
    scen_spec::Union{DataFrame,Nothing}=nothing
)
    mkpath(out_dir)

    # 1. Extract dimensions
    n_timesteps = length(rs.outcomes[:relative_cover].timesteps)
    n_locs = length(rs.outcomes[:relative_cover].locations)
    n_scens = length(rs.outcomes[:relative_cover].scenarios)

    # Load PDP region information if geopackage is provided
    loc_idx_to_region = if isnothing(pdp_geopackage_path)
        fill("Domain", n_locs)
    else
        pdp_df = GDF.read(pdp_geopackage_path)
        id_to_region = Dict(
            string(row.UNIQUE_ID) => string(row.PDP_region) for row in eachrow(pdp_df)
        )
        [get(id_to_region, id, "outside_pdp") for id in rs.loc_ids]
    end

    total_cover_data = Array(ltmp_cover(rs).data) .* 100.0

    # Juveniles: ADRIA relative_juveniles is density (individuals/m2)
    juveniles_data = Array(rs.outcomes[:relative_juveniles].data)

    # Relative shelter volume
    shelter_volume_data = Array(
        rs.outcomes[:relative_shelter_volume].data
    )

    # Placeholders for cots and rubble
    cots_data = zeros(Float32, n_timesteps, n_locs, n_scens)
    rubble_data = zeros(Float32, n_timesteps, n_locs, n_scens)

    # Save to NetCDF
    nc_path = joinpath(out_dir, "results.nc")
    isfile(nc_path) && rm(nc_path)

    # 1. Define dimensions
    # Note: we omit the 'coral_id' coordinate variable as requested.
    years = Int32.(parse.(Int, string.(lookup(rs.outcomes[:relative_cover], :timesteps))))
    loc_ids = String.(rs.loc_ids)
    scen_indices = Int32.(1:n_scens)

    nccreate(
        nc_path,
        "total_cover",
        "timesteps",
        n_timesteps,
        "locations",
        n_locs,
        "scenarios",
        n_scens;
        t=NC_FLOAT
    )
    nccreate(
        nc_path,
        "relative_juveniles",
        "timesteps",
        n_timesteps,
        "locations",
        n_locs,
        "scenarios",
        n_scens;
        t=NC_FLOAT
    )
    nccreate(
        nc_path,
        "relative",
        "timesteps",
        n_timesteps,
        "locations",
        n_locs,
        "scenarios",
        n_scens;
        t=NC_FLOAT
    )
    nccreate(
        nc_path,
        "relative_shelter_volume",
        "timesteps",
        n_timesteps,
        "locations",
        n_locs,
        "scenarios",
        n_scens;
        t=NC_FLOAT
    )
    nccreate(
        nc_path,
        "rubble",
        "timesteps",
        n_timesteps,
        "locations",
        n_locs,
        "scenarios",
        n_scens;
        t=NC_FLOAT
    )
    nccreate(
        nc_path,
        "cots",
        "timesteps",
        n_timesteps,
        "locations",
        n_locs,
        "scenarios",
        n_scens;
        t=NC_FLOAT
    )

    # 2. Define coordinate variables
    nccreate(nc_path, "timesteps", "timesteps", n_timesteps; t=NC_INT)
    nccreate(nc_path, "locations", "locations", n_locs; t=NC_STRING)
    nccreate(nc_path, "scenarios", "scenarios", n_scens; t=NC_INT)

    # 3. Write data
    ncwrite(total_cover_data, nc_path, "total_cover")
    ncwrite(juveniles_data, nc_path, "relative_juveniles")
    ncwrite(shelter_volume_data, nc_path, "relative")
    ncwrite(shelter_volume_data, nc_path, "relative_shelter_volume")
    ncwrite(rubble_data, nc_path, "rubble")
    ncwrite(cots_data, nc_path, "cots")

    # 4. Write coordinates
    ncwrite(years, nc_path, "timesteps")
    ncwrite(loc_ids, nc_path, "locations")
    ncwrite(scen_indices, nc_path, "scenarios")
    # 3. Dynamic Reef set and IV Scenario synthesis
    # A scenario is a counterfactual if no interventions occur
    # Check all seeding and intervention factors

    # Use provided scenario specification if available, otherwise fallback to rs.inputs
    inputs = isnothing(scen_spec) ? rs.inputs : scen_spec

    seed_factors = [
        :N_seed_TA, :N_seed_CA, :N_seed_CNA, :N_seed_SM, :N_seed_LM, :N_mc_settlers
    ]

    # Filter to only those present in inputs
    seed_cols = intersect(seed_factors, propertynames(inputs))
    intensity = zeros(n_scens)
    for c in seed_cols
        intensity .+= inputs[!, c]
    end

    is_counterfactual = (intensity .== 0)

    # Also check fogging and shading if they exist
    if hasproperty(inputs, :fogging)
        is_counterfactual = is_counterfactual .& (inputs.fogging .== 0)
    end
    if hasproperty(inputs, :SRM)
        is_counterfactual = is_counterfactual .& (inputs.SRM .== 0)
    end

    # Also consider guided status (<= 0 includes unguided and counterfactual)
    if hasproperty(inputs, :guided)
        is_counterfactual = is_counterfactual .& (inputs.guided .<= 0)
    end

    # Group scenarios by intervention settings to assign IDs and repetitions
    env_cols = intersect(
        [:dhw_scenario, :wave_scenario, :cyclone_scenario, :RCP], propertynames(inputs)
    )
    iv_cols = setdiff(propertynames(inputs), env_cols)
    iv_groups = groupby(inputs, iv_cols)
    iv_id_map = zeros(Int, n_scens)
    rep_map = zeros(Int, n_scens)
    for (g_idx, group) in enumerate(iv_groups)
        orig_indices = parentindices(group)[1]
        for (r_idx, scen_idx) in enumerate(orig_indices)
            iv_id_map[scen_idx] = g_idx
            rep_map[scen_idx] = r_idx
        end
    end

    start_year = parse(Int, string(rs.outcomes[:relative_cover].timesteps[1]))

    # Linker expects: intervention id, GCM name, type, reefset, year, rep, number of corals, corals per m2, intervention area km2
    iv_df = DataFrame(
        "intervention id" => Int[],
        "GCM name" => String[],
        "type" => String[],
        "reefset" => String[],
        "year" => Int[],
        "rep" => Int[],
        "number of corals" => Float64[],
        "corals per m2" => Float64[],
        "intervention area km2" => Float64[]
    )

    # Registry of unique location sets to avoid redundancy in scenario_info.json
    # Map: (Region, Set(location_indices)) -> reefset_name
    reefset_registry = Dict{Tuple{String,Set{Int}},String}()
    reefset_counter = 1

    # Process seed_log: (timesteps, coral_id, locations, scenarios)
    # Sum over coral_id (2nd dim) to find intervened locations
    for scen in 1:n_scens
        if is_counterfactual[scen]
            continue
        end
        for t in 1:n_timesteps
            # Find indices of locations where seeding occurred in this (t, scen)
            # Sum over coral_id dimension
            loc_seeding = sum(rs.seed_log[t, :, :, scen]; dims=1)
            intervened_loc_indices = getindex.(findall(loc_seeding .> 0), 2)

            gcm_name::String = extract_GCM_from_results(dom, rs, scen)

            if !isempty(intervened_loc_indices)
                # Group by region
                regional_groups = Dict{String,Vector{Int}}()
                for idx in intervened_loc_indices
                    region = loc_idx_to_region[idx]
                    push!(get!(regional_groups, region, Int[]), idx)
                end

                for (region, indices) in regional_groups
                    loc_set = Set(indices)
                    key = (region, loc_set)

                    # Assign or retrieve reefset name
                    if !haskey(reefset_registry, key)
                        rs_name = "reefset_$(region)_$(reefset_counter)"
                        reefset_registry[key] = rs_name
                        reefset_counter += 1
                    end
                    rs_name = reefset_registry[key]

                    year = start_year + t - 1

                    # Get total corals for this regional group
                    regional_total_corals = sum(loc_seeding[indices])

                    # Get the actual seeding density for this scenario
                    density = inputs[scen, :seeding_devices_per_m2]
                    # Calculate the realized area in km2 (m2 / 1e6)
                    area_km2 = (regional_total_corals / density) / 1e6

                    push!(
                        iv_df,
                        (
                            iv_id_map[scen], gcm_name, "outplant", rs_name, year, rep_map[scen],
                            Float64(regional_total_corals),
                            Float64(density),
                            Float64(area_km2)
                        )
                    )
                end
            end

            # Process mc_log: (timesteps, coral_id, locations, scenarios)
            loc_mc = sum(rs.mc_log[t, :, :, scen]; dims=1)
            mc_loc_indices = getindex.(findall(loc_mc .> 0), 2)

            if !isempty(mc_loc_indices)
                # Group by region
                regional_groups_mc = Dict{String,Vector{Int}}()
                for idx in mc_loc_indices
                    region = loc_idx_to_region[idx]
                    push!(get!(regional_groups_mc, region, Int[]), idx)
                end

                for (region, indices) in regional_groups_mc
                    loc_set_mc = Set(indices)
                    key = (region, loc_set_mc)

                    if !haskey(reefset_registry, key)
                        rs_name_mc = "reefset_$(region)_$(reefset_counter)"
                        reefset_registry[key] = rs_name_mc
                        reefset_counter += 1
                    end
                    rs_name_mc = reefset_registry[key]

                    year = start_year + t - 1

                    # Get total corals for this regional group
                    regional_total_mc_corals = sum(loc_mc[indices])

                    # Calculate actual density and area for MC
                    # loc_area is in m^2
                    total_loc_area_m2 = sum(rs.loc_area[indices])
                    mc_density = regional_total_mc_corals / total_loc_area_m2
                    mc_area_km2 = total_loc_area_m2 / 1e6

                    push!(
                        iv_df,
                        (
                            iv_id_map[scen], gcm_name, "lm", rs_name_mc, year, rep_map[scen],
                            Float64(regional_total_mc_corals),
                            Float64(mc_density),
                            Float64(mc_area_km2)
                        )
                    )
                end
            end
        end
    end

    # 4. Build scenario_info.json
    # Map reefset names back to location IDs
    scenario_info = Dict{String,Any}(
        "counterfactual" => Int.(is_counterfactual),
        "dhw_tolerance" => inputs.a_adapt
    )

    if map_counterfactuals
        env_cols = intersect(
            [:dhw_scenario, :wave_scenario, :cyclone_scenario, :RCP],
            propertynames(inputs)
        )
        cf_indices = findall(is_counterfactual)

        if !isempty(cf_indices)
            cf_profiles = inputs[cf_indices, env_cols]
            mapping = zeros(Int, n_scens)
            for i in 1:n_scens
                current_profile = Vector(inputs[i, env_cols])
                # Find matching CF index
                match_row_idx = findfirst(
                    r -> Vector(r) == current_profile, eachrow(cf_profiles)
                )
                if !isnothing(match_row_idx)
                    mapping[i] = cf_indices[match_row_idx]
                end
            end
            scenario_info["counterfactual_mapping"] = mapping
        end
    end

    for ((region, loc_set), rs_name) in reefset_registry
        # loc_set contains indices, convert to IDs
        scenario_info[rs_name] = rs.loc_data.GBRMPA_ID[collect(loc_set)]
    end

    # Fallback if no interventions occurred
    if isempty(reefset_registry)
        scenario_info["reefset_empty"] = [rs.loc_ids[1]]
    end

    open(joinpath(out_dir, "scenario_info.json"), "w") do f
        write(f, JSON.json(scenario_info))
    end

    CSV.write(joinpath(out_dir, "iv_yearly_scenarios.csv"), iv_df)

    # 5. Build reef_information.csv
    reef_info = DataFrame(
        "reef_id" => rs.loc_ids,
        "area_km2" => rs.loc_area ./ 1e6
    )
    CSV.write(joinpath(out_dir, "reef_information.csv"), reef_info)

    @info "Successfully exported ADRIA ResultSet to RME format with dynamic yearly reefsets at $out_dir"
end

end # module
