module RMEExport

using DataFrames, NetCDF, CSV, JSON, YAXArrays, DimensionalData
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
    export_to_rme(rs::ResultSet, out_dir::String)

Exports an ADRIA `ResultSet` into the format expected by `cost-eco-model-linker`,
emulating the output of ReefModEngine.jl.

# Arguments
- `rs` : ResultSet to export
- `out_dir` : Directory to save the exported files
"""
function export_to_rme(dom::Domain, rs::ResultSet, out_dir::String)
    mkpath(out_dir)

    # 1. Extract dimensions
    n_timesteps = length(rs.outcomes[:relative_cover].timesteps)
    n_locs = length(rs.outcomes[:relative_cover].locations)
    n_scens = length(rs.outcomes[:relative_cover].scenarios)

    total_cover_data = Array(ltmp_cover(rs).data) .* 100.0

    # Juveniles: ADRIA relative_juveniles is density (individuals/m2)
    juveniles_data = Array(rs.outcomes[:relative_juveniles].data)

    # Relative shelter volume
    msg = "Exporting relative shelter volume as half of the ADRIAIndicators metric."
    msg *= "ADRIAIndicators uses a maximum shelter volume 50% lower than used in RME."
    msg *= "Halving the metric aligns them."
    @info msg
    shelter_volume_data = Array(rs.outcomes[:relative_shelter_volume].data) ./ 2.0

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
    seed_factors = [
        :N_seed_TA, :N_seed_CA, :N_seed_CNA, :N_seed_SM, :N_seed_LM, :N_mc_settlers
    ]

    # Filter to only those present in inputs
    seed_cols = intersect(seed_factors, propertynames(rs.inputs))
    intensity = zeros(n_scens)
    for c in seed_cols
        intensity .+= rs.inputs[!, c]
    end

    is_counterfactual = (intensity .== 0)

    # Also check fogging and shading if they exist
    if hasproperty(rs.inputs, :fogging)
        is_counterfactual = is_counterfactual .& (rs.inputs.fogging .== 0)
    end
    if hasproperty(rs.inputs, :SRM)
        is_counterfactual = is_counterfactual .& (rs.inputs.SRM .== 0)
    end

    # Also consider guided status (<= 0 includes unguided and counterfactual)
    if hasproperty(rs.inputs, :guided)
        is_counterfactual = is_counterfactual .& (rs.inputs.guided .<= 0)
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
    # Map: Set(location_indices) -> reefset_name
    reefset_registry = Dict{Set{Int},String}()
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

            total_corals = sum(loc_seeding)
            gcm_name::String = extract_GCM_from_results(dom, rs, scen)

            if !isempty(intervened_loc_indices)
                loc_set = Set(intervened_loc_indices)

                # Assign or retrieve reefset name
                if !haskey(reefset_registry, loc_set)
                    rs_name = "reefset_$(reefset_counter)"
                    reefset_registry[loc_set] = rs_name
                    reefset_counter += 1
                end
                rs_name = reefset_registry[loc_set]

                year = start_year + t - 1

                # Get the actual seeding density for this scenario
                density = rs.inputs[scen, :seeding_devices_per_m2]
                # Calculate the realized area in km2 (m2 / 1e6)
                area_km2 = (total_corals / density) / 1e6

                push!(
                    iv_df,
                    (
                        scen, gcm_name, "outplant", rs_name, year, 1,
                        Float64(total_corals),
                        Float64(density),
                        Float64(area_km2)
                    )
                )
            end

            # Process mc_log: (timesteps, coral_id, locations, scenarios)
            loc_mc = sum(rs.mc_log[t, :, :, scen]; dims=1)
            mc_loc_indices = getindex.(findall(loc_mc .> 0), 2)
            total_mc_corals = sum(loc_mc)

            if !isempty(mc_loc_indices)
                loc_set_mc = Set(mc_loc_indices)
                if !haskey(reefset_registry, loc_set_mc)
                    rs_name_mc = "reefset_$(reefset_counter)"
                    reefset_registry[loc_set_mc] = rs_name_mc
                    reefset_counter += 1
                end
                rs_name_mc = reefset_registry[loc_set_mc]

                year = start_year + t - 1
                push!(
                    iv_df,
                    (
                        scen, gcm_name, "lm", rs_name_mc, year, 1,
                        Float64(total_mc_corals),
                        1.0,  # dummy density
                        0.01  # dummy area
                    )
                )
            end
        end
    end

    # 4. Build scenario_info.json
    # Map reefset names back to location IDs
    scenario_info = Dict{String,Any}(
        "counterfactual" => Int.(is_counterfactual),
        "dhw_tolerance" => rs.inputs.a_adapt
    )

    for (loc_set, rs_name) in reefset_registry
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
