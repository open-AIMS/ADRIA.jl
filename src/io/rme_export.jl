module RMEExport

using DataFrames, NetCDF, CSV, JSON, YAXArrays, DimensionalData
using ADRIA: ResultSet, timesteps, total_absolute_cover, relative_cover, relative_loc_taxa_cover, relative_shelter_volume, relative_juveniles
using ADRIA: to_coral_spec

export export_to_rme

"""
    export_to_rme(rs::ResultSet, out_dir::String)

Exports an ADRIA `ResultSet` into the format expected by `cost-eco-model-linker`,
emulating the output of ReefModEngine.jl.

# Arguments
- `rs` : ResultSet to export
- `out_dir` : Directory to save the exported files
"""
function export_to_rme(rs::ResultSet, out_dir::String)
    mkpath(out_dir)
    
    # 1. Extract dimensions
    n_timesteps = length(rs.outcomes[:relative_cover].timesteps)
    n_locs = length(rs.outcomes[:relative_cover].locations)
    n_scens = length(rs.outcomes[:relative_cover].scenarios)
    
    # Check if relative_loc_taxa_cover is available (we just added it)
    if !haskey(rs.outcomes, :relative_loc_taxa_cover)
        @warn "Result set is missing `relative_loc_taxa_cover`. Please re-run the scenarios with the updated ADRIA code."
        return
    end
    n_taxa = length(rs.outcomes[:relative_loc_taxa_cover].groups)

    # 2. Build the NetCDF data
    # Pre-allocate RME-formatted YAXArrays
    # RME expects dimensions (timesteps, locations, scenarios) mostly
    # and (timesteps, locations, taxa, scenarios) for taxa cover
    
    # Convert ADRIA's proportions (0-1) to RME's percentages (0-100)
    total_cover_data = Array(rs.outcomes[:relative_cover].data) .* 100.0f0
    
    # Reorder relative_loc_taxa_cover from (timesteps, groups, locations, scenarios)
    # to (timesteps, locations, taxa, scenarios) and multiply by 100
    taxa_data_orig = Array(rs.outcomes[:relative_loc_taxa_cover].data)
    total_taxa_cover_data = permutedims(taxa_data_orig, (1, 3, 2, 4)) .* 100.0f0
    
    # Juveniles: ADRIA relative_juveniles is density (individuals/m2)
    juveniles_data = Array(rs.outcomes[:relative_juveniles].data)
    
    # Relative shelter volume
    shelter_volume_data = Array(rs.outcomes[:relative_shelter_volume].data)
    
    # Placeholders for cots and rubble, required by the linker
    cots_data = zeros(Float32, n_timesteps, n_locs, n_scens)
    rubble_data = zeros(Float32, n_timesteps, n_locs, n_scens)

    # Note: Linker expects these variable names explicitly
    vars = [
        "total_cover" => total_cover_data,
        "coral_juv_m2" => juveniles_data,
        "relative_shelter_volume" => shelter_volume_data,
        "rubble" => rubble_data,
        "cots" => cots_data
    ]

    # Save to NetCDF
    nc_path = joinpath(out_dir, "results.nc")
    # Clean up existing file to prevent NetCDF dimension mismatch errors
    isfile(nc_path) && rm(nc_path)
    
    # Create the variables using NetCDF.jl directly to perfectly emulate RME
    nccreate(
        nc_path, "total_cover",
        "timesteps", n_timesteps,
        "locations", n_locs,
        "scenarios", n_scens,
        t=NC_FLOAT
    )
    nccreate(
        nc_path, "coral_juv_m2",
        "timesteps", n_timesteps,
        "locations", n_locs,
        "scenarios", n_scens,
        t=NC_FLOAT
    )
    nccreate(
        nc_path, "relative_shelter_volume",
        "timesteps", n_timesteps,
        "locations", n_locs,
        "scenarios", n_scens,
        t=NC_FLOAT
    )
    nccreate(
        nc_path, "rubble",
        "timesteps", n_timesteps,
        "locations", n_locs,
        "scenarios", n_scens,
        t=NC_FLOAT
    )
    nccreate(
        nc_path, "cots",
        "timesteps", n_timesteps,
        "locations", n_locs,
        "scenarios", n_scens,
        t=NC_FLOAT
    )
    nccreate(
        nc_path, "total_taxa_cover",
        "timesteps", n_timesteps,
        "locations", n_locs,
        "taxa", n_taxa,
        "scenarios", n_scens,
        t=NC_FLOAT
    )

    ncwrite(total_cover_data, nc_path, "total_cover")
    ncwrite(juveniles_data, nc_path, "coral_juv_m2")
    ncwrite(shelter_volume_data, nc_path, "relative_shelter_volume")
    ncwrite(rubble_data, nc_path, "rubble")
    ncwrite(cots_data, nc_path, "cots")
    ncwrite(total_taxa_cover_data, nc_path, "total_taxa_cover")

    # 3. Build scenario_info.json
    # The linker expects "counterfactual" array and reefsets based on the scenario
    # We will identify unguided/non-intervened scenarios as counterfactuals.
    
    # In ADRIA, a scenario is a counterfactual if no outplanting/enrichment occurs
    is_counterfactual = (rs.inputs.N_seed_TA .+ rs.inputs.N_seed_CA .+ rs.inputs.N_seed_SM) .== 0
    # Also considering fogging and shading
    is_counterfactual = is_counterfactual .& (rs.inputs.fogging .== 0) .& (rs.inputs.SRM .== 0)
    
    # Generate a single dummy reefset or extract real reefsets
    # The linker looks at all keys not in ("counterfactual", "dhw_tolerance") as reefsets
    # We'll just define "reefset_all" containing all locations that got intervened anywhere
    intervened_locs = findall(sum(rs.seed_log, dims=(1,2,4))[1,1,:,1] .> 0)
    reefset_ids = rs.loc_ids[intervened_locs]
    if isempty(reefset_ids)
        # fallback if no intervention
        reefset_ids = [rs.loc_ids[1]]
    end

    scenario_info = Dict(
        "counterfactual" => Int.(is_counterfactual),
        "dhw_tolerance" => zeros(Float64, n_scens), # placeholder
        "reefset_adria" => reefset_ids
    )

    open(joinpath(out_dir, "scenario_info.json"), "w") do f
        write(f, JSON.json(scenario_info))
    end

    # 4. Build iv_yearly_scenarios.csv
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
    
    start_year = parse(Int, string(rs.outcomes[:relative_cover].timesteps[1]))
    
    # ADRIA's seed_log is in m² of coral coverage, or number of corals?
    # Actually, seed_log stores "Estimated number of corals seeded" (per scenario.jl)
    # Let's sum across species and locations
    for scen in 1:n_scens
        if is_counterfactual[scen]
            continue
        end
        for t in 1:n_timesteps
            # total corals seeded this year
            total_corals = sum(rs.seed_log[t, :, :, scen])
            if total_corals > 0
                year = start_year + t - 1
                
                # Approximate area (km²) and density
                # ADRIA usually seeds into a small portion, we can rough this out
                # or derive it exactly if needed. 
                # For Linker, "number of corals" is the primary driver of cost.
                push!(iv_df, (
                    1, "ADRIA_GCM", "outplant", "reefset_adria", year, scen,
                    total_corals, 
                    1.0, # dummy density
                    (total_corals * 0.01) / 1e6 # dummy area
                ))
            end
        end
    end
    
    CSV.write(joinpath(out_dir, "iv_yearly_scenarios.csv"), iv_df)

    # 5. Build reef_information.csv
    # Expected by Linker? Actually Linker loads reefmod_gbr.gpkg, but RME saves reef_information.csv
    reef_info = DataFrame(
        "reef_id" => rs.loc_ids,
        "area_km2" => rs.loc_area ./ 1e6
    )
    CSV.write(joinpath(out_dir, "reef_information.csv"), reef_info)

    @info "Successfully exported ADRIA ResultSet to RME format at $out_dir"
end

end # module
