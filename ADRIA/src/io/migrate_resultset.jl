# migrate_resultset.jl
#
# One-off utility to bring a ResultSet saved before the intervention-parameter rename
# (seed -> CAq, mc -> LvM, fog -> Fog, shade/SRM -> Shd, mcb -> MCB, and the iv_<abbrev>_
# prefix pass) up to date, so it can be re-loaded and its `inputs` columns matched by the
# renamed code paths.
#
# What does NOT need migrating (left as-is on disk, confirmed by tracing result_io.jl /
# ResultSet.jl):
#   - The log arrays themselves (`seed`, `moving_corals`, `shading_log`, `rankings`, ...)
#     are opened by their on-disk literal Zarr names, which were never renamed — only the
#     Julia struct field that receives them (`rs.CAq_log`, `rs.LvM_log`, ...) changed. No
#     data or folder needs to move.
#   - `rs.ranks`'s intervention-axis labels are derived live from `interventions()` at
#     load time, not stored on disk. Order was preserved across the rename, so old data
#     already reads correctly with the new labels.
#
# What DOES need migrating:
#   - `inputs.attrs["columns"]` — the scenario table's column names, stored as a plain
#     list of strings in the `inputs` Zarr array's `.zattrs`. This is the only place the
#     *old* Factor/criteria-weight names are actually persisted.

const _RENAMED_COLUMNS = Dict{String,String}(
    "N_seed_TA" => "iv_CAq_N_TA",
    "N_seed_CA" => "iv_CAq_N_CA",
    "N_seed_CNA" => "iv_CAq_N_CNA",
    "N_seed_SM" => "iv_CAq_N_SM",
    "N_seed_LM" => "iv_CAq_N_LM",
    "N_mc_settlers" => "iv_LvM_N_settlers",
    "seeding_devices_per_m2" => "iv_CAq_devices_per_m2",
    "mc_min_iv_locations" => "iv_LvM_min_iv_locations",
    "fogging" => "iv_Fog",
    "SRM" => "iv_Shd",
    "a_adapt" => "iv_CAq_a_adapt",
    "a_adapt_ref" => "iv_CAq_a_adapt_ref",
    "seed_years" => "iv_CAq_years",
    "shade_years" => "iv_Shd_years",
    "fog_years" => "iv_Fog_years",
    "seed_deployment_freq" => "iv_CAq_deployment_freq",
    "seed_revisit_cadence" => "iv_CAq_revisit_cadence",
    "fog_deployment_freq" => "iv_Fog_deployment_freq",
    "fog_revisit_cadence" => "iv_Fog_revisit_cadence",
    "shade_deployment_freq" => "iv_Shd_deployment_freq",
    "mc_deployment_freq" => "iv_LvM_deployment_freq",
    "mc_revisit_cadence" => "iv_LvM_revisit_cadence",
    "seed_year_start" => "iv_CAq_year_start",
    "shade_year_start" => "iv_Shd_year_start",
    "fog_year_start" => "iv_Fog_year_start",
    "mc_year_start" => "iv_LvM_year_start",
    "mc_years" => "iv_LvM_years",
    "mcb_albedo" => "iv_MCB_albedo",
    "mcb_duration" => "iv_MCB_duration",
    "mcb_deployment_freq" => "iv_MCB_deployment_freq",
    "seed_strategy" => "iv_CAq_strategy",
    "fog_strategy" => "iv_Fog_strategy",
    "mc_strategy" => "iv_LvM_strategy",
    "seed_heat_stress" => "iv_CAq_heat_stress",
    "seed_wave_stress" => "iv_CAq_wave_stress",
    "seed_in_connectivity" => "iv_CAq_in_connectivity",
    "seed_out_connectivity" => "iv_CAq_out_connectivity",
    "seed_depth" => "iv_CAq_depth",
    "seed_coral_cover" => "iv_CAq_coral_cover",
    "seed_cluster_diversity" => "iv_CAq_cluster_diversity",
    "seed_geographic_separation" => "iv_CAq_geographic_separation",
    "fog_heat_stress" => "iv_Fog_heat_stress",
    "fog_wave_stress" => "iv_Fog_wave_stress",
    "fog_in_connectivity" => "iv_Fog_in_connectivity",
    "fog_out_connectivity" => "iv_Fog_out_connectivity",
    "fog_depth" => "iv_Fog_depth",
    "fog_coral_cover" => "iv_Fog_coral_cover",
    "fog_cluster_diversity" => "iv_Fog_cluster_diversity",
    "fog_geographic_separation" => "iv_Fog_geographic_separation",
    "mc_heat_stress" => "iv_LvM_heat_stress",
    "mc_wave_stress" => "iv_LvM_wave_stress",
    "mc_in_connectivity" => "iv_LvM_in_connectivity",
    "mc_out_connectivity" => "iv_LvM_out_connectivity",
    "mc_depth" => "iv_LvM_depth",
    "mc_coral_cover" => "iv_LvM_coral_cover",
    "mc_cluster_diversity" => "iv_LvM_cluster_diversity",
    "mc_geographic_separation" => "iv_LvM_geographic_separation"
    # srm_* criteria weights deliberately excluded: ShdCriteriaWeights was never wired
    # into the module (see DecisionWeights.jl), so no ResultSet's `inputs` can contain
    # those columns.
)

"""
    migrate_resultset_columns!(result_loc::String; dry_run::Bool=true)::Vector{String}

Rewrite the `inputs.attrs["columns"]` list of a saved `ResultSet` (an ADRIA-format,
`DirectoryStore`-backed Zarr result location) from pre-rename intervention parameter
names to their current `iv_<abbrev>_<name>` equivalents. See `_RENAMED_COLUMNS` for the
full old -> new mapping.

Only rewrites the small `.zattrs` metadata for the `inputs` array — no bulk numeric data,
and no other array in the result store, is touched or copied. Safe to run against a
result set produced by any ADRIA version, including ones that already use the current
naming (in which case nothing changes).

# Arguments
- `result_loc` : Path to a saved ResultSet directory (same path you would pass to
  `ADRIA.load_results`).
- `dry_run` : If `true` (default), reports what would change without writing anything.
  Set to `false` to actually rewrite the attrs.

# Returns
The list of old column names that were found and would be (or were) renamed. Empty if
the result set already uses current naming.

# Example
```julia
# See what would change, without writing:
ADRIA.migrate_resultset_columns!("path/to/old/result/set")

# Actually apply the rename:
ADRIA.migrate_resultset_columns!("path/to/old/result/set"; dry_run=false)
```
"""
function migrate_resultset_columns!(result_loc::String; dry_run::Bool=true)::Vector{String}
    input_path = joinpath(result_loc, INPUTS)
    isdir(input_path) ||
        error("No `inputs` array found at $(input_path) — is this an ADRIA result location?")

    # dry_run reads read-only ("r", zopen's default); actually writing needs mode="w" —
    # Zarr.jl only marks a ZArray `writeable` when opened that way (see ZArray.jl:126).
    input_set = zopen(input_path, dry_run ? "r" : "w"; fill_as_missing=false)
    old_columns::Vector{String} = input_set.attrs["columns"]

    to_rename = filter(c -> haskey(_RENAMED_COLUMNS, c), old_columns)

    if isempty(to_rename)
        @info "No columns need renaming at $(result_loc) — already current, or not an intervention-parameter column set."
        return String[]
    end

    @info "Found $(length(to_rename)) column(s) to rename:" to_rename
    if dry_run
        @info "dry_run=true — no changes written. Re-run with dry_run=false to apply."
        return to_rename
    end

    new_columns = [get(_RENAMED_COLUMNS, c, c) for c in old_columns]
    input_set.attrs["columns"] = new_columns
    # Mutating .attrs in place does NOT persist to disk on its own — Zarr.jl only writes
    # attrs at array-creation time (see zcreate -> writeattrs in ZArray.jl). Explicitly
    # flush the updated dict back to inputs/.zattrs.
    Zarr.writeattrs(
        Zarr.ZarrFormat(Zarr.zarr_format(input_set)),
        input_set.storage,
        input_set.path,
        input_set.attrs
    )
    @info "Rewrote inputs.attrs[\"columns\"] at $(result_loc)."

    return to_rename
end
