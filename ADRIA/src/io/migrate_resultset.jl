# migrate_resultset.jl
#
# One-off utility to bring a ResultSet saved before the intervention-parameter rename
# (seed -> CAq, mc -> LvM, fog -> Fog, shade/SRM -> Shd, mcb -> MCB, and the iv_<abbrev>_
# prefix pass) up to date, so it can be re-loaded and its `inputs` columns matched by the
# renamed code paths.
#
# What does NOT need migrating (left as-is on disk, confirmed by tracing result_io.jl /
# ResultSet.jl):
#   - `rs.ranks`'s intervention-axis labels are derived live from `interventions()` at
#     load time, not stored on disk. Order was preserved across the rename, so old data
#     already reads correctly with the new labels.
#   - `moving_corals` and `shading_log`/`rankings` on-disk array names are untouched (only
#     `seed` was renamed on disk — see below), so the Julia struct fields that receive
#     them (`rs.LvM_log`, ...) are still just relabelling an unchanged on-disk name.
#
# What DOES need migrating:
#   - `inputs.attrs["columns"]` — the scenario table's column names, stored as a plain
#     list of strings in the `inputs` Zarr array's `.zattrs`. This is the only place the
#     *old* Factor/criteria-weight names are actually persisted.
#   - The `seed` log array's on-disk name. Originally this was left alone deliberately
#     (the Julia field name and the on-disk Zarr name are independent, so renaming
#     `rs.seed_log` -> `rs.CAq_log` needed no data migration). That decision was reversed
#     on request: the on-disk name now also reads `coral_aquaculture` (result_io.jl's
#     `zcreate(...; name="coral_aquaculture", ...)`, `ResultSet.jl`'s
#     `log_set["coral_aquaculture"]`), so old result stores need their `logs/seed`
#     directory physically renamed to `logs/coral_aquaculture` — see
#     `migrate_resultset_logs!` below. `moving_corals` was NOT renamed on disk (not asked
#     for yet) — if it is, add an entry to `_RENAMED_LOG_GROUPS`.

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

const _RENAMED_LOG_GROUPS = Dict{String,String}(
    "seed" => "coral_aquaculture"
    # "moving_corals" => "larval_methods",  # not renamed on disk yet -- add here if it is
)

"""
    migrate_resultset_logs!(result_loc::String; dry_run::Bool=true)::Vector{String}

Rename the on-disk log array directories of a saved `ResultSet` from their pre-rename
names to their current equivalents. See `_RENAMED_LOG_GROUPS` for the full old -> new
mapping (currently just `seed -> coral_aquaculture`).

Each rename is a plain filesystem directory move (`logs/seed` -> `logs/coral_aquaculture`)
— Zarr v2 array metadata does not store the array's own name anywhere inside itself (the
name is purely the directory it lives in), so no chunk data or `.zarray`/`.zattrs` content
needs to change, only the containing folder's name.

# Arguments
- `result_loc` : Path to a saved ResultSet directory (same path you would pass to
  `ADRIA.load_results`).
- `dry_run` : If `true` (default), reports what would move without touching the
  filesystem. Set to `false` to actually rename the directories.

# Returns
The list of old log group names that were found and would be (or were) renamed. Empty if
the result set already uses current naming.
"""
function migrate_resultset_logs!(result_loc::String; dry_run::Bool=true)::Vector{String}
    log_dir = joinpath(result_loc, LOG_GRP)
    isdir(log_dir) ||
        error("No `logs` directory found at $(log_dir) — is this an ADRIA result location?")

    to_rename = String[]
    for (old_name, new_name) in _RENAMED_LOG_GROUPS
        old_path = joinpath(log_dir, old_name)
        isdir(old_path) || continue
        new_path = joinpath(log_dir, new_name)
        isdir(new_path) &&
            error("Both $(old_path) and $(new_path) exist — refusing to overwrite. Resolve manually.")
        push!(to_rename, old_name)
    end

    if isempty(to_rename)
        @info "No log directories need renaming at $(log_dir) — already current, or not an ADRIA log group."
        return String[]
    end

    @info "Found $(length(to_rename)) log director(y/ies) to rename:" to_rename
    if dry_run
        @info "dry_run=true — nothing moved. Re-run with dry_run=false to apply."
        return to_rename
    end

    for old_name in to_rename
        old_path = joinpath(log_dir, old_name)
        new_path = joinpath(log_dir, _RENAMED_LOG_GROUPS[old_name])
        mv(old_path, new_path)
        @info "Renamed $(old_path) -> $(new_path)"
    end

    return to_rename
end

"""
    migrate_resultset!(result_loc::String; dry_run::Bool=true)::Nothing

Bring a ResultSet saved before the intervention-parameter rename fully up to date: renames
the `seed` log directory to `coral_aquaculture` (`migrate_resultset_logs!`) and rewrites
`inputs.attrs["columns"]` to the current `iv_<abbrev>_<name>` names
(`migrate_resultset_columns!`). Convenience wrapper — call the two functions separately if
you want to apply/inspect them independently.

# Example
```julia
ADRIA.migrate_resultset!("path/to/old/result/set")             # dry run, reports only
ADRIA.migrate_resultset!("path/to/old/result/set"; dry_run=false)  # actually migrate
```
"""
function migrate_resultset!(result_loc::String; dry_run::Bool=true)::Nothing
    migrate_resultset_logs!(result_loc; dry_run=dry_run)
    migrate_resultset_columns!(result_loc; dry_run=dry_run)
    return nothing
end
