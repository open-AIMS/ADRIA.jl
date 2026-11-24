using Dates: now

using Zarr

using Setfield: @set!
using DataFrames: DataFrame, colmetadata!, nrow
using ADRIA: Domain, EnvLayer

const CONNECTIVITY = "connectivity"
const ENV_STATS = "env_stats"
const INPUTS = "inputs"
const LOG_GRP = "logs"
const MODEL_SPEC = "model_spec"
const RESULTS = "results"
const SPATIAL_DATA = "spatial"

# Combined store name for the 6 location-based outcome metrics
const LOC_METRICS = "loc_outcomes"
const LOC_METRIC_NAMES = [
    :relative_cover, :relative_shelter_volume, :absolute_shelter_volume,
    :relative_juveniles, :juvenile_indicator, :coral_evenness
]

abstract type ResultSet end

struct ADRIAResultSet{T1,T2,A,B,C,D,G,D1,D2,D3,D4,LG,DF} <: ResultSet
    name::String
    RCP::String
    invoke_time::String
    ADRIA_VERSION::String

    loc_ids::T1
    loc_area::Vector{Float64}
    loc_max_coral_cover::Vector{Float64}
    loc_centroids::T2
    env_layer_md::EnvLayer
    dhw_stats::D
    wave_stats::D
    connectivity_data::D
    loc_data::G

    inputs::G
    sim_constants::D1
    model_spec::DF
    model_specs::Vector{DF}  # model spec of each source result set (one if not merged)
    source_sizes::Vector{Int}  # number of scenarios taken from each source, in order
    provenance::Vector{Any}  # provenance record of each source (`nothing` if not stored)

    # raw::AbstractArray
    outcomes::D2
    ranks::A
    iv_LvM_log::B  # Number of individuals deployed per functional group per location
    iv_CAq_log::B  # Number of individuals deployed per functional group per location
    iv_Shd_log::C  # Fog and shade intervention log, dims: (timesteps, locations, intervention, scenarios)
    coral_dhw_tol_log::D3
    coral_cover_log::D4
    logs::LG  # submodel logs by submodel, e.g. `logs.cots.pop`; empty if none
end

"""
Names of the arrays in the log group that are part of the core result set. Submodel logs
are in subgroups of the log group, so these names are reserved as submodel names. `fog` and
`shade` are the old (pre-`shading_log`) separate shading arrays.
"""
const CORE_LOG_NAMES = (
    "rankings", "larval_methods", "coral_aquaculture", "shading_log", "fog", "shade",
    "coral_dhw_log", "coral_cover_log"
)

"""
Result sets stored before logs were grouped by submodel hold the COTS logs as flat arrays
of the log group, under their `LogSpec` names. Maps those array names to
`(submodel, log name)`, the explicit mapping being `COTS_LOG_NAMES`. Any other flat array
that is not a core log is not recognised (such stores never held one).
"""
const LEGACY_LOG_NAMES = NamedTuple{keys(COTS_LOG_NAMES)}(
    map(v -> (COTS_LOG_GROUP, v), values(COTS_LOG_NAMES))
)

"""
Submodel logs that older code exposed as fields of the result set (`rs.cots_pop_log`),
as `(submodel, log name)`; deprecated in favour of `rs.logs`. `nothing` when absent.
"""
const DEPRECATED_LOG_ACCESSORS = LEGACY_LOG_NAMES

"""
    Base.getproperty(rs::ADRIAResultSet, name::Symbol)

Fields, and the deprecated accessors of COTS logs (`rs.cots_pop_log`, `nothing` when
absent as the former fields were). Logs are otherwise reached through `rs.logs`.
"""
@inline function Base.getproperty(rs::ADRIAResultSet, name::Symbol)
    hasfield(ADRIAResultSet, name) && return getfield(rs, name)
    if haskey(DEPRECATED_LOG_ACCESSORS, name)
        group, log = DEPRECATED_LOG_ACCESSORS[name]
        logs = getfield(rs, :logs)
        return haskey(logs, group) && haskey(logs[group], log) ? logs[group][log] : nothing
    end

    throw(
        ArgumentError(
            "ADRIAResultSet has no property `$name` (submodel logs are in `rs.logs`)"
        )
    )
end

"""
    Base.propertynames(rs::ADRIAResultSet, private::Bool=false)

Only the fields (the deprecated accessors are omitted), so that Setfield's `@set!` works.
"""
Base.propertynames(rs::ADRIAResultSet, private::Bool=false) = fieldnames(ADRIAResultSet)

"""
    _load_submodel_logs(log_set::Zarr.ZGroup)::NamedTuple

Load the submodel logs of a log group, as `(submodel = (log = array, ...), ...)` with axes
named by each array's `structure` attribute, sorted by name; empty when there are none.
Every subgroup of the log group is a submodel and its arrays are its logs. Stores written
before logs were grouped hold the COTS logs as flat arrays: those named in
`LEGACY_LOG_NAMES` load as `cots.pop` and `cots.condition`. A store holding a log both flat
and in its subgroup is rejected with an `ArgumentError`.
"""
function _load_submodel_logs(log_set::Zarr.ZGroup)::NamedTuple
    found = Dict{Symbol,Dict{Symbol,Any}}()
    stored_as = Dict{Tuple{Symbol,Symbol},String}()  # path in the log group, for errors
    for g in keys(log_set.groups)
        grp = log_set.groups[g]
        found[Symbol(g)] = Dict{Symbol,Any}(Symbol(n) => grp[n] for n in keys(grp.arrays))
        for n in keys(grp.arrays)
            stored_as[(Symbol(g), Symbol(n))] = "$g/$n"
        end
    end

    for n in keys(log_set.arrays)
        haskey(LEGACY_LOG_NAMES, Symbol(n)) || continue
        group, log = LEGACY_LOG_NAMES[Symbol(n)]
        haskey(stored_as, (group, log)) && throw(
            ArgumentError(
                "Log `$group.$log` is stored twice in the log group, as " *
                "`$(stored_as[(group, log)])` and as the flat array `$n`"
            )
        )

        get!(found, group, Dict{Symbol,Any}())[log] = log_set[n]
        stored_as[(group, log)] = n
    end

    groups = sort!(collect(keys(found)))
    loaded = map(groups) do g
        logs = sort!(collect(keys(found[g])))
        for l in logs
            haskey(found[g][l].attrs, "structure") || throw(
                ArgumentError(
                    "Log `$g.$l` (array `$(stored_as[(g, l)])` of the log group) has no " *
                    "`structure` attribute naming its axes"
                )
            )
        end

        arrs = [
            DataCube(found[g][l], Symbol.(Tuple(found[g][l].attrs["structure"])))
            for l in logs
        ]
        NamedTuple{Tuple(logs)}(Tuple(arrs))
    end
    return NamedTuple{Tuple(groups)}(Tuple(loaded))
end

"""
    _log_attrs_to_props(attrs)::Dict{Symbol,Any}

Extract `units` and `description` from a Zarr array's `.attrs` dict into a
`Dict{Symbol,Any}` suitable for passing as `properties` to `DataCube`.
"""
function _log_attrs_to_props(attrs)::Dict{Symbol,Any}
    props = Dict{Symbol,Any}()
    haskey(attrs, "units") && (props[:units] = attrs["units"])
    haskey(attrs, "description") && (props[:description] = attrs["description"])
    return props
end

"""
Load shading log from a Zarr log group, handling both the new combined format
(`shading_log`) and old stores that kept `fog` and `shade` as separate arrays.
"""
function _load_Shd_log(log_set::Zarr.ZGroup)::YAXArray
    if haskey(log_set, "shading_log")
        arr = log_set["shading_log"]
        ax_names = Symbol.(Tuple(arr.attrs["structure"]))
        ax_labels = [1:s for s in size(arr)]
        return DataCube(
            arr;
            properties=_log_attrs_to_props(arr.attrs),
            NamedTuple{ax_names}(ax_labels)...
        )
    end

    # Old format: two separate (T, L, N) arrays — load and stack into (T, L, 2, N)
    fog_arr = log_set["fog"][:, :, :]
    shade_arr = log_set["shade"][:, :, :]
    tf, n_locs, n_scens = size(fog_arr)
    combined = Array{Float32}(undef, tf, n_locs, 2, n_scens)
    combined[:, :, 1, :] .= fog_arr
    combined[:, :, 2, :] .= shade_arr
    return DataCube(
        combined;
        timesteps=1:tf, locations=1:n_locs,
        intervention=["fog", "shade"], scenarios=1:n_scens
    )
end

"""
    _tag_inputs_metadata!(inputs::DataFrame, model_spec::DataFrame)::DataFrame

Tag each column of `inputs` that has a matching entry (by `fieldname`) in `model_spec`
with `"ptype"` and `"label"` column metadata (sourced from `model_spec`'s `ptype` and
`name` columns respectively). Uses `style=:note` so the metadata survives `copy()` and
column-/row-subsetting operations performed by downstream analysis code.
"""
function _tag_inputs_metadata!(inputs::DataFrame, model_spec::DataFrame)::DataFrame
    for row in eachrow(model_spec)
        col = row.fieldname
        if hasproperty(inputs, col)
            colmetadata!(inputs, col, "ptype", row.ptype; style=:note)
            colmetadata!(inputs, col, "label", row.name; style=:note)
        end
    end

    return inputs
end

function ResultSet(
    input_set::AbstractArray,
    env_layer_md::EnvLayer,
    inputs_used::DataFrame,
    outcomes::Dict,
    log_set::Zarr.ZGroup,
    dhw_stats_set::Dict,
    wave_stats_set::Dict,
    conn_data::Dict,
    loc_data::DataFrame,
    model_spec::DataFrame,
    model_specs::Vector{DataFrame}=[model_spec],
    source_sizes::Vector{Int}=[nrow(inputs_used)],
    provenance::Vector{Any}=Any[nothing]
)::ResultSet
    rcp = "RCP" in keys(input_set.attrs) ? input_set.attrs["RCP"] : input_set.attrs["rcp"]

    _tag_inputs_metadata!(inputs_used, model_spec)

    return ADRIAResultSet(
        input_set.attrs["name"],
        string(rcp),
        input_set.attrs["invoke_time"],
        input_set.attrs["ADRIA_VERSION"],
        input_set.attrs["loc_ids"],
        convert.(Float64, input_set.attrs["loc_area"]),
        convert.(Float64, input_set.attrs["loc_max_coral_cover"]),
        input_set.attrs["loc_centroids"],
        env_layer_md,
        dhw_stats_set,
        wave_stats_set,
        conn_data,
        loc_data,
        inputs_used,
        input_set.attrs["sim_constants"],
        model_spec,
        model_specs,
        source_sizes,
        provenance,
        outcomes,
        map(Float64, _rankings_data(log_set["rankings"])),
        let arr = log_set["larval_methods"]
            ax = Symbol.(Tuple(arr.attrs["structure"]))
            map(
                Float64,
                DataCube(
                    arr;
                    properties=_log_attrs_to_props(arr.attrs),
                    NamedTuple{ax}([1:s for s in size(arr)])...
                )
            )
        end,
        let arr = log_set["coral_aquaculture"]
            ax = Symbol.(Tuple(arr.attrs["structure"]))
            map(
                Float64,
                DataCube(
                    arr;
                    properties=_log_attrs_to_props(arr.attrs),
                    NamedTuple{ax}([1:s for s in size(arr)])...
                )
            )
        end,
        map(Float64, _load_Shd_log(log_set)),
        let arr = log_set["coral_dhw_log"]
            ax = Symbol.(Tuple(arr.attrs["structure"]))
            map(
                Float64,
                DataCube(
                    arr;
                    properties=_log_attrs_to_props(arr.attrs),
                    NamedTuple{ax}([1:s for s in size(arr)])...
                )
            )
        end,
        haskey(log_set, "coral_cover_log") ?
        let arr = log_set["coral_cover_log"]
            ax = Symbol.(Tuple(arr.attrs["structure"]))
            map(
                x -> Float64(x) / 65535.0,
                DataCube(
                    arr;
                    properties=_log_attrs_to_props(arr.attrs),
                    NamedTuple{ax}([1:s for s in size(arr)])...
                )
            )
        end : nothing,
        _load_submodel_logs(log_set)
    )
end

function _rankings_data(rankings_set::ZArray{T})::YAXArray{T} where {T}
    ax_names = Symbol.(Tuple(rankings_set.attrs["structure"]))
    ax_labels::Vector{Union{UnitRange{Int64},Vector{Symbol}}} = range.(
        [1], size(rankings_set)
    )

    # Replace interventions axis with named labels
    intervention_idx = findfirst(x -> x == :interventions, ax_names)
    !isnothing(intervention_idx) && (ax_labels[intervention_idx] = interventions())

    return DataCube(
        rankings_set;
        properties=_log_attrs_to_props(rankings_set.attrs),
        NamedTuple{ax_names}(ax_labels)...
    )
end

"""
    _copy_env_data(src::String, dst::String, subdir::String)::Nothing

Helper function to copy environmental data layer statistics from data store.
"""
function _copy_env_data(
    src::String, dst::String, folder_name::String, subdir=""::String
)::Nothing
    src_dir = joinpath(src, folder_name, subdir)
    dst_dir = joinpath(dst, folder_name, subdir)
    mkpath(dst_dir)
    src_ds = filter(d -> isdir(joinpath(src_dir, d)), readdir(src_dir))
    for ds in src_ds
        cp(joinpath(src_dir, ds), joinpath(dst_dir, ds); force=true)
    end

    return nothing
end

"""
    _union_inputs(inputs::Vector{<:DataFrame})::DataFrame

Stack scenario tables over the union of their columns (order of first appearance). The
stored scenario table is all `Float64`, so columns a table lacks are filled with `NaN`,
never `missing`: the result has only `Float64` columns, and an absent factor is found with
`isnan`, not `ismissing`.
"""
function _union_inputs(inputs::Vector{<:DataFrame})::DataFrame
    length(inputs) == 1 && return copy(inputs[1])

    all_inputs = reduce((a, b) -> vcat(a, b; cols=:union), inputs)
    for col in names(all_inputs)
        all_inputs[!, col] = Float64.(coalesce.(all_inputs[!, col], NaN))
    end
    return all_inputs
end

"""
    _union_model_specs(specs::Vector{<:DataFrame}, inputs::DataFrame)::DataFrame

Union of model specifications: the rows of the first, then rows of later specs whose
(`component`, `fieldname`) it does not have yet. The first source wins for the value,
distribution and description of a factor several share. The bounds and `is_constant` are
recomputed over the sources, so that no factor that varies in the merged scenario table
`inputs` is marked constant: `lower_bound` and `upper_bound` are the smallest and largest
bound of the sources holding the factor, and a factor is constant only if these are equal
and its column in `inputs` holds at most one value (`NaN`, for sources lacking the factor,
is not counted). Specifications without these columns, or a factor without bounds in any
source (`missing` after the union), are left as the first source has them.
"""
function _union_model_specs(specs::Vector{<:DataFrame}, inputs::DataFrame)::DataFrame
    length(specs) == 1 && return specs[1]

    all_specs = reduce((a, b) -> vcat(a, b; cols=:union), specs)
    merged = unique(all_specs, [:component, :fieldname])
    all(in(names(merged)), ("lower_bound", "upper_bound", "is_constant")) || return merged

    for (i, row) in enumerate(eachrow(merged))
        key = (row.component, row.fieldname)
        holding = [k == key for k in zip(all_specs.component, all_specs.fieldname)]
        lbs = collect(skipmissing(all_specs.lower_bound[holding]))
        ubs = collect(skipmissing(all_specs.upper_bound[holding]))
        (isempty(lbs) || isempty(ubs)) && continue

        lb, ub = minimum(lbs), maximum(ubs)
        merged[i, :lower_bound] = lb
        merged[i, :upper_bound] = ub

        col = string(row.fieldname)
        n_vals = col in names(inputs) ? length(unique(filter(!isnan, inputs[!, col]))) : 0
        merged[i, :is_constant] = lb == ub && n_vals <= 1
    end

    return merged
end

"""
    _has_log(rs::ADRIAResultSet, group::Symbol, log::Symbol)::Bool

Whether result set `rs` holds log `log` of submodel `group` (`rs.logs.<group>.<log>`).
"""
_has_log(rs::ADRIAResultSet, group::Symbol, log::Symbol)::Bool =
    haskey(rs.logs, group) && haskey(rs.logs[group], log)

"""
    _source_logspec(rs::ADRIAResultSet, group::Symbol, log::Symbol)::LogSpec

`LogSpec` of a submodel log of the stored result set `rs`, rebuilt from its store (a flat
array for result sets stored before logs were grouped, see `LEGACY_LOG_NAMES`).
"""
function _source_logspec(rs::ADRIAResultSet, group::Symbol, log::Symbol)::LogSpec
    log_set = zopen(joinpath(rs.env_layer_md.dpkg_path, LOG_GRP))
    haskey(log_set.groups, String(group)) &&
        return _logspec_from_store(log_set.groups[String(group)], log)

    legacy = only(k for (k, v) in pairs(LEGACY_LOG_NAMES) if v == (group, log))
    spec = _logspec_from_store(log_set, legacy)
    return LogSpec(
        log, spec.dims, spec.size;
        labels=spec.labels, units=spec.units, description=spec.description
    )
end

"""
    _merged_logspec(result_sets, group::Symbol, log::Symbol)::LogSpec

`LogSpec` of log `log` of submodel `group` in a merged result set: that of the first source
holding the log. Every source holding it must give it the same axes and shape (a
`DimensionMismatch` is thrown otherwise) and the same axis labels and units (an
`ArgumentError` otherwise). The description is taken from the first source.
"""
function _merged_logspec(result_sets, group::Symbol, log::Symbol)::LogSpec
    specs = [
        _source_logspec(rs, group, log) for rs in result_sets if _has_log(rs, group, log)
    ]
    s1 = first(specs)
    for s in specs[2:end]
        (s.dims == s1.dims && s.size == s1.size) || throw(
            DimensionMismatch(
                "Log `$group.$log` has axes $(s1.dims) of size $(s1.size) in one source " *
                "and axes $(s.dims) of size $(s.size) in another, so the sources cannot " *
                "be merged"
            )
        )
        (s.labels == s1.labels && s.units == s1.units) || throw(
            ArgumentError(
                "Log `$group.$log` has different axis labels or units in different " *
                "sources, so the sources cannot be merged"
            )
        )
    end

    return s1
end

"""
    _check_sources(rs::ADRIAResultSet)::Nothing

Check that the sources of `rs` account for its scenarios: `rs.source_sizes` adds up to the
number of rows of the scenario table, with one model specification and one provenance
record per source. Throws a `DimensionMismatch` otherwise, before anything is written.
"""
function _check_sources(rs::ADRIAResultSet)::Nothing
    n = nrow(rs.inputs)
    sum(rs.source_sizes) == n || throw(
        DimensionMismatch(
            "Result set `$(rs.name)` is inconsistent: its sources hold " *
            "$(sum(rs.source_sizes)) scenarios but the scenario table has $n"
        )
    )

    n_src = length(rs.source_sizes)
    (length(rs.model_specs) == n_src && length(rs.provenance) == n_src) || throw(
        DimensionMismatch(
            "Result set `$(rs.name)` has $n_src sources but " *
            "$(length(rs.model_specs)) model specifications and " *
            "$(length(rs.provenance)) provenance records"
        )
    )

    return nothing
end

"""
    combine_results(result_sets...)::ResultSet
    combine_results(result_set_locs::Array{String})::ResultSet

Combine arbitrary number of ADRIA result sets into a single data store.

Note: Results are stored in Zarr format. Combining data sets can be
  expected to double total disk space requirement.

# Sources with different factors or logs

- Scenario tables are combined over the union of their columns, in order of first
  appearance. A factor a source does not have is filled with `NaN` for its scenarios. The
  fill value is always `NaN`, never `missing`: the stored scenario table is all `Float64`,
  so the merged table keeps `Float64` columns and an absent factor is found with `isnan`.
- Submodel logs (`rs.logs`) are combined over the union of the (submodel, log) pairs of the
  sources and concatenated along scenarios. Logs of different submodels never mix, even if
  they have the same name. Scenarios from a source that has no such log are filled with
  `NaN`, which tells them apart from scenarios that ran with the model off (zeros). If
  every source has every log, the logs are zero-filled as usual. The fill value is chosen
  once for all logs: if any source lacks any log, every submodel log array is created with
  a `NaN` fill value, which only affects scenarios never written.
- The model specification of each source is kept (`ADRIA.model_specs(rs)`), with
  `ADRIA.scenario_sources(rs)` giving the source of each scenario. `ADRIA.model_spec(rs)` is the
  union of the specifications, the first source winning for a parameter several share,
  except for the bounds and `is_constant`: the bounds span those of every source, and a
  factor that varies in the merged scenario table is never marked constant (see
  `_union_model_specs`). The previous behaviour of keeping only the first set's
  specification is gone.
- The provenance record of each source is kept as well (`ADRIA.provenance(rs)`, aligned
  with `ADRIA.model_specs(rs)`; `nothing` for a source stored without one).

!!! warning
    A merged scenario table can hold `NaN` columns. `ADRIAanalysis` functions that take
    scenario inputs do not handle them: `rsa` (Mann-Whitney path), `pawn` and `tsa` throw an
    `ArgumentError`, while `rsa` with `method=:cramers_v` or on unordered categorical columns
    silently treats `NaN` as a level, and a constant `NaN` column is not treated as constant.
    `find_pareto_optimal`, `find_robust` and `data_envelopment_analysis` also take the
    scenario table and should be checked likewise. Subset the table (and matching outcomes)
    to the factors all sources share, or to the scenarios of one source
    (`scenario_sources`), first.
"""
function combine_results(result_sets...)::ResultSet
    return _combine_results(collect(ResultSet, result_sets))
end

"""
    _combine_results(result_sets::Vector{ResultSet})::ResultSet

[`combine_results`](@ref), for the result sets as an abstractly typed vector. The type of
an `ADRIAResultSet` depends on the submodel logs it holds, so specialising on the argument
types would compile this (large) function again for every mix of sources.
"""
function _combine_results(result_sets::Vector{ResultSet})::ResultSet
    # ── Validation ────────────────────────────────────────────────────────────
    # Check first vs last only — comparing all N*(N-1)/2 pairs is O(N) dict comparisons
    # and at ~200ms each (large dicts) adds 20+ seconds for 94 result sets.
    if length(result_sets) > 1
        @assert result_sets[1].sim_constants == result_sets[end].sim_constants
    end
    versions = Set(rs.ADRIA_VERSION for rs in result_sets)
    if length(versions) != 1
        @warn "Results were created with different versions of ADRIA so errors may occur!"
        @warn "Results from model runs < 0.9 are no longer compatible!"
    end

    foreach(_check_sources, result_sets)

    # ── Store init + inputs ───────────────────────────────────────────────────
    # Collect scalar fields from all result sets up front to avoid repeated
    # lazy-field or GC-pressure hits later. This single pass also catches any
    # per-rs metadata needed downstream.
    rs1 = result_sets[1]
    canonical_name = rs1.name
    combined_time = replace(string(now()), "T" => "_", ":" => "_", "." => "_")
    rcps = join(unique(rs.RCP for rs in result_sets), "_")
    new_loc = joinpath(
        ENV["ADRIA_OUTPUT_DIR"], "$(rs1.name)__RCPs$(rcps)__$(combined_time)"
    )
    z_store = DirectoryStore(new_loc)

    # ── Build inputs + scenario attributes ───────────────────────────────────
    n_scenarios = sum(size(rs.inputs, 1) for rs in result_sets)
    envlayer = rs1.env_layer_md
    env_md = EnvLayer(
        envlayer.dpkg_path,
        envlayer.loc_data_fn,
        envlayer.loc_id_col,
        envlayer.cluster_id_col,
        envlayer.init_coral_cov_fn,
        envlayer.connectivity_fn,
        dirname(envlayer.DHW_fn),
        dirname(envlayer.wave_fn),
        envlayer.timeframe
    )
    all_inputs = _union_inputs([getfield(rs, :inputs) for rs in result_sets])
    source_sizes = reduce(vcat, [rs.source_sizes for rs in result_sets])
    source_specs = reduce(vcat, [rs.model_specs for rs in result_sets])
    source_prov = reduce(vcat, [rs.provenance for rs in result_sets])
    input_dims = size(all_inputs)
    attrs = scenario_attributes(
        canonical_name,
        rcps,
        names(all_inputs),
        combined_time,
        env_md,
        rs1.sim_constants,
        rs1.loc_ids,
        rs1.loc_area,
        rs1.loc_max_coral_cover,
        rs1.loc_centroids
    )

    length(source_sizes) > 1 && (attrs[:source_sizes] = source_sizes)

    # One provenance record per source, aligned with `source_sizes`
    all(isnothing, source_prov) || (attrs[:provenance] = source_prov)

    # ── Copy spatial data and model spec ──────────────────────────────────────
    mkdir(joinpath(new_loc, SPATIAL_DATA))
    cp(
        attrs[:loc_data_file],
        joinpath(new_loc, SPATIAL_DATA, basename(attrs[:loc_data_file]));
        force=true
    )
    mkdir(joinpath(new_loc, MODEL_SPEC))
    CSV.write(
        joinpath(new_loc, MODEL_SPEC, "model_spec.csv"),
        _union_model_specs([rs.model_spec for rs in result_sets], all_inputs)
    )

    if length(source_sizes) > 1
        for (i, ms) in enumerate(source_specs)
            CSV.write(joinpath(new_loc, MODEL_SPEC, "model_spec_source_$(i).csv"), ms)
        end
    end

    # ── Write inputs Zarr + compute chunk size ────────────────────────────────
    input_loc::String = joinpath(z_store.folder, INPUTS)
    input_set = zcreate(
        Float64,
        input_dims...;
        fill_value=-9999.0,
        fill_as_missing=false,
        path=input_loc,
        chunks=(1, input_dims[2]),
        attrs=attrs
    )
    input_set[:, :] = Matrix(all_inputs)
    n_scenarios = nrow(all_inputs)
    env_batch::Int = parse(Int, get(ENV, "ADRIA_BATCH_SIZE", "0"))
    # Chunk size must match write granularity (scenarios per source result set).
    # Using total_scenarios/nthreads causes chunk write amplification: each small
    # write reads, decompresses, modifies, and rewrites the full chunk — O(n_sets×).
    per_rs_len::Int = minimum(size(rs.inputs, 1) for rs in result_sets)
    batch_size::Int = if env_batch > 0
        min(env_batch, n_scenarios)
    else
        per_rs_len
    end

    # ── Setup log Zarr stores ─────────────────────────────────────────────────
    n_tf = size(rs1.iv_CAq_log, :timesteps)
    n_loc = size(rs1.iv_CAq_log, :locations)

    # Submodel logs present in any source set, as (submodel, log) pairs with specs rebuilt
    # from the first source holding each log. Sources lacking a log are filled with NaN, so
    # the fill value is NaN only when some source lacks one.
    sub_pairs = sort!(
        unique([
            (g, l) for rs in result_sets for g in keys(rs.logs) for l in keys(rs.logs[g])
        ])
    )
    sub_specs = Pair{Symbol,Vector{LogSpec}}[]
    for g in unique(first.(sub_pairs))
        specs = LogSpec[
            _merged_logspec(result_sets, g, l) for (g2, l) in sub_pairs if g2 == g
        ]
        push!(sub_specs, g => specs)
    end

    sub_fill =
        any(rs -> !all(((g, l),) -> _has_log(rs, g, l), sub_pairs), result_sets) ?
        Float32(NaN) : Float32(0)
    logs = (;
        zip(
            [
                :ranks,
                :iv_LvM_log,
                :iv_CAq_log,
                :iv_Shd_log,
                :coral_dhw_tol_log,
                :coral_cover_log,
                :submodel_logs
            ],
            setup_logs(
                z_store,
                rs1.loc_ids,
                nrow(all_inputs),
                n_tf,
                n_loc,
                size(rs1.iv_CAq_log, :coral_id),
                size(rs1.coral_dhw_tol_log, :species) ÷ size(rs1.iv_CAq_log, :coral_id),
                batch_size;
                sub_specs=sub_specs,
                sub_fill=sub_fill
            )
        )...
    )

    # ── Copy logs ─────────────────────────────────────────────────────────────
    # Pre-read all source data in parallel (reads from separate Zarr stores are
    # independent), then write sequentially to each destination ZArray.
    n_rs = length(result_sets)

    # Destination arrays, with the source log each is read from: a core log, or a
    # (submodel, log) pair
    dest_logs = Pair{Any,Any}[
        k => getfield(logs, k) for k in keys(logs) if k != :submodel_logs
    ]
    for (g, arrs) in pairs(logs.submodel_logs), (l, arr) in pairs(arrs)
        push!(dest_logs, (g, l) => arr)
    end

    for (log, n_log) in dest_logs
        isnothing(n_log) && continue

        src = Vector{Union{Nothing,Array}}(undef, n_rs)
        Threads.@threads for i = 1:n_rs
            s =
                log isa Tuple ?
                (
                    _has_log(result_sets[i], log...) ? result_sets[i].logs[log[1]][log[2]] :
                    nothing
                ) :
                getfield(result_sets[i], log)
            if isnothing(s)
                src[i] = nothing
            elseif log == :coral_cover_log
                src[i] = round.(UInt16, clamp.(Array(s), 0.0, 1.0) .* 65535.0)
            else
                src[i] = Array(s)
            end
        end

        scen_id = 1
        for i = 1:n_rs
            rs_scen_len =
                isnothing(src[i]) ?
                size(result_sets[i].iv_CAq_log, :scenarios) : size(src[i], ndims(src[i]))
            if !isnothing(src[i])
                # Scenarios are the last dimension; the COTS condition log is 3-D
                lead = ntuple(_ -> Colon(), ndims(n_log) - 1)
                n_log[lead..., scen_id:(scen_id + rs_scen_len - 1)] = src[i]
            end
            scen_id += rs_scen_len
        end
    end

    # ── Copy outcome metrics ──────────────────────────────────────────────────
    # Write all metrics together per result set (one full-chunk write each) to avoid
    # the partial-chunk write amplification of the per-metric outer loop:
    # with 6 metrics in one chunk, writing metric-by-metric triggered 5 read-modify-write
    # cycles per result set per chunk.
    # Reads are parallelised; the sequential write loop is kept because loc_store is shared.
    tf_len = size(rs1.outcomes[LOC_METRIC_NAMES[1]], :timesteps)
    n_locs = size(rs1.outcomes[LOC_METRIC_NAMES[1]], :locations)
    n_loc_metrics = length(LOC_METRIC_NAMES)
    loc_dims = (tf_len, n_locs, n_loc_metrics, n_scenarios)
    loc_store = zcreate(
        Float32,
        loc_dims...;
        fill_value=Float32(0),
        fill_as_missing=false,
        path=joinpath(z_store.folder, RESULTS, LOC_METRICS),
        chunks=(loc_dims[1:3]..., batch_size),
        attrs=Dict{Symbol,Any}(
            :structure => ("timesteps", "locations", "metrics", "scenarios"),
            :unique_loc_ids => rs1.loc_ids,
            :metrics => string.(LOC_METRIC_NAMES)
        ),
        compressor=COMPRESSOR
    )
    loc_bufs = Vector{Array{Float32,4}}(undef, n_rs)
    Threads.@threads for i = 1:n_rs
        rs = result_sets[i]
        rs_scen_len = size(rs.outcomes[LOC_METRIC_NAMES[1]], :scenarios)
        buf = Array{Float32}(undef, tf_len, n_locs, n_loc_metrics, rs_scen_len)
        for (m_idx, m_name) in enumerate(LOC_METRIC_NAMES)
            buf[:, :, m_idx, :] .= Array{Float32}(rs.outcomes[m_name])
        end
        loc_bufs[i] = buf
    end
    scen_id = 1
    for i = 1:n_rs
        rs_scen_len = size(loc_bufs[i], 4)
        loc_store[:, :, :, scen_id:(scen_id + rs_scen_len - 1)] = loc_bufs[i]
        scen_id += rs_scen_len
    end

    taxa_name = :relative_taxa_cover
    taxa_dims = (size(rs1.outcomes[taxa_name])[1:(end - 1)]..., n_scenarios)
    taxa_store = zcreate(
        Float32,
        taxa_dims...;
        fill_value=Float32(0),
        fill_as_missing=false,
        path=joinpath(z_store.folder, RESULTS, string(taxa_name)),
        chunks=(taxa_dims[1:(end - 1)]..., batch_size),
        attrs=Dict{Symbol,Any}(
            :structure => ("timesteps", "groups", "scenarios")
        ),
        compressor=COMPRESSOR
    )
    taxa_bufs = Vector{Array{Float32,3}}(undef, n_rs)
    Threads.@threads for i = 1:n_rs
        taxa_bufs[i] = Array{Float32}(result_sets[i].outcomes[taxa_name])
    end
    scen_id = 1
    for i = 1:n_rs
        rs_scen_len = size(taxa_bufs[i], 3)
        taxa_store[:, :, scen_id:(scen_id + rs_scen_len - 1)] = taxa_bufs[i]
        scen_id += rs_scen_len
    end

    # ── Copy env stats + final load ───────────────────────────────────────────
    mkdir(joinpath(new_loc, ENV_STATS))
    for rs in result_sets
        loc::String = rs.env_layer_md.dpkg_path
        _copy_env_data(loc, new_loc, ENV_STATS, "dhw")
        _copy_env_data(loc, new_loc, ENV_STATS, "wave")
        _copy_env_data(loc, new_loc, CONNECTIVITY)
    end

    result = load_results(z_store.folder)
    return result
end
function combine_results(result_set_locs::Array{String})::ResultSet
    return combine_results(load_results.(result_set_locs)...)
end

"""
    env_stats(rs::ResultSet, s_name::String, rcp::String)
    env_stats(rs::ResultSet, s_name::String, rcp::String, scenario::Int)
    env_stats(rs::ResultSet, s_name::String, stat::String, rcp::String, scenario::Int)

Extract statistics for a given environmental layer ("DHW" or "wave")
"""
function env_stats(rs::ResultSet, s_name::String, rcp::String)
    return getfield(rs, Symbol("$(s_name)_stats"))[rcp]
end
function env_stats(rs::ResultSet, s_name::String, rcp::String, scenario::Int)
    return getfield(rs, Symbol("$(s_name)_stats"))[rcp][:, scenario]
end
function env_stats(rs::ResultSet, s_name::String, stat::String, rcp::String, scenario::Int)
    return getfield(rs, Symbol("$(s_name)_stats"))[rcp][stat, scenario]
end

"""
    store_name(rs::ResultSet)::String

Get name of result set.
"""
function store_name(rs::ResultSet)::String
    return "$(rs.name)__RCPs_$(rs.RCP)__$(rs.invoke_time)"
end

"""
    result_location(rs::ResultSet)::String

Get location of result set.
"""
function result_location(rs::ResultSet)::String
    store = ""
    try
        store = joinpath(ENV["ADRIA_OUTPUT_DIR"], store_name(rs))
    catch err
        if isa(err, KeyError)
            @warn "Output directory not yet set. Displaying result directory instead."
            store = store_name(rs)
        else
            rethrow(err)
        end
    end

    return replace(store, "\\" => "/")
end

"""
    select(r::ResultSet, op::String)

Hacky scenario filtering - to be replaced with more robust approach.

Only supports filtering by single attribute.
Should be expanded to support filtering metric results too.

# Examples
```julia
select(result, "guided .> 0.0")

# Above expands to:
# result.inputs.guided .> 0.0
```
"""
function select(r::ResultSet, op::String)
    scens = r.inputs

    col, qry = split(op, " "; limit=2)
    col = Symbol(col)

    df_ss = getproperty(scens, col)

    return eval(Meta.parse("$df_ss $qry"))
end

"""
    timesteps(rs::ResultSet)

Retrieve the time steps represented in the result set.

# Arguments
- `rs` : ResultSet
"""
function timesteps(rs::ResultSet)
    return rs.env_layer_md.timeframe
end

"""
    timesteps(outcomes::YAXArray)::Vector{Int64}

Extract time step labels from a YAXArray. Returns an empty `Vector{Int64}` if
the array has no `:timesteps` dimension.
"""
function timesteps(outcomes::YAXArray)::Vector{Int64}
    return :timesteps in axes_names(outcomes) ? Int.(outcomes.timesteps) : Int64[]
end

@deprecate site_k_area(rs::ResultSet) loc_k_area(rs)
@deprecate site_k(rs::ResultSet) loc_k(rs)

"""
    loc_k_area(rs::ResultSet)::Vector{Float64}

Extract vector of a location's coral carrying capacity in terms of absolute area.
"""
function loc_k_area(rs::ResultSet)::Vector{Float64}
    return rs.loc_max_coral_cover .* rs.loc_area
end

"""
    loc_k(rs::ResultSet)::Vector{Float64}

Extract vector of a location's coral carrying capacity in as a proportion relative to the
location's total area.
"""
function loc_k(rs::ResultSet)::Vector{Float64}
    return rs.loc_max_coral_cover
end

"""
    loc_area(rs::ResultSet)::Vector{Float64}

Extract vector of a location's total area in its areal unit (m², km², etc).
"""
function loc_area(rs::ResultSet)::Vector{Float64}
    return rs.loc_area
end

"""
    n_locations(rs::ResultSet)::Int64

Retrieve the number of locations represented in the result set.

# Arguments
- `rs` : ResultSet
"""
function n_locations(rs::ResultSet)::Int64
    return size(rs.loc_ids, 1)
end

"""
    n_scenarios(rs::ResultSet)::Int64

Number of scenarios stored in the given result set.

# Arguments
- `rs` : ResultSet
"""
function n_scenarios(rs::ResultSet)::Int64
    return size(rs.inputs, 1)
end

"""
    component_params(spec::DataFrame, component::Type)::DataFrame

Extract parameters for a specific model component from exported model specification.
"""
function component_params(rs::ResultSet, component::T)::DataFrame where {T}
    spec = rs.model_spec
    return spec[spec.component .== string(component), :]
end
function component_params(rs::ResultSet, components::Vector{T})::DataFrame where {T}
    spec = rs.model_spec
    return spec[spec.component .∈ [replace.(string.(components), "ADRIA." => "")], :]
end

"""
    model_spec(rs::ResultSet)::DataFrame

Extract model specification from Result Set.
"""
function model_spec(rs::ResultSet)::DataFrame
    return rs.model_spec
end

function model_spec(rs::ResultSet, param_names::Vector{Symbol})::DataFrame
    return rs.model_spec[rs.model_spec.fieldname .∈ Ref(param_names), :]
end

"""
    model_specs(rs::ADRIAResultSet)::Vector{DataFrame}

A result set made by [`combine_results`](@ref) keeps the model specification of every source
(`model_specs`, in the order the sources were given); [`model_spec`](@ref) is then their
union (first source wins for a parameter that several have). A result set that is not
merged has one source.
"""
model_specs(rs::ADRIAResultSet)::Vector{DataFrame} = rs.model_specs

"""
    scenario_sources(rs::ADRIAResultSet)::Vector{Int}

For each scenario of the result set, the index of the source it came from (see
[`model_specs`](@ref)). All ones for a result set that is not merged.
"""
scenario_sources(rs::ADRIAResultSet)::Vector{Int} = reduce(
    vcat, [fill(i, n) for (i, n) in enumerate(rs.source_sizes)]; init=Int[]
)

"""
    provenance(rs::ADRIAResultSet)::Vector{Any}

What produced each source of the result set, aligned with [`model_specs`](@ref): a
`Dict` with `adria_version` (leading `v`, as the `ADRIA_VERSION` attr), `julia_version`,
`interface_version`, `models` (attached submodels in attach order, then COTS if it was on;
each with `name`, `type`, `module`, `version`, `revision`) and `overrides` (always empty
today).

A result set that is not merged has one record. A set made by [`combine_results`](@ref)
has one record per source, in source order; merging a merged set adds its records, not one
record for the whole set. An entry is `nothing` for a source stored before provenance was
recorded.
"""
provenance(rs::ADRIAResultSet)::Vector{Any} = rs.provenance

function Base.show(io::IO, mime::MIME"text/plain", rs::ResultSet)
    vers_id = rs.ADRIA_VERSION

    tf, locs, scens = size(rs.outcomes[:relative_cover])
    # Species/size groups represented: $(species)

    rcps = join(split(rs.RCP, "_"), ", ")

    println("""
    Domain: $(rs.name)

    Run with ADRIA $(vers_id) on $(rs.invoke_time)
    Results stored at: $(result_location(rs))

    RCP(s) represented: $(rcps)
    Scenarios run: $(scens)
    Number of locations: $(locs)
    Timesteps: $(tf)

    Input layers
    ------------""")

    for fn in fieldnames(typeof(rs.env_layer_md))
        if fn == :timeframe
            tf = getfield(rs.env_layer_md, fn)
            println("$(fn) : $(tf[1]) - $(tf[end])")
            continue
        end

        println("$(fn) : $(getfield(rs.env_layer_md, fn))")
    end
end

"""
    _show_external_resultset(io::IO, rs::ResultSet; count_label::String, count::Int, n_locs::Int)::Nothing

Shared `Base.show` body for the external-model result-set interfaces (`RMEResultSet`,
`CScapeResultSet`). `count_label`/`count` describe the interface-specific scenario count.
"""
function _show_external_resultset(
    io::IO, rs::ResultSet; count_label::String, count::Int, n_locs::Int
)::Nothing
    print(
        io,
        """
        Name: $(rs.name)

        Results stored at: $(rs.env_layer_md.dpkg_path)

        RCP(s) represented: $(rs.RCP)
        $(count_label): $(count)
        Number of locations: $(n_locs)
        Timesteps: $(rs.env_layer_md.timeframe)
        """
    )
    return nothing
end
