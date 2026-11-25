"""
    _collect_scenario_results(domain, scenario, functional_groups)

Run one scenario and return all computed metric arrays without writing to disk.
Called by `run_scenario` (single write) and by the batched pmap path in `run_scenarios`.
"""
function _collect_scenario_results(
    domain::Domain,
    scenario::Union{AbstractVector,DataFrameRow},
    functional_groups::Vector{Vector{FunctionalGroup}};
    sim_cache=nothing,
    submodel_caches=nothing
)::NamedTuple
    result_set = run_model(domain, scenario, functional_groups; sim_cache, submodel_caches)
    threshold = parse(Float32, ENV["ADRIA_THRESHOLD"])

    rs_raw::Array{Float64} = result_set.raw
    lk = loc_k_area(domain)
    coral_spec::DataFrame = to_coral_spec(scenario)

    # 6 location-based metrics — order must match LOC_METRIC_NAMES
    loc_out = Array{Float32}(
        undef, size(rs_raw, 1), size(rs_raw, 4), length(LOC_METRIC_NAMES)
    )
    vals = relative_cover(rs_raw)
    vals[vals .< threshold] .= 0.0
    loc_out[:, :, 1] .= vals
    vals = relative_shelter_volume(rs_raw, lk, scenario)
    vals[vals .< threshold] .= 0.0
    loc_out[:, :, 2] .= vals
    vals = absolute_shelter_volume(rs_raw, lk, scenario)
    vals[vals .< threshold] .= 0.0
    loc_out[:, :, 3] .= vals
    vals = relative_juveniles(rs_raw)
    vals[vals .< threshold] .= 0.0
    loc_out[:, :, 4] .= vals
    vals = juvenile_indicator(rs_raw, coral_spec, lk)
    vals[vals .< threshold] .= 0.0
    loc_out[:, :, 5] .= vals
    rtc_vals = relative_loc_taxa_cover(rs_raw)
    vals = coral_evenness(rtc_vals.data)
    vals[vals .< threshold] .= 0.0
    loc_out[:, :, 6] .= vals

    # Taxa cover (groups axis, not locations)
    taxa_vals = relative_taxa_cover(rs_raw, lk)
    taxa_vals[taxa_vals .< threshold] .= 0.0

    # Fog and shading combined into a single (tf, n_locs, 2) buffer
    Shd_buf = Array{Float32}(undef, size(rs_raw, 1), size(rs_raw, 4), 2)
    Fog_vals = Matrix{Float32}(result_set.Fog_log)
    Shd_vals = Matrix{Float32}(result_set.Shd_log)
    Fog_vals[Fog_vals .< threshold] .= 0.0f0
    Shd_vals[Shd_vals .< threshold] .= 0.0f0
    Shd_buf[:, :, 1] .= Fog_vals
    Shd_buf[:, :, 2] .= Shd_vals

    # Remaining log arrays — apply threshold where possible
    site_ranks_vals = result_set.site_ranks
    try
        site_ranks_vals[site_ranks_vals .< threshold] .= Float32(0.0)
    catch err
        err isa MethodError ? nothing : rethrow(err)
    end

    LvM_vals = result_set.LvM_log
    LvM_vals[LvM_vals .< threshold] .= Float32(0.0)

    CAq_vals = result_set.CAq_log
    CAq_vals[CAq_vals .< threshold] .= Float32(0.0)

    # coral_dhw_log / coral_cover_log are false when their env var is disabled
    dhw_vals = result_set.coral_dhw_log
    cover_vals = result_set.coral_cover_log

    return (
        loc_out=loc_out,
        taxa_cover=taxa_vals,
        Shd=Shd_buf,
        site_ranks=site_ranks_vals,
        LvM_log=LvM_vals,
        CAq_log=CAq_vals,
        coral_dhw_log=dhw_vals,
        coral_cover_log=cover_vals,
        submodel_logs=_cots_log_results(result_set)
    )
end

"""
    _log_bytes_per_scen(domain::Domain, sub_specs)::NamedTuple

Uncompressed bytes one scenario adds to the result store, assuming 4 bytes per value:
`total` is the sum over the always-present arrays plus the submodel logs in `sub_specs`
(`LogSpec`s by submodel, see [`_submodel_log_specs`](@ref))
(the per-thread batch budget; the optional coral_dhw_log / coral_cover_log arrays are not
counted, so it is a lower bound) and `largest` is the largest single array, including the
optional logs when enabled by environment variable.
"""
function _log_bytes_per_scen(
    domain::Domain, sub_specs::AbstractVector{<:Pair{Symbol,<:AbstractVector{LogSpec}}}
)::NamedTuple
    n_locs = domain.coral_growth.n_locs
    n_groups = domain.coral_growth.n_groups
    n_sizes = domain.coral_growth.n_sizes
    tf = length(timesteps(domain))
    n_ivs = length(interventions())

    arrays = [
        tf * n_locs * length(LOC_METRIC_NAMES),  # loc_outcomes
        tf * n_groups,                           # relative_taxa_cover
        tf * n_locs * 2,                         # shading_log (fog + shade channels)
        tf * n_locs * n_ivs,                     # site_ranks
        tf * n_groups * n_locs,                  # LvM_log
        tf * n_groups * n_locs                   # CAq_log
    ]
    sub_sizes = [prod(s.size) for (_, specs) in sub_specs for s in specs]
    total = 4 * (sum(arrays) + sum(sub_sizes; init=0))

    # Largest single array, with the element sizes actually stored: Float32 for the
    # outcome and submodel logs, UInt16 ranks, Float16 shading, and the optional DHW
    # (Float16) and cover (UInt16) logs, which hold a single location unless enabled by
    # env var.
    sized = [
        4 * tf * n_locs * length(LOC_METRIC_NAMES),
        4 * tf * n_groups,
        2 * tf * n_locs * 2,
        2 * tf * n_locs * n_ivs,
        4 * tf * n_groups * n_locs,
        4 * tf * n_groups * n_locs
    ]
    sub_bytes = [4 * n for n in sub_sizes]

    dhw_locs = parse(Bool, get(ENV, "ADRIA_LOG_DHW_TOLS", "false")) ? n_locs : 1
    cover_locs = parse(Bool, get(ENV, "ADRIA_LOG_COVER", "false")) ? n_locs : 1
    optional = [
        2 * tf * n_groups * n_sizes * dhw_locs,
        2 * tf * n_groups * n_sizes * cover_locs
    ]
    largest = maximum(vcat(sized, sub_bytes, optional))

    return (; total, largest)
end

"""
    _max_chunk_by_size_for(domain::Domain, sub_specs)::Int

Largest number of scenarios per chunk (and per in-memory batch) within a 1.8 GB working
budget, below the hard 2 GB limit per array: see [`_chunk_cap_by_bytes`](@ref).
"""
function _max_chunk_by_size_for(
    domain::Domain, sub_specs::AbstractVector{<:Pair{Symbol,<:AbstractVector{LogSpec}}}
)::Int
    b = _log_bytes_per_scen(domain, sub_specs)
    return _chunk_cap_by_bytes(b.total, b.largest)
end

"""
    _chunk_cap_by_bytes(total::Int, largest::Int)::Int

Scenarios per chunk within the 1.8 GB working budget: the cap from the summed bytes per
scenario (the batch budget), and the separate cap from the largest array (2 GB per array).
As every counted array is part of the sum, the second binds only for arrays that are not
counted in `total` (the optional logs), and only when one really nears the limit.
"""
function _chunk_cap_by_bytes(total::Int, largest::Int)::Int
    budget = 1_800_000_000
    return max(1, min(floor(Int, budget / total), floor(Int, budget / largest)))
end

"""
    _cots_log_results(result_set)::NamedTuple

Submodel logs of a model run as Float32 arrays, grouped by submodel and named as in the
result set: `(cots = (pop = ..., condition = ...),)` (see `COTS_LOG_NAMES`). Until COTS is a
submodel its two logs are taken from the run's COTS buffers.

This is the `submodel_logs` entry of the collected scenario results, which is written to
the result store. It is not the `submodel_logs` that `run_model` returns: those are the
log buffers of the attached (hook) submodels (see [`_submodel_logs`](@ref)), which are not
stored yet.
"""
function _cots_log_results(result_set)::NamedTuple
    return (;
        COTS_LOG_GROUP => (
            pop=Float32.(result_set.cots_log),
            condition=Float32.(result_set.cots_condition_log)
        )
    )
end

"""
    _write_batch!(data_store, start_idx, results)

Write a contiguous batch of collected scenario results to the Zarr data store.
All arrays in the batch are stacked into a single in-memory buffer and written in one
Zarr call per array, so the entire batch lands in a single chunk file per array.
"""
function _write_batch!(
    data_store::NamedTuple, start_idx::Int, results::AbstractVector
)::Nothing
    n = length(results)
    idx_range = start_idx:(start_idx + n - 1)

    # loc_outcomes: (tf, n_locs, n_metrics, n)
    tf, n_locs, n_metrics = size(results[1].loc_out)
    loc_batch = Array{Float32}(undef, tf, n_locs, n_metrics, n)
    for (i, r) in enumerate(results)
        loc_batch[:, :, :, i] .= r.loc_out
    end
    data_store.loc_outcomes[:, :, :, idx_range] .= loc_batch

    # relative_taxa_cover: (tf, n_groups, n)
    _, n_groups = size(results[1].taxa_cover)
    taxa_batch = Array{Float32}(undef, tf, n_groups, n)
    for (i, r) in enumerate(results)
        taxa_batch[:, :, i] .= r.taxa_cover
    end
    data_store.relative_taxa_cover[:, :, idx_range] .= taxa_batch

    # shading_log: (tf, n_locs, 2, n)
    Shd_batch = Array{Float32}(undef, tf, n_locs, 2, n)
    for (i, r) in enumerate(results)
        Shd_batch[:, :, :, i] .= r.Shd
    end
    data_store.Shd_log[:, :, :, idx_range] .= Shd_batch

    # site_ranks: (tf, n_locs, n_interventions, n)
    if !isnothing(data_store.site_ranks)
        _, _, n_ivs = size(results[1].site_ranks)
        ranks_batch = Array{Float32}(undef, tf, n_locs, n_ivs, n)
        for (i, r) in enumerate(results)
            ranks_batch[:, :, :, i] .= r.site_ranks
        end
        data_store.site_ranks[:, :, :, idx_range] .= ranks_batch
    end

    # LvM_log and CAq_log: (tf, n_groups, n_locs, n)
    _, n_g, n_l = size(results[1].LvM_log)
    LvM_batch = Array{Float32}(undef, tf, n_g, n_l, n)
    CAq_batch = Array{Float32}(undef, tf, n_g, n_l, n)
    for (i, r) in enumerate(results)
        LvM_batch[:, :, :, i] .= r.LvM_log
        CAq_batch[:, :, :, i] .= r.CAq_log
    end
    data_store.LvM_log[:, :, :, idx_range] .= LvM_batch
    data_store.CAq_log[:, :, :, idx_range] .= CAq_batch

    # coral_dhw_log (conditional on ADRIA_LOG_DHW_TOLS)
    if parse(Bool, get(ENV, "ADRIA_LOG_DHW_TOLS", "false")) == true
        dhw_r = results[1].coral_dhw_log
        if dhw_r !== false && !isnothing(dhw_r)
            _, n_sp, n_ld = size(dhw_r)
            dhw_batch = Array{Float32}(undef, tf, n_sp, n_ld, n)
            for (i, r) in enumerate(results)
                dhw_batch[:, :, :, i] .= r.coral_dhw_log
            end
            data_store.coral_dhw_log[:, :, :, idx_range] .= dhw_batch
        end
    end

    # coral_cover_log (conditional on ADRIA_LOG_COVER)
    if parse(Bool, get(ENV, "ADRIA_LOG_COVER", "false")) == true
        cc_r = results[1].coral_cover_log
        if cc_r !== false && !isnothing(cc_r)
            _, n_sp_c, n_lc = size(cc_r)
            cc_batch = Array{Float32}(undef, tf, n_sp_c, n_lc, n)
            for (i, r) in enumerate(results)
                cc_batch[:, :, :, i] .= r.coral_cover_log
            end
            data_store.coral_cover_log[:, :, :, idx_range] .= round.(
                UInt16, clamp.(cc_batch, 0.0f0, 1.0f0) .* 65535.0f0
            )
        end
    end

    # Submodel logs (those declared for this run; scenarios axis is last)
    for (group, arrs) in pairs(data_store.submodel_logs), (name, arr) in pairs(arrs)
        r1 = results[1].submodel_logs[group][name]
        batch = Array{Float32}(undef, size(r1)..., n)
        lead = ntuple(_ -> Colon(), ndims(r1))
        for (i, r) in enumerate(results)
            batch[lead..., i] .= r.submodel_logs[group][name]
        end
        arr[lead..., idx_range] .= batch
    end

    return nothing
end

"""
    _run_batch!(batch_args, functional_groups, data_store)

Compute and write a contiguous batch of scenarios entirely on the calling worker.
Each pmap task gets one batch, so result arrays never travel back to the main process.
Workers write to non-overlapping chunk ranges, so concurrent writes are safe.
"""
function _run_batch!(
    batch_args::AbstractVector,
    functional_groups::Vector{Vector{FunctionalGroup}},
    data_store::NamedTuple
)::Nothing
    start_idx = batch_args[1][2]
    results = [
        _collect_scenario_results(args[1], args[3], functional_groups)
        for args in batch_args
    ]
    _write_batch!(data_store, start_idx, results)
    return nothing
end
