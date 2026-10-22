using CoralBlox

import CoralBlox:
    SizeClass,
    FunctionalGroup,
    reuse_buffers!,
    apply_mortality!,
    timestep!,
    coral_cover,
    max_projected_cover,
    linear_extension_scale_factors

using .metrics:
    relative_cover,
    relative_loc_taxa_cover,
    total_absolute_cover,
    absolute_shelter_volume,
    relative_shelter_volume,
    relative_juveniles,
    relative_taxa_cover,
    juvenile_indicator,
    coral_evenness

using SparseArrays
using DataStructures: CircularBuffer

using .decision

"""
    setup_cache(domain::Domain)::NamedTuple

Establish tuple of matrices/vectors for use as reusable data stores to avoid repeated memory allocations.
"""
function setup_cache(domain::Domain)::NamedTuple

    # Simulation constants
    n_locs::Int64 = domain.coral_growth.n_locs
    n_sizes::Int64 = domain.coral_growth.n_sizes
    n_groups::Int64 = domain.coral_growth.n_groups
    tf = length(timesteps(domain))

    cache = (
        # sf=zeros(n_groups, n_locs),  # stressed fecundity, commented out as it is disabled
        fec_scope=zeros(n_groups, n_locs),  # fecundity scope
        recruitment=zeros(n_groups, n_locs),  # coral recruitment
        dhw_step=zeros(n_locs),  # DHW for each time step
        eff_dhw_step=zeros(n_locs),  # effective DHW at each location's depth
        depth_coeff=zeros(n_locs),  # store for depth coefficient
        loc_area=Matrix{Float64}(loc_area(domain)'),  # area of locations
        habitable_area=Matrix{Float64}(loc_k_area(domain)'),  # location carrying capacity
        dhw_tol_mean_log=zeros(tf, n_sizes * n_groups, n_locs),  # tmp log for mean dhw tolerances
        cover_log=zeros(tf, n_sizes * n_groups, n_locs),  # tmp log for per-species cover
        # Per-scenario computation buffers: pre-allocated once per-thread and reused
        # across scenarios to avoid repeated GC-triggering allocations.
        C_cover=zeros(tf, n_groups, n_sizes, n_locs),     # main coral cover output
        bleach_dhw=zeros(tf, n_groups, n_sizes, n_locs),  # per-location bleaching DHW history
        sim_C_cover_t=zeros(n_groups, n_sizes, n_locs),   # coral cover at previous timestep
        ΔC_cover_t=zeros(n_groups, n_sizes, n_locs),      # cover change this timestep
        survival_rate_cache=ones(n_groups, n_sizes, n_locs)  # per-group/size/loc survival rate
    )

    return cache
end

"""
    _reshape_init_cover(data::AbstractMatrix{<:Union{Float32, Float64}})

Reshape initial coral cover of shape [n_groups * n_sizes ⋅ n_locations] to shape
[groups ⋅ sizes ⋅ locations]
"""
function _reshape_init_cover(
    C_cover::AbstractMatrix{T},
    dims::NTuple{3,Int64}
)::Array{T} where {T<:Union{Float32,Float64}}
    return permutedims(
        reshape(
            C_cover,
            dims
        ),
        (2, 1, 3)
    )
end

"""
    apply_survival_scaling(survival_rate::AbstractMatrix, scaling_param::AbstractVector)::AbstractMatrix
"""
function apply_survival_scaling(
    survival_rate::AbstractMatrix,
    scaling_param::AbstractVector
)::AbstractMatrix
    # If `scaling_param` is a very low negative number this function will return
    # negative numbers for lower values of `survival_rate`. For instance, when
    # `scaling_param = -1` that will happen whenever `survival_rate < 1/3`.
    # In practice, the lowest `survival_rate` value we have is 0.716, so that issue never
    # happens.
    return survival_rate .+ 0.5 .* (1 .- survival_rate) .* scaling_param
end

"""
    _to_group_size(growth_spec::CoralGrowth, data::AbstractVector{<:Union{Float32, Float64, Tuple}})::Matrix{<:Union{Float32, Float64}}

Reshape vector to shape [functional_groups ⋅ sizes]
"""
function _to_group_size(
    growth_spec::CoralGrowth, data::AbstractVector{T}
)::Matrix{T} where {T<:Union{Float32,Float64,Bool}}
    # Data is reshaped to size ⋅ groups then transposed to maintain expected order
    return Matrix(reshape(data, (growth_spec.n_sizes, growth_spec.n_groups))')
end

"""
    run_scenarios(dom::Domain, scens::DataFrame, RCP::String; show_progress=true, remove_workers=true)
    run_scenarios(dom::Domain, scens::DataFrame, RCP::Vector{String}; show_progress=true, remove_workers=true)

Run scenarios defined by the parameter table storing results to disk.
Scenarios are run in parallel where the number of scenarios > 256 for base Domain data
packages and > 4 for GBR-scale datasets.

# Notes
- Returned `Domain` holds scenario invoke time used as unique result set identifier.

# Examples
```julia
...
julia> rs_45 = ADRIA.run_scenarios(dom, scens, "45")
julia> rs_45_60 = ADRIA.run_scenarios(dom, scens, ["45", "60"])
```

# Arguments
- `dom` : Domain, to run scenarios with
- `scens` : DataFrame of scenarios to run
- `RCP` : ID or list of RCP(s) to run scenarios under.
- `show_progress` : Display progress
- `remove_workers` : If running in parallel, removes workers after completion

# Returns
ResultSet
"""
function run_scenarios(
    dom::Domain, scens::DataFrame, RCP::String; show_progress=true
)::ResultSet
    return run_scenarios(dom, scens, [RCP]; show_progress)
end
function run_scenarios(
    dom::Domain,
    scens::DataFrame,
    RCP::Vector{String};
    show_progress=true
)::ResultSet
    # Initialize ADRIA configuration options
    setup()

    # Sort RCPs so the dataframe order match the output filepath
    RCP = sort(RCP)

    @info "Running $(nrow(scens)) scenarios over $(length(RCP)) RCPs: $RCP"

    env_debug = parse(Bool, get(ENV, "ADRIA_DEBUG", "false")) == true
    @info "ADRIA_DEBUG = $env_debug"

    parallel::Bool = !env_debug && (Threads.nthreads() > 1)

    # chunk_size sets the Zarr chunk size along the scenario axis and controls how many
    # results the writer accumulates before flushing to disk. It also sets the compute
    # batching granularity (one @threads :static unit per chunk). Coarser batching means
    # fewer channel touches, which reduces channel-lock contention; finer batching
    # reduces peak in-memory results but increases channel-lock frequency. The default
    # (one chunk per thread) is the coarsest useful granularity and minimises channel
    # contention. In serial mode: fall back to the env var (default 1).
    active_threads::Int = Threads.nthreads()

    # Estimate the uncompressed size of one scenario's worth of Float32 output arrays to
    # derive the largest chunk count that keeps every Zarr chunk file below the 2 GB limit
    # enforced by BloscCompressor.  Only the always-present arrays are counted; the
    # optional coral_dhw_log / coral_cover_log arrays are excluded, so the estimate is
    # a lower bound and the cap is slightly conservative.
    _n_locs_cz = dom.coral_growth.n_locs
    _n_groups_cz = dom.coral_growth.n_groups
    _tf_cz = length(timesteps(dom))
    _n_ivs_cz = length(interventions())
    _n_metrics_cz = length(LOC_METRIC_NAMES)
    _bytes_per_scen::Int =
        4 * _tf_cz *
        (
            _n_locs_cz * _n_metrics_cz +   # loc_out
            _n_groups_cz +                  # relative_taxa_cover
            _n_locs_cz * 2 +               # shading_log (fog + Shd channels)
            _n_locs_cz * _n_ivs_cz +       # site_ranks
            _n_groups_cz * _n_locs_cz * 2  # LvM_log + CAq_log
        )
    # Stay below the hard 2 GB limit with a conservative 1.8 GB working budget
    _max_chunk_by_size::Int = max(1, floor(Int, 1_800_000_000 / _bytes_per_scen))

    env_chunk = parse(Int, get(ENV, "ADRIA_CHUNK_SIZE", "0"))
    chunk_size::Int = if env_chunk > 0
        env_chunk  # explicit override always wins
    elseif parallel
        # One chunk per thread keeps scheduling balanced; size budget is the hard cap
        min(cld(nrow(scens), active_threads), _max_chunk_by_size)
    else
        1
    end
    @info "ADRIA_CHUNK_SIZE (effective) = $chunk_size"

    # Cross product between rcps and param_df to have every row of param_df for each rcp
    rcps_df = DataFrame(; RCP=parse.(Int64, RCP))
    scenarios_df = crossjoin(scens, rcps_df)
    sort!(scenarios_df, :RCP)

    @info "Setting up Result Set"
    dom, data_store = ADRIA.setup_result_store!(dom, scenarios_df, chunk_size)

    # Convert DataFrame to named matrix for faster iteration
    scenarios_matrix::YAXArray = DataCube(
        Matrix(scenarios_df); scenarios=1:nrow(scenarios_df), factors=names(scenarios_df)
    )

    # Pre-allocate FunctionalGroup buffers so each scenario run can reuse memory.
    # Buffer contents grow over time as size classes accumulate; reusing avoids repeated
    # large allocations. In the parallel path one buffer set is created per thread to
    # avoid data races on the mutable CircularBuffer and TerminalClass fields.
    n_locs::Int64 = dom.coral_growth.n_locs
    n_sizes::Int64 = dom.coral_growth.n_sizes
    n_groups::Int64 = dom.coral_growth.n_groups
    _bin_edges::Matrix{Float64} = bin_edges(; unit=:m)

    _make_fg_buffers() = [
        FunctionalGroup.(
            eachrow(_bin_edges[:, 1:(end - 1)]),
            eachrow(_bin_edges[:, 2:end]),
            eachrow(zeros(n_groups, n_sizes))
        )
        for _ = 1:n_locs
    ]

    if parallel
        # One buffer set per thread, indexed directly by threadid() — no Channel/lock
        # involved. This is only safe because the scenario loop below also uses
        # :static scheduling, which pins each iteration to a fixed thread for its
        # whole lifetime (unlike :dynamic, whose tasks may migrate between threads,
        # which is why a Channel-based checkout pool was needed previously). Buffer
        # sets are allocated in parallel since each set is large (n_locs × n_groups
        # × n_sizes CircularBuffers) and independent.
        # Use maxthreadid() so the buffer vector covers all thread pools
        # (default + interactive). Threads.nthreads() only counts the default
        # pool, but Threads.threadid() returns a global ID that can exceed that
        # count when interactive-pool threads are present (e.g. --threads=6,2).
        n_thread_bufs = Threads.maxthreadid()
        fg_bufs = Vector{typeof(_make_fg_buffers())}(undef, n_thread_bufs)
        sim_caches = Vector{NamedTuple}(undef, n_thread_bufs)
        Threads.@threads :static for i = 1:Threads.nthreads()
            fg_bufs[i] = _make_fg_buffers()
            sim_caches[i] = setup_cache(dom)
        end
        # Initialize remaining slots (interactive-pool threads) serially.
        for i = (Threads.nthreads() + 1):n_thread_bufs
            fg_bufs[i] = _make_fg_buffers()
            sim_caches[i] = setup_cache(dom)
        end

        # Prevent BLAS internal threads from competing with Julia threads.
        BLAS.set_num_threads(1)

        for rcp in RCP
            run_msg = "Running $(nrow(scens)) scenarios for RCP $rcp"

            dom = switch_RCPs!(dom, rcp)
            target_rows = findall(
                scenarios_matrix[factors = At("RCP")] .== parse(Float64, rcp)
            )
            scen_args = collect(
                _scenario_args(dom, scenarios_matrix, rcp, length(target_rows))
            )
            N_rcp = length(scen_args)
            n_chunks = cld(N_rcp, chunk_size)

            prog =
                show_progress ?
                ProgressMeter.Progress(n_chunks; desc=run_msg, dt=4) :
                nothing

            # Batch accumulation: each :static thread collects chunk_size results
            # serially before a single channel put. This keeps channel-lock touches at
            # their minimum (1 put per thread = n_chunks total) and avoids any
            # back-pressure stalls that would arise from per-scenario puts filling a
            # bounded channel buffer. The writer task is the sole consumer; with the
            # channel buffer sized at 2×nthreads all chunks can queue without blocking.
            ch = Channel{Tuple{Int,Vector{NamedTuple}}}(2 * Threads.nthreads())

            chunk_len(c) =
                c < n_chunks - 1 ? chunk_size : N_rcp - (n_chunks - 1) * chunk_size

            writer = Threads.@spawn begin
                for (start_idx, results) in ch
                    _write_batch!(data_store, start_idx, results)
                    isnothing(prog) || ProgressMeter.next!(prog)
                end
            end

            try
                Threads.@threads :static for c = 0:(n_chunks - 1)
                    lo = c * chunk_size + 1
                    n_c = chunk_len(c)
                    tid = Threads.threadid()
                    fg = fg_bufs[tid]
                    sc = sim_caches[tid]
                    results = [
                        _collect_scenario_results(
                            scen_args[lo + k][1], scen_args[lo + k][3], fg;
                            sim_cache=sc
                        ) for k = 0:(n_c - 1)
                    ]
                    put!(ch, (scen_args[lo][2], results))
                end
            finally
                close(ch)
            end
            wait(writer)
        end
    else
        # Serial: one shared buffer set, collect batch then write.
        functional_groups = _make_fg_buffers()

        for rcp in RCP
            run_msg = "Running $(nrow(scens)) scenarios for RCP $rcp"

            dom = switch_RCPs!(dom, rcp)
            scen_args = collect(
                _scenario_args(dom, scenarios_matrix, rcp, size(scenarios_matrix, 1))
            )
            all_batches = [
                collect(b) for b in Iterators.partition(scen_args, chunk_size)
            ]

            if show_progress
                @showprogress desc = run_msg dt = 4 for batch in all_batches
                    batch_results = map(
                        dfx -> _collect_scenario_results(dfx[1], dfx[3], functional_groups),
                        batch
                    )
                    _write_batch!(data_store, first(batch)[2], batch_results)
                end
            else
                for batch in all_batches
                    batch_results = map(
                        dfx -> _collect_scenario_results(dfx[1], dfx[3], functional_groups),
                        batch
                    )
                    _write_batch!(data_store, first(batch)[2], batch_results)
                end
            end
        end
    end

    return load_results(_result_location(dom, RCP))
end

function growth_acceleration(
    height::Float64, midpoint::Float64, steepness::Float64, available_space::Float64
)
    return height / (1 + exp(-steepness * (available_space - midpoint))) + 1.0
end

function _scenario_args(dom, scenarios_matrix::YAXArray, rcp::String, n::Int)
    target_rows = findall(
        collect(scenarios_matrix[factors = At("RCP")]) .== parse(Float64, rcp)
    )
    rep_doms = Iterators.repeated(dom, n)
    return zip(
        rep_doms, target_rows, eachrow(scenarios_matrix[target_rows, :])
    )
end

"""
    _collect_scenario_results(domain, scenario, functional_groups)

Run one scenario and return all computed metric arrays without writing to disk.
Called by `run_scenario` (single write) and by the batched pmap path in `run_scenarios`.
"""
function _collect_scenario_results(
    domain::Domain,
    scenario::Union{AbstractVector,DataFrameRow},
    functional_groups::Vector{Vector{FunctionalGroup}};
    sim_cache=nothing
)::NamedTuple
    result_set = run_model(domain, scenario, functional_groups; sim_cache)
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
        cots_pop_log=Float32.(result_set.cots_log),
        cots_condition_log=Float32.(result_set.cots_condition_log)
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

    # COTS population log: (tf, 3, n_locs, n)
    cots_r1 = results[1].cots_pop_log
    cots_tf, cots_ages, cots_nlocs = size(cots_r1)
    cots_pop_batch = Array{Float32}(undef, cots_tf, cots_ages, cots_nlocs, n)
    for (i, r) in enumerate(results)
        cots_pop_batch[:, :, :, i] .= r.cots_pop_log
    end
    data_store.cots_pop_log[:, :, :, idx_range] .= cots_pop_batch

    # COTS body condition log: (tf, n_locs, n)
    cots_bc_r1 = results[1].cots_condition_log
    cots_bc_tf, cots_bc_nlocs = size(cots_bc_r1)
    cots_bc_batch = Array{Float32}(undef, cots_bc_tf, cots_bc_nlocs, n)
    for (i, r) in enumerate(results)
        cots_bc_batch[:, :, i] .= r.cots_condition_log
    end
    data_store.cots_condition_log[:, :, idx_range] .= cots_bc_batch

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

"""
    run_scenario(domain::Domain, idx::Int64, scenario::Union{AbstractVector,DataFrameRow}, functional_groups::Vector{Vector{FunctionalGroup}}, data_store::NamedTuple)::Nothing
    run_scenario(domain::Domain, scenario::Union{AbstractVector,DataFrameRow})::NamedTuple
    run_scenario(domain::Domain, scenario::Union{AbstractVector,DataFrameRow}, RCP::String)::NamedTuple

Run individual scenarios for a given domain, saving results to a Zarr data store.
Results are stored in Zarr format at a pre-configured location.
Sets up a new `cache` if not provided.

# Arguments
- `domain` : Simulation domain (may be modified via `switch_RCPs!`).
- `idx` : Scenario index, to store results into `data_store`.
- `scenario` : Parameter row describing the scenario.
- `functional_groups` : Preallocated functional group buffers.
- `data_store` : Pre-opened store with arrays to write results into.

# Returns
Nothing
"""
function run_scenario(
    domain::Domain,
    idx::Int64,
    scenario::Union{AbstractVector,DataFrameRow},
    functional_groups::Vector{Vector{FunctionalGroup}},
    data_store::NamedTuple
)::Nothing
    if domain.RCP == ""
        local rcp
        try
            rcp = string(Int64(scenario[At("RCP")]))
        catch err
            if !(err isa MethodError)
                rethrow(err)
            end

            rcp = scenario.RCP  # Extract from dataframe
        end

        domain = switch_RCPs!(domain, rcp)
    end

    result = _collect_scenario_results(domain, scenario, functional_groups)
    _write_batch!(data_store, idx, [result])
    return nothing
end
function run_scenario(
    domain::Domain, scenario::Union{AbstractVector,DataFrameRow}
)::NamedTuple
    return run_model(domain, scenario)
end
function run_scenario(
    domain::Domain, scenario::Union{AbstractVector,DataFrameRow}, RCP::String
)::NamedTuple
    domain = switch_RCPs!(domain, RCP)
    return run_scenario(domain, scenario)
end

"""
    run_model(domain::Domain, param_set::Union{DataFrameRow,YAXArray})::NamedTuple
    run_model(domain::Domain, param_set::DataFrameRow, functional_groups::Vector{Vector{FunctionalGroup}})::NamedTuple
    run_model(domain::Domain, param_set::YAXArray, functional_groups::Vector{Vector{FunctionalGroup}})::NamedTuple

Core scenario running function. When called with only `domain` and `param_set` arguments,
`ADRIA.setup()` must be run beforehand.

# Arguments
- `apply_allee_effect` : Whether fecundity is suppressed at low density (see `allee_effect()`
    in `ecosystem/corals/growth.jl`). Defaults to `true`.

# Returns
NamedTuple of collated results
- `raw` : Array, Coral cover relative to k area
- `CAq_log` : Array, Log of coral aquaculture locations
- `LvM_log` : Array, Log of larval methods locations
- `Fog_log` : Array, Log of fogged locations
- `Shd_log` : Array, Log of shaded locations
- `site_ranks` : Array, Log of location rankings
- `bleaching_mortality` : Array, Log of bleaching DHW (DHW-weeks) per location/group/size class — used as the lower bound of the tolerance distribution in the following timestep
- `coral_dhw_log` : Array, Log of DHW tolerances / adaptation over time (only logged in debug mode)
"""
function run_model(
    domain::Domain,
    param_set::Union{DataFrameRow,YAXArray};
    apply_allee_effect::Bool=true
)
    n_locs::Int64 = domain.coral_growth.n_locs
    n_sizes::Int64 = domain.coral_growth.n_sizes
    n_groups::Int64 = domain.coral_growth.n_groups

    edges = bin_edges()
    bin_start = edges[:, 1:(end - 1)]
    bin_end = edges[:, 2:end]
    functional_groups::Vector{Vector{FunctionalGroup}} = Vector(undef, n_locs)
    group_info = zeros(n_groups, n_sizes)
    for i = 1:n_locs
        fg::Vector{FunctionalGroup} = Vector(undef, n_groups)
        for g = 1:n_groups
            fg[g] = FunctionalGroup(
                bin_start[g, :],
                bin_end[g, :],
                group_info[g, :]
            )
        end

        functional_groups[i] = fg
    end

    return run_model(
        domain, param_set, functional_groups; apply_allee_effect=apply_allee_effect
    )
end
function run_model(
    domain::Domain,
    param_set::DataFrameRow,
    functional_groups::Vector{Vector{FunctionalGroup}};
    sim_cache=nothing,
    apply_allee_effect::Bool=true
)::NamedTuple
    setup()
    ps = DataCube(Vector(param_set); factors=names(param_set))
    return run_model(
        domain,
        ps,
        functional_groups;
        sim_cache=sim_cache,
        apply_allee_effect=apply_allee_effect
    )
end
function run_model(
    domain::Domain,
    param_set::YAXArray,
    functional_groups::Vector{Vector{FunctionalGroup}};
    sim_cache=nothing,
    apply_allee_effect::Bool=true
)::NamedTuple
    corals = to_coral_spec(param_set)
    cache = isnothing(sim_cache) ? setup_cache(domain) : sim_cache
    if !isnothing(sim_cache)
        # Reset log accumulators for this new scenario (they are written to each timestep)
        fill!(cache.dhw_tol_mean_log, 0.0)
        fill!(cache.cover_log, 0.0)
    end

    factor_names::Vector{String} = collect(param_set.factors.val)

    # Sim constants
    sim_params = domain.sim_constants
    n_locs::Int64 = domain.coral_growth.n_locs
    n_groups::Int64 = domain.coral_growth.n_groups
    n_sizes::Int64 = domain.coral_growth.n_sizes

    rng = set_random_seed(param_set)

    # Extract environmental data
    dhw_idx::Int64 = Int64(param_set[At("dhw_scenario")])
    if dhw_idx > 0.0
        if has_MCB_scenarios(domain.dhw_scens)
            # Cast to axis type to prevent Float32/Float64 mismatch SelectorErrors
            MCB_albedo = eltype(domain.dhw_scens.albedo)(param_set[At("iv_MCB_albedo")])
            MCB_duration = eltype(domain.dhw_scens.mcb_durations)(
                param_set[At("iv_MCB_duration")]
            )
            MCB_freq = Int64(param_set[At("iv_MCB_deployment_freq")])

            # Hardcode MCB start year to 2035
            MCB_start_year = findfirst(domain.env_layer_md.timeframe .== 2035)
            if isnothing(MCB_start_year)
                MCB_start_year = 1
                @warn "MCB start year 2035 not found in timeframe. Defaulting to first year."
            end

            # Get baseline (0-day) and treated slices
            # Slicing results in (timesteps, locations)
            dhw_baseline = @view(
                domain.dhw_scens[
                    scenarios = dhw_idx,
                    mcb_durations = 1,
                    albedo = At(domain.dhw_scens.albedo[1])
                ]
            )

            # Default treated to baseline. If MCB is active, slice the treated array.
            dhw_treated = dhw_baseline
            if MCB_duration > 0.0 && MCB_albedo > 0.0
                dhw_treated = @view(
                    domain.dhw_scens[
                        scenarios = dhw_idx,
                        mcb_durations = At(MCB_duration),
                        albedo = At(MCB_albedo)
                    ]
                )
            end

            tf::Int64 = size(dhw_baseline, 1)

            # Create hybrid DHW environment (Temporal Splicing)
            dhw_scen = copy(dhw_baseline)
            MCB_years = MCB_start_year:tf
            MCB_loc_mask = domain.loc_data.UNIQUE_ID .∈ [domain.Shd_target_locations]
            if !isempty(MCB_years)
                MCB_active_years = decision_frequency(
                    MCB_start_year, tf, length(MCB_years), MCB_freq
                )
                for t = 1:tf
                    if MCB_active_years[t]
                        dhw_scen[t, MCB_loc_mask] .= dhw_treated[t, MCB_loc_mask]
                    end
                end
            end
        else
            dhw_scen = @view(domain.dhw_scens[:, :, dhw_idx])
            tf = size(dhw_scen, 1)
        end
    else
        # Run with no DHW disturbances
        dhw_scen = if has_MCB_scenarios(domain.dhw_scens)
            copy(domain.dhw_scens[:, :, 1, 1, 1])
        else
            copy(domain.dhw_scens[:, :, 1])
        end
        dhw_scen .= 0.0
        tf = size(dhw_scen, 1)
    end

    wave_idx::Int64 = Int64(param_set[At("wave_scenario")])
    if wave_idx > 0.0
        # TODO: Better conversion of Ub to wave mortality
        #       Currently scaling significant wave height by its max to non-dimensionalize values
        wave_scen = copy(domain.wave_scens[:, :, wave_idx])
        wave_scen .= wave_scen ./ maximum(wave_scen)
        replace!(wave_scen, Inf => 0.0, NaN => 0.0)
    else
        # Run with no wave disturbances
        wave_scen = copy(domain.wave_scens[:, :, 1])
        wave_scen .= 0.0
    end

    cyclone_mortality_idx::Int64 = Int64(param_set[At("cyclone_mortality_scenario")])
    if cyclone_mortality_idx > 0.0
        cyclone_mortality_scen = @view(
            domain.cyclone_mortality_scens[:, :, :, cyclone_mortality_idx]
        )
    else
        # Run with no cyclone disturbances
        cyclone_mortality_scen = copy(domain.cyclone_mortality_scens[:, :, :, 1])
        cyclone_mortality_scen .= 0.0
    end
    # Environment variables are stored as strings, so convert to bool for use
    in_debug_mode = parse(Bool, get(ENV, "ADRIA_DEBUG", "false")) == true
    log_dhw_tols = parse(Bool, get(ENV, "ADRIA_LOG_DHW_TOLS", "true")) == true
    log_cover = parse(Bool, get(ENV, "ADRIA_LOG_COVER", "false")) == true
    # Opt-in minimum bleaching susceptibility threshold (DHW-weeks). Unset keeps HEAT_LB,
    # the established behaviour; see bleaching_mortality! and the sandbox calibration
    # protocol (2026-10-14).
    bleaching_dhw_floor = parse(Float64, get(ENV, "ADRIA_BLEACHING_DHW_FLOOR", string(HEAT_LB)))
    isfinite(bleaching_dhw_floor) && bleaching_dhw_floor >= HEAT_LB ||
        error("ADRIA_BLEACHING_DHW_FLOOR must be finite and >= HEAT_LB ($HEAT_LB)")
    # Opt-in scaling of bleaching mortality below a DHW threshold (unset: established path).
    bleaching_low_dhw_scale = parse(Float64, get(ENV, "ADRIA_BLEACHING_LOW_DHW_SCALE", "1.0"))
    bleaching_low_dhw_threshold = parse(Float64, get(ENV, "ADRIA_BLEACHING_LOW_DHW_THRESHOLD", "4.0"))
    0.0 <= bleaching_low_dhw_scale <= 1.0 && isfinite(bleaching_low_dhw_threshold) ||
        error("ADRIA_BLEACHING_LOW_DHW_SCALE must lie in [0, 1] with a finite threshold")

    # Initialize cover loss tracking for reactive strategies
    max_lookback = Int64(param_set[At("reactive_response_delay")])
    if max_lookback < 1
        recent_cover_losses = CircularBuffer{Vector{Float64}}(1)
    else
        recent_cover_losses = CircularBuffer{Vector{Float64}}(max_lookback)
    end
    push!(recent_cover_losses, zeros(n_locs))

    # Locations to intervene
    min_iv_locs::Int64 = param_set[At("min_iv_locations")]
    LvM_min_iv_locs::Int64 = param_set[At("iv_LvM_min_iv_locations")]

    fogging::Float64 = param_set[At("iv_Fog")]
    Shd::Float64 = param_set[At("iv_Shd")]  # DHW equivalents reduced by some shading mechanism
    Shd_years::Int64 = param_set[At("iv_Shd_years")]  # number of years with shading

    # Years to start coral aquaculture/shading/fogging
    Shd_start_year::Int64 = param_set[At("iv_Shd_year_start")]

    colony_areas = _to_group_size(
        domain.coral_growth, colony_mean_area(corals.mean_colony_diameter_m)
    )

    habitable_areas::Matrix{Float64} = cache.habitable_area
    fecundity_per_m²::Matrix{Float64} = _to_group_size(
        domain.coral_growth, corals.fecundity
    ) # number of larvae produced per m²

    # Caches
    conn = sparse(domain.conn.data)

    # Determine contribution of each source to a sink location
    # i.e., columns should sum to 1!
    TP_data = conn ./ sum(conn; dims=1)
    replace!(nonzeros(TP_data), NaN => 0.0)
    dropzeros!(TP_data)

    # sf = cache.sf  # unused as it is currently deactivated
    fec_scope::Matrix{Float64} = cache.fec_scope
    fill!(fec_scope, 0.0)
    recruitment::Matrix{Float64} = cache.recruitment
    fill!(recruitment, 0.0)
    dhw_t::Vector{Float64} = cache.dhw_step
    eff_dhw_t::Vector{Float64} = cache.eff_dhw_step
    depth_coeff::Vector{Float64} = cache.depth_coeff

    # Used to distribute larval methods settlers
    prop_fecundity::Matrix{Float64} = copy(cache.fec_scope)

    loc_data = domain.loc_data
    depth_coeff .= depth_coefficient.(loc_data.depth_med)

    # Coral cover relative to available area (i.e., 1.0 == site is filled to max capacity)
    C_cover = cache.C_cover
    fill!(C_cover, 0.0)
    C_cover[1, :, :, :] .= _reshape_init_cover(
        domain.init_coral_cover, (n_sizes, n_groups, n_locs)
    )

    # Locations that can support corals
    vec_abs_k = loc_k_area(domain)
    habitable_locs::BitVector = location_k(domain) .> 0.0
    habitable_loc_areas = vec_abs_k[habitable_locs]
    habitable_loc_areas′ = reshape(habitable_loc_areas, (1, 1, length(habitable_locs)))
    n_habitable_locs::Int64 = length(habitable_locs)

    # Avoid placing importance on sites that were not considered
    # Lower values are higher importance/ranks.
    # Values of n_locs+1 indicate locations that were not considered in rankings.
    log_location_ranks = ZeroDataCube(;     # log coral aquaculture/fogging ranks
        T=Float64,
        timesteps=1:tf,
        locations=domain.loc_ids,
        intervention=interventions()
    )

    YShd = spzeros(tf, n_locs)
    Yfog = spzeros(tf, n_locs)
    YCAq = zeros(tf, n_groups, n_locs)  # 3 = the number of coral aquaculture coral types
    YLvM = zeros(tf, n_groups, n_locs)

    # Prep scenario-specific flags/values
    # Intervention strategy: < 0 is no intervention, 0 is random location selection, > 0 is guided
    is_guided = param_set[At("guided")] > 0
    if is_guided
        MCDA_approach = mcda_methods()[Int64(param_set[At("mcda_method")])]
    end

    # Number of time steps in environmental layers to look ahead when making decisions
    plan_horizon::Int64 = Int64(param_set[At("plan_horizon")])

    # Decisions should place more weight on environmental conditions closer to
    # the decision point. `projection_confidence` shapes how aggressively that
    # weighting decays with lead time (see decision.build_decay).
    projection_confidence::Float64 = param_set[At("projection_confidence")]
    decay = build_decay(plan_horizon, projection_confidence)

    # Years at which intervention locations are re-evaluated and deployed
    # CAq_decision_years = decision_frequency(
    #     CAq_start_year, tf, CAq_years, param_set[At("iv_CAq_deployment_freq")]
    # )
    # fog_decision_years = decision_frequency(
    #     fog_start_year, tf, fog_years, param_set[At("iv_Fog_deployment_freq")]
    # )
    iv_CAq_last_deployment = zeros(Int64, n_locs)
    last_Fog_deployment = zeros(Int64, n_locs)
    last_LvM_deployment = zeros(Int64, n_locs)
    Shd_decision_years = decision_frequency(
        Shd_start_year, tf, Shd_years, param_set[At("iv_Shd_deployment_freq")]
    )

    # Define taxa and size class for coral aquaculture, and identify their factor names
    # TODO: Coral-aquaculture-deploy 1-year old corals!!! If this is the 1st size class, that's fine but needs
    # to be confirmed with ecoRRAP
    _CAq_size_groups::BitMatrix = CAq_size_groups(n_groups, n_sizes)

    # Set up assisted adaptation values
    a_adapt::Vector{Float64} = fill(param_set[At("iv_CAq_a_adapt")], n_groups)

    CAq_idx = findall(n -> contains(n, "iv_CAq_N"), factor_names)
    CAq_volume = view(param_set.data, CAq_idx)
    CAq_devices_per_m2::Float64 = param_set[At("iv_CAq_devices_per_m2")]

    is_unguided = param_set[At("guided")] == 0.0
    is_CAq = any(>(0), CAq_volume)

    # Flag indicating whether to deploy coral aquaculture when unguided
    unguided_CAq = is_unguided && is_CAq

    # Flag indicating whether to fog or not fog
    is_fogging = fogging > 0.0
    unguided_fogging = is_unguided && is_fogging

    # Larval methods flag
    iv_LvM_N_settlers = param_set[At("iv_LvM_N_settlers")]
    is_LvM = iv_LvM_N_settlers > 0.0
    unguided_LvM = is_unguided && is_LvM

    # Flag indicating whether to apply shading
    apply_Shd = Shd > 0.0

    depth_criteria = identify_within_depth_bounds(
        loc_data.depth_med, param_set[At("depth_min")], param_set[At("depth_offset")]
    )

    CAq_t_locs = eltype(domain.CAq_target_locations) <: AbstractString ? String[] :
        vcat(getproperty.(domain.CAq_target_locations, :target_locs)...)
    LvM_t_locs = eltype(domain.LvM_target_locations) <: AbstractString ? String[] :
        vcat(getproperty.(domain.LvM_target_locations, :target_locs)...)
    Shd_locs_mask::BitVector = domain.loc_ids .∈ [domain.Shd_target_locations]

    if is_guided
        CAq_pref, CAq_decision_mat, CAq_strategy = setup_guided_intervention(
            domain, param_set, depth_criteria, CAqPreferences, CAq_t_locs, is_CAq,
            build_CAq_strategy
        )
        Fog_pref, Fog_decision_mat, Fog_strategy = setup_guided_intervention(
            domain, param_set, depth_criteria, FogPreferences, domain.Fog_target_locations,
            is_fogging, build_Fog_strategy
        )
        LvM_pref, LvM_decision_mat, LvM_strategy = setup_guided_intervention(
            domain, param_set, depth_criteria, LvMPreferences, LvM_t_locs,
            is_LvM, build_LvM_strategy
        )
    else
        CAq_strategy =
            unguided_CAq ? build_CAq_strategy(param_set, domain, CAq_t_locs) : nothing
        Fog_strategy =
            unguided_fogging ?
            build_Fog_strategy(param_set, domain, domain.Fog_target_locations) :
            nothing
        LvM_strategy =
            unguided_LvM ? build_LvM_strategy(param_set, domain, LvM_t_locs) : nothing
    end

    dhw_projection::Vector{Float64} = zeros(Float64, n_locs)
    wave_projection::Vector{Float64} = zeros(Float64, n_locs)

    # Set up distributions for natural adaptation/heritability
    c_mean_t::Array{Float64,3} = repeat(
        _to_group_size(domain.coral_growth, corals.dist_mean),
        1,
        1,
        n_locs
    )
    c_std::Array{Float64,2} = _to_group_size(
        domain.coral_growth, corals.dist_std
    )

    # Store means before mortality and growth events to use during update
    c_mean_t_1 = copy(c_mean_t)

    # Hard ceiling: tolerance cannot increase more than HEAT_UB DHW-weeks above initial values
    c_mean_tol_ceil = c_mean_t .+ HEAT_UB

    # Projection coral switches (crash_mechanism_protocol.md, 2026-10-21). Unset = established path.
    # M4: mean DHW tolerance of every group and size class rises by `rate` DHW-weeks per year from
    # the start year, under the HEAT_UB ceiling. M5: from the start year every habitable location
    # receives `fraction` x the domain-mean larval density per group of the earlier years.
    thermal_acclimation_rate = parse(Float64, get(ENV, "ADRIA_THERMAL_ACCLIMATION_RATE", "0.0"))
    thermal_acclimation_start = parse(Int, get(ENV, "ADRIA_THERMAL_ACCLIMATION_START_YEAR", "0"))
    ext_coral_fraction = parse(Float64, get(ENV, "ADRIA_EXTERNAL_CORAL_LARVAE_FRACTION", "0.0"))
    ext_coral_start = parse(Int, get(ENV, "ADRIA_EXTERNAL_CORAL_LARVAE_START_YEAR", "0"))
    coral_calendar_years = Int.(domain.env_layer_md.timeframe[1:tf])
    thermal_acclimation_rate >= 0.0 && ext_coral_fraction >= 0.0 ||
        error("Thermal acclimation rate and external coral larvae fraction must be non-negative")
    (thermal_acclimation_rate == 0.0 || thermal_acclimation_start >= 1000) &&
        (ext_coral_fraction == 0.0 || ext_coral_start >= 1000) &&
        (thermal_acclimation_rate == 0.0 && ext_coral_fraction == 0.0 || all(>=(1000), coral_calendar_years)) ||
        error("Coral projection switches need a calendar start year and a calendar-year timeframe")
    ext_coral_sum = zeros(Float64, n_groups)  # domain-mean larvae m^-2 per group, summed over years
    ext_coral_n = Ref(0)
    ext_coral_larvae = ext_coral_fraction > 0.0 ? zeros(Float64, n_groups, n_locs) : nothing

    # Snapshot of juvenile tolerance AFTER natural adaptation but BEFORE coral aquaculture.
    # Used as the reference base for a_adapt enhancement (c_mean_reference).
    # Storing the post-coral-aquaculture state would cause a_adapt to compound year-over-year,
    # driving unrealistic tolerance values (>100 DHW-weeks over long runs).
    c_mean_pre_CAq = zeros(n_groups, n_locs)

    # Log of distributions
    dhw_tol_mean_log = cache.dhw_tol_mean_log  # tmp log for mean dhw tolerances
    cover_log = cache.cover_log  # tmp log for per-species cover

    # Cache for per-location bleaching DHW (DHW-weeks) used as the lower bound of the
    # tolerance distribution in the following timestep. Sliding window over two timesteps:
    # dim-1 index 1 = previous timestep's bleaching DHW, 2 = current (being written).
    bleach_dhw = cache.bleach_dhw
    fill!(bleach_dhw, 0.0)
    #### End coral constants

    ## Update ecological parameters based on intervention option
    # Pre-calculate proportion of survivers from wave stress
    # Sw_t = wave_damage!(cache.wave_damage, wave_scen, corals.wavemort90, n_species)

    area_weighted_conn = sparse(conn .* vec_abs_k)
    conn_cache = copy(area_weighted_conn)

    # basal_area_per_settler is the area in m^2 of a size class one coral
    basal_area_per_settler = colony_mean_area(
        corals.mean_colony_diameter_m[corals.class_id .== 1]
    )

    # Dummy vars to fill/replace with ranks of selected locations
    selected_CAq_ranks = []
    selected_Fog_ranks = []
    selected_LvM_ranks = []

    # Cache matrix to store potential settlers
    potential_settlers = zeros(size(fec_scope)...)
    _linear_extensions = _to_group_size(domain.coral_growth, corals.linear_extension)
    _bin_edges = bin_edges()
    survival_rate = 1.0 .- _to_group_size(domain.coral_growth, corals.mb_rate)

    # Empty the old contents of the buffers and add the new blocks
    cover_view = [@view C_cover[1, :, :, loc] for loc = 1:n_locs]
    functional_groups::Vector{Vector{FunctionalGroup}} = reuse_buffers!.(
        functional_groups, (cover_view .* vec_abs_k)
    )

    # Reuse pre-allocated buffers from cache to avoid per-scenario GC pressure
    survival_rate_cache = cache.survival_rate_cache
    fill!(survival_rate_cache, 1.0)
    C_cover_t = cache.sim_C_cover_t
    fill!(C_cover_t, 0.0)
    ΔC_cover_t = cache.ΔC_cover_t
    fill!(ΔC_cover_t, 0.0)

    # Max projected cover is used on linear extensions scale factors
    habitable_max_projected_cover = max_projected_cover(
        _linear_extensions,
        _bin_edges,
        habitable_loc_areas
    )

    # Extract unique cb_calib_groups and create location masks for each biogroup
    unique_cb_calib_groups::Vector{Int64} = sort(unique(domain.loc_data.CB_CALIB_GROUPS))
    
    # Get total calibration groups from params to ensure consistent sizes
    lin_ext_data = param_set[occursin.("linear_extension_scale", factor_names)].data
    total_calib_groups = div(length(lin_ext_data), n_groups)

    cb_calib_group_masks::BitMatrix = falses(n_locs, total_calib_groups)
    for cb_calib_group in unique_cb_calib_groups
        cb_calib_group_masks[:, cb_calib_group] .= domain.loc_data.CB_CALIB_GROUPS .== cb_calib_group
    end

    # Index into unique_biogroups for each location
    loc_cb_calib_group_idxs::Vector{Int64} = domain.loc_data.CB_CALIB_GROUPS

    is_growth_acc_mask = occursin.("growth_acceleration", factor_names)
    growth_acc_steepness::Vector{Float64} =
        param_set[is_growth_acc_mask .&& occursin.("steepness", factor_names)].data
    growth_acc_height::Vector{Float64} =
        param_set[is_growth_acc_mask .&& occursin.("height", factor_names)].data
    growth_acc_midpoint::Vector{Float64} =
        param_set[is_growth_acc_mask .&& occursin.("midpoint", factor_names)].data

    # linear_extension scale factors with dimensions (cb_calib_groups ⋅ functional_groups)
    lin_ext_data = param_set[occursin.("linear_extension_scale", factor_names)].data
    total_calib_groups = div(length(lin_ext_data), n_groups)
    _linear_extension_scale_factors = reshape(
        lin_ext_data,
        total_calib_groups,
        n_groups
    )

    # mb_rate scale factors with dimensions (cb_calib_groups ⋅ functional_groups)
    mb_rate_data = param_set[occursin.("mb_rate_scale", factor_names)].data
    _mb_rate_scale_factors = reshape(
        mb_rate_data,
        total_calib_groups,
        n_groups
    )

    biogrp_lin_ext::Array{Float64,3} = repeat(_linear_extensions, 1, 1, total_calib_groups)
    biogrp_survival::Array{Float64,3} = repeat(survival_rate, 1, 1, total_calib_groups)
    for i in 1:total_calib_groups
        biogrp_lin_ext[:, :, i] .*= _linear_extension_scale_factors[i, :]
        biogrp_survival[:, :, i] .= apply_survival_scaling(
            biogrp_survival[:, :, i], _mb_rate_scale_factors[i, :]
        )
    end

    # Preallocate vector for growth constraints
    growth_constraints::Vector{Float64} = zeros(Float64, n_locs)

    habitable_loc_idxs = findall(habitable_locs)

    # Growth constraints and acceleration caches and masks
    relative_habitable_cover_cache = zeros(n_locs)
    growth_threshold_mask_cache::BitVector = trues(n_locs)
    cover_threshold_mask::BitVector = falses(n_locs)
    cache_habitable_max_projected_cover = copy(habitable_max_projected_cover)
    agg_cover_above_threshold_mask::BitVector = falses(n_habitable_locs)
    net_growth_rates = zeros(n_groups, n_sizes, n_locs)

    # Assume heat tolerance enhancement is based on population from `a_adapt_ref` years ago
    # If a_adapt == 0, use always first year as reference for a_adapt
    a_adapt_ref::Int64 = param_set[At("iv_CAq_a_adapt_ref")]

    # Depth attenuation of surface DHW (see `effective_dhw_at_depth`). Read once here rather
    # than inside the time loop - they are scenario-level constants.
    eff_dhw_base::Float64 = param_set[At("eff_dhw_base")]
    eff_dhw_mix::Float64 = param_set[At("eff_dhw_mix")]

    # c_mean tolerance of the "natural" population used as the reference for coral aquaculture corals
    # heat tolerance distribution
    c_mean_reference::Array{Float64,3} = if a_adapt_ref == 0
        # If a_adapt_ref == 0, 3rd dim holds c_mean at t-1 (idx 1), updated yearly, and
        # c_mean at t=1
        zeros(n_groups, n_locs, 2)
    else
        # If a_adapt_ref > 0, updates c_mean_reference yearly moving elements 1:end-1 to
        # 2:end and adds c_mean at t-1 (idx 1).
        zeros(n_groups, n_locs, a_adapt_ref)
    end

    # Keeps track of the heat tolerance means of juvenilles to use as a base for thermal
    # tolerance enhancement
    c_mean_reference .= c_mean_t[:, 1, :]

    cover_transition_lb = 0.05
    cover_transition_ub = 0.80

    # Preallocate per-timestep buffers to avoid repeated allocation in the hot loop
    Δcover_loss_proportion = zeros(n_locs)
    _is_reactive =
        is_reactive(param_set[At("iv_CAq_strategy")]) ||
        is_reactive(param_set[At("iv_Fog_strategy")]) ||
        is_reactive(param_set[At("iv_LvM_strategy")])
    dhw_p = is_guided ? similar(dhw_scen) : nothing
    current_loc_cover = zeros(n_locs)
    _loc_coral_cover = zeros(n_locs)
    leftover_space_m² = zeros(n_locs)
    _recruitment_col_sum = zeros(n_locs)
    _recruitment_col_sum_2d = reshape(_recruitment_col_sum, 1, n_locs)
    _permuted_buf = zeros(n_sizes, n_groups, n_locs)

    # --- COTS Submodel Initialization ---
    # Toggle COTS via ENV["ADRIA_COTS_ENABLED"]. Off by default so that runs without the
    # switch reproduce the coral-only model; set "true" to run the COTS submodel.
    cots_enabled = get(ENV, "ADRIA_COTS_ENABLED", "false") == "true"

    cots_prey_map = CotsPreyMap([1, 2, 3], [4, 5])  # Fast: Acropora spp; Slow: Massives

    # COTS parameters - read from the sampled param_set and converted to COTSMod runtime params
    cots_params = COTSMod.COTSParams(
        a = param_set[At("a")],
        b = param_set[At("b")],
        IMM = cots_enabled ? param_set[At("IMM")] : 0.0,
        p_tilde = param_set[At("p_tilde")],
        C_max = param_set[At("C_max")],
        m1 = param_set[At("m1")],
        m2 = param_set[At("m2")],
        m3 = param_set[At("m3")],
        a_F = param_set[At("a_F")],
        a_S = param_set[At("a_S")],
        h = param_set[At("h")],
        eta_F = param_set[At("eta_F")],
        eta_S = param_set[At("eta_S")],
        eta_starve = param_set[At("eta_starve")],
        eta_imm = param_set[At("eta_imm")],
        imm_threshold = param_set[At("imm_threshold")],
        fecundity_gate = param_set[At("fecundity_gate")] == 1,
        a_ricker = param_set[At("a_ricker")],
        b_ricker = param_set[At("b_ricker")],
        tau_condition = param_set[At("tau_condition")],
        allee_threshold = parse(
            Float64,
            get(ENV, "COTS_ALLEE_THRESHOLD", string(param_set[At("allee_threshold")]))
        ),
        separated_recruitment = lowercase(
            get(ENV, "COTS_SEPARATED_RECRUITMENT", "false")
        ) == "true",
        larval_fecundity = parse(
            Float64,
            get(ENV, "COTS_LARVAL_FECUNDITY", string(param_set[At("a_ricker")]))
        ),
        settlement_probability = parse(
            Float64, get(ENV, "COTS_SETTLEMENT_PROBABILITY", "1.0")
        ),
        juvenile_storage = lowercase(get(ENV, "COTS_JUVENILE_STORAGE", "false")) == "true",
        juvenile_maturation_min = parse(Float64, get(ENV, "COTS_JUVENILE_MATURATION_MIN", "0.05")),
        juvenile_maturation_max = parse(Float64, get(ENV, "COTS_JUVENILE_MATURATION_MAX", "1.0")),
        juvenile_maturation_cover = parse(Float64, get(ENV, "COTS_JUVENILE_MATURATION_COVER", "0.4")),
        habitat_mediation = lowercase(get(ENV, "COTS_HABITAT_MEDIATION", "false")) == "true",
        settlement_floor = parse(Float64, get(ENV, "COTS_SETTLEMENT_FLOOR", "0.05")),
        low_cover_mortality = lowercase(get(ENV, "COTS_LOW_COVER_MORTALITY", "false")) == "true",
        low_cover_threshold = parse(Float64, get(ENV, "COTS_LOW_COVER_THRESHOLD", "0.10")),
        low_cover_strength = parse(Float64, get(ENV, "COTS_LOW_COVER_STRENGTH", "0.0")),
        size_weighted_fecundity = lowercase(get(ENV, "COTS_SIZE_WEIGHTED_FECUNDITY", "false")) == "true",
        size_small_diameter_mm = parse(Float64, get(ENV, "COTS_SIZE_SMALL_DIAMETER_MM", "200.0")),
        size_large_diameter_mm = parse(Float64, get(ENV, "COTS_SIZE_LARGE_DIAMETER_MM", "300.0")),
        size_fecundity_slope_per_mm = parse(Float64, get(ENV, "COTS_SIZE_FECUNDITY_SLOPE_PER_MM", "0.0115")),
        size_initial_large_fraction = parse(Float64, get(ENV, "COTS_SIZE_INITIAL_LARGE_FRACTION", "0.5")),
        size_growth_max = parse(Float64, get(ENV, "COTS_SIZE_GROWTH_MAX", "0.5")),
        size_growth_cover = parse(Float64, get(ENV, "COTS_SIZE_GROWTH_COVER", "0.4")),
        size_growth_food_mediated = lowercase(get(ENV, "COTS_SIZE_GROWTH_FOOD_MEDIATED", "false")) == "true",
        lagged_adult_hazard = lowercase(get(ENV, "COTS_LAGGED_ADULT_HAZARD", "false")) == "true",
        hazard_food_coupled = lowercase(get(ENV, "COTS_HAZARD_FOOD_COUPLED", "false")) == "true",
        hazard_rho = parse(Float64, get(ENV, "COTS_HAZARD_RHO", "0.5")),
        hazard_threshold_ha = parse(Float64, get(ENV, "COTS_HAZARD_THRESHOLD_HA", "1.5")),
        hazard_width_ha = parse(Float64, get(ENV, "COTS_HAZARD_WIDTH_HA", "0.5")),
        hazard_coral_threshold = parse(Float64, get(ENV, "COTS_HAZARD_CORAL_THRESHOLD", "0.15")),
        hazard_coral_width = parse(Float64, get(ENV, "COTS_HAZARD_CORAL_WIDTH", "0.03")),
        hazard_max = parse(Float64, get(ENV, "COTS_HAZARD_MAX", "0.0")),
        preferred_prey_starvation = lowercase(get(ENV, "COTS_PREFERRED_PREY_STARVATION", "false")) == "true",
        starvation_preferred_threshold = parse(Float64, get(ENV, "COTS_STARVATION_PREFERRED_THRESHOLD", "0.05")),
        starvation_memory_years = parse(Float64, get(ENV, "COTS_STARVATION_MEMORY_YEARS", "0.0")),
        adult_senescence = lowercase(get(ENV, "COTS_ADULT_SENESCENCE", "false")) == "true",
        senescence_age = parse(Int, get(ENV, "COTS_SENESCENCE_AGE", "6")),
        senescence_mortality = parse(Float64, get(ENV, "COTS_SENESCENCE_MORTALITY", "0.8")),
        per_capita_consumption = lowercase(get(ENV, "COTS_PER_CAPITA_CONSUMPTION", "false")) == "true",
        consumption_m2_per_adult_year = parse(Float64, get(ENV, "COTS_CONSUMPTION_M2_PER_ADULT_YEAR", "10.0")),
        prey_switching = lowercase(get(ENV, "COTS_PREY_SWITCHING", "false")) == "true",
        switch_midpoint = parse(Float64, get(ENV, "COTS_SWITCH_MIDPOINT", "0.10")),
        switch_steepness = parse(Float64, get(ENV, "COTS_SWITCH_STEEPNESS", "70.0")),
        preferred_condition = lowercase(get(ENV, "COTS_PREFERRED_CONDITION", "false")) == "true",
        condition_saturation = parse(Float64, get(ENV, "COTS_CONDITION_SATURATION", "0.10")),
        condition_memory_years = parse(Float64, get(ENV, "COTS_CONDITION_MEMORY_YEARS", "1.0")),
        diet_linked_starvation = lowercase(get(ENV, "COTS_DIET_LINKED_STARVATION", "false")) == "true",
        nonpreferred_feeding_efficiency = parse(Float64,
            get(ENV, "COTS_NONPREFERRED_FEEDING_EFFICIENCY", "1.0")),
        juvenile_predation = lowercase(get(ENV, "COTS_JUVENILE_PREDATION", "false")) == "true",
        juvenile_predation_capacity = parse(Float64,
            get(ENV, "COTS_JUVENILE_PREDATION_CAPACITY", "20.0")),
        preferred_maturation = lowercase(get(ENV, "COTS_PREFERRED_MATURATION", "false")) == "true",
        starvation_bank = lowercase(get(ENV, "COTS_JUVENILE_STARVATION_BANK", "false")) == "true",
        allee_stretched = lowercase(get(ENV, "COTS_ALLEE_STRETCHED", "false")) == "true",
        allee_exponent = parse(Float64, get(ENV, "COTS_ALLEE_EXPONENT", "0.5")),
        starvation_post_grazing = lowercase(get(ENV, "COTS_STARVATION_POST_GRAZING", "false")) == "true",
        starvation_grace_years = parse(Int, get(ENV, "COTS_STARVATION_GRACE_YEARS", "0")),
        outbreak_collapse = lowercase(get(ENV, "COTS_OUTBREAK_COLLAPSE", "false")) == "true",
        collapse_threshold_ha = parse(Float64, get(ENV, "COTS_COLLAPSE_THRESHOLD_HA", "67.5")),
        collapse_duration_years = parse(Int, get(ENV, "COTS_COLLAPSE_DURATION_YEARS", "4")),
        collapse_background_ha = parse(Float64, get(ENV, "COTS_COLLAPSE_BACKGROUND_HA", "0.25")),
        collapse_juveniles = lowercase(get(ENV, "COTS_COLLAPSE_JUVENILES", "false")) == "true"
    )

    # Initialize COTS populations (zero if disabled)
    # TODO: Replace with empirical COTS distribution data when available
    # Parse debug overrides for COTS seeding (ENV → function args)
    _cots_seed_locs_str = get(ENV, "ADRIA_DEBUG_SEED_LOCATIONS", "")
    _cots_seed_n_str = get(ENV, "ADRIA_DEBUG_SEED_FIRST_N", "")
    _cots_init_density = parse(Float64, get(ENV, "ADRIA_DEBUG_INIT_DENSITY", "0.1"))
    _cots_pulse_start = parse(Int, get(ENV, "COTS_PULSE_START", "1"))
    _cots_pulse_duration = max(1, parse(Int, get(ENV, "COTS_PULSE_DURATION", "1")))
    _cots_pulse_repeat_interval = parse(Int, get(ENV, "COTS_PULSE_REPEAT_INTERVAL", "0"))
    _cots_pulse_relative_magnitude = parse(Float64, get(ENV, "COTS_PULSE_RELATIVE_MAGNITUDE", "0.0"))
    _cots_seed_locs = if _cots_seed_locs_str != ""
        Set(parse.(Int, split(_cots_seed_locs_str, ",")))
    elseif _cots_seed_n_str != ""
        Set(1:min(n_locs, parse(Int, _cots_seed_n_str)))
    else
        nothing  # let init_cots_populations use its default random seeding
    end

    _cots_spatial_initial_density = if hasproperty(domain, :cots_init_density) && !isnothing(domain.cots_init_density)
        domain.cots_init_density
    else
        nothing
    end

    cots_models = initialize_cots(
        n_locs,
        cots_params;
        enabled=cots_enabled,
        spatial_initial_density=_cots_spatial_initial_density,
        seed_locs=_cots_seed_locs,
        init_density=_cots_init_density,
        rng=rng
    )
    # Explicit historical stage densities are an opt-in initial-condition
    # treatment. The legacy probability-based initializer remains the default.
    _cots_initial_state_csv = get(ENV, "COTS_INITIAL_STATE_CSV", "")
    if !isempty(_cots_initial_state_csv)
        cots_enabled || error("COTS_INITIAL_STATE_CSV requires COTS to be enabled")
        stage_state = load_cots_initial_state_csv(_cots_initial_state_csv, domain.loc_ids)
        apply_cots_initial_state!(cots_models, stage_state)
    end

    # Use COTS-specific larval connectivity when the domain provides it, while retaining
    # coral connectivity as the backwards-compatible fallback for existing domains.
    cots_conn = sparse(cots_connectivity(domain).data)
    cots_max_larval_supply = zeros(Float64, n_locs)
    cots_temporal_mode = lowercase(get(ENV, "COTS_CONNECTIVITY_TEMPORAL_MODE", "mean"))
    cots_temporal_mode in ["mean", "cycle", "sample"] || error(
        "COTS_CONNECTIVITY_TEMPORAL_MODE must be mean, cycle, or sample"
    )
    cots_forcing = hasproperty(domain, :cots_forcing) ? domain.cots_forcing : nothing
    cots_temporal_mode == "mean" || !isnothing(cots_forcing) || error(
        "Temporal COTS connectivity requested for a domain without annual forcing"
    )
    cots_connectivity_seed = parse(Int, get(ENV, "COTS_CONNECTIVITY_SEED", "20260929"))
    cots_apply_larval_survival = lowercase(get(ENV, "COTS_APPLY_LARVAL_SURVIVAL", "false")) == "true"
    cots_separated_recruitment = cots_params.separated_recruitment
    blackout_start = get(ENV, "COTS_INTERNAL_SETTLEMENT_BLACKOUT_START_YEAR", "")
    blackout_end = get(ENV, "COTS_INTERNAL_SETTLEMENT_BLACKOUT_END_YEAR", "")
    isempty(blackout_start) == isempty(blackout_end) || error(
        "Both COTS internal-settlement blackout years must be set together"
    )
    cots_internal_blackout = if isempty(blackout_start)
        nothing
    else
        start_year = parse(Int, blackout_start)
        end_year = parse(Int, blackout_end)
        start_year <= end_year || error("COTS internal-settlement blackout years are reversed")
        cots_separated_recruitment || error(
            "COTS internal-settlement blackout requires separated recruitment"
        )
        start_year:end_year
    end
    grazing_off_start = get(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_START_YEAR", "")
    grazing_off_end = get(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_END_YEAR", "")
    isempty(grazing_off_start) == isempty(grazing_off_end) || error(
        "Both COTS diagnostic grazing-off years must be set together"
    )
    cots_grazing_off_years = if isempty(grazing_off_start)
        nothing
    else
        start_year = parse(Int, grazing_off_start)
        end_year = parse(Int, grazing_off_end)
        start_year <= end_year || error("COTS diagnostic grazing-off years are reversed")
        start_year:end_year
    end
    cots_params.larval_fecundity >= 0.0 ||
        error("COTS_LARVAL_FECUNDITY must be non-negative")
    0.0 <= cots_params.settlement_probability <= 1.0 ||
        error("COTS_SETTLEMENT_PROBABILITY must lie in [0, 1]")
    cots_params.low_cover_threshold > 0.0 ||
        error("COTS_LOW_COVER_THRESHOLD must be positive")
    0.0 <= cots_params.low_cover_strength <= 1.0 ||
        error("COTS_LOW_COVER_STRENGTH must lie in [0, 1]")
    if cots_params.lagged_adult_hazard
        isfinite(cots_params.hazard_rho) && 0.0 <= cots_params.hazard_rho < 1.0 ||
            error("COTS_HAZARD_RHO must lie in [0, 1)")
        isfinite(cots_params.hazard_threshold_ha) && cots_params.hazard_threshold_ha >= 0.0 ||
            error("COTS_HAZARD_THRESHOLD_HA must be nonnegative")
        isfinite(cots_params.hazard_width_ha) && cots_params.hazard_width_ha > 0.0 ||
            error("COTS_HAZARD_WIDTH_HA must be positive")
        isfinite(cots_params.hazard_coral_threshold) &&
            0.0 <= cots_params.hazard_coral_threshold <= 1.0 ||
            error("COTS_HAZARD_CORAL_THRESHOLD must lie in [0, 1]")
        isfinite(cots_params.hazard_coral_width) && cots_params.hazard_coral_width > 0.0 ||
            error("COTS_HAZARD_CORAL_WIDTH must be positive")
        isfinite(cots_params.hazard_max) && cots_params.hazard_max >= 0.0 ||
            error("COTS_HAZARD_MAX must be nonnegative")
    end
    if cots_params.preferred_prey_starvation
        0.0 < cots_params.starvation_preferred_threshold <= 1.0 ||
            error("COTS_STARVATION_PREFERRED_THRESHOLD must lie in (0, 1]")
        cots_params.starvation_memory_years >= 0.0 ||
            error("COTS_STARVATION_MEMORY_YEARS must be nonnegative")
    end
    if cots_params.adult_senescence
        cots_params.senescence_age >= 3 || error("COTS_SENESCENCE_AGE must be at least 3")
        0.0 < cots_params.senescence_mortality <= 1.0 ||
            error("COTS_SENESCENCE_MORTALITY must lie in (0, 1]")
    end
    cots_params.consumption_m2_per_adult_year >= 0.0 ||
        error("COTS_CONSUMPTION_M2_PER_ADULT_YEAR must be nonnegative")
    if cots_params.prey_switching
        cots_params.per_capita_consumption ||
            error("COTS_PREY_SWITCHING requires COTS_PER_CAPITA_CONSUMPTION=true")
        0.0 <= cots_params.switch_midpoint <= 1.0 ||
            error("COTS_SWITCH_MIDPOINT must lie in [0, 1]")
        cots_params.switch_steepness > 0.0 || error("COTS_SWITCH_STEEPNESS must be positive")
    end
    if cots_params.preferred_condition
        0.0 < cots_params.condition_saturation <= 1.0 ||
            error("COTS_CONDITION_SATURATION must lie in (0, 1]")
        cots_params.condition_memory_years > 0.0 ||
            error("COTS_CONDITION_MEMORY_YEARS must be positive")
    end
    if cots_params.diet_linked_starvation
        cots_params.prey_switching && cots_params.preferred_prey_starvation || error(
            "COTS_DIET_LINKED_STARVATION requires COTS_PREY_SWITCHING and COTS_PREFERRED_PREY_STARVATION")
    end
    0.0 <= cots_params.nonpreferred_feeding_efficiency <= 1.0 ||
        error("COTS_NONPREFERRED_FEEDING_EFFICIENCY must lie in [0, 1]")
    cots_params.juvenile_predation && !(cots_params.juvenile_predation_capacity > 0.0) &&
        error("COTS_JUVENILE_PREDATION_CAPACITY must be positive")
    cots_params.preferred_maturation && cots_params.juvenile_storage &&
        error("COTS_PREFERRED_MATURATION and COTS_JUVENILE_STORAGE are exclusive")
    cots_params.allee_stretched && !(cots_params.allee_exponent > 0.0) &&
        error("COTS_ALLEE_EXPONENT must be positive")
    cots_params.starvation_grace_years >= 0 ||
        error("COTS_STARVATION_GRACE_YEARS must be nonnegative")
    cots_params.starvation_grace_years > 0 && !(cots_params.preferred_prey_starvation &&
        !cots_params.diet_linked_starvation) &&
        error("COTS_STARVATION_GRACE_YEARS requires COTS_PREFERRED_PREY_STARVATION without " *
              "COTS_DIET_LINKED_STARVATION")
    cots_params.starvation_post_grazing && !(cots_params.per_capita_consumption &&
        cots_params.preferred_prey_starvation && !cots_params.prey_switching) &&
        error("COTS_STARVATION_POST_GRAZING requires COTS_PER_CAPITA_CONSUMPTION and " *
              "COTS_PREFERRED_PREY_STARVATION without COTS_PREY_SWITCHING")
    if cots_params.outbreak_collapse
        isfinite(cots_params.collapse_threshold_ha) && cots_params.collapse_threshold_ha > 0.0 ||
            error("COTS_COLLAPSE_THRESHOLD_HA must be positive")
        cots_params.collapse_duration_years >= 1 ||
            error("COTS_COLLAPSE_DURATION_YEARS must be at least 1")
        isfinite(cots_params.collapse_background_ha) &&
            0.0 <= cots_params.collapse_background_ha < cots_params.collapse_threshold_ha ||
            error("COTS_COLLAPSE_BACKGROUND_HA must lie in [0, COTS_COLLAPSE_THRESHOLD_HA)")
    end
    # Per-capita grazing converts adults per site hectare into habitable-cover
    # fractions; the legacy response is already a fraction and ignores this.
    cots_habitable_fraction = cots_params.per_capita_consumption ?
        Float64.(location_k(domain)) : nothing
    cots_separated_recruitment && isnothing(cots_forcing) && error(
        "Separated COTS recruitment requires reef-scale forcing and reef areas"
    )
    cots_mean_reef_conn = isnothing(cots_forcing) ? nothing :
        reduce(+, cots_forcing.reef_connectivity) ./ length(cots_forcing.reef_connectivity)
    cots_mean_source_survival = isnothing(cots_forcing) ? nothing :
        vec(sum(cots_forcing.source_survival; dims=2)) ./ size(cots_forcing.source_survival, 2)
    # Kept as COTS_EXTERNAL_SOURCE_DENSITY for compatibility, but this is an
    # absolute pre-dispersal recruit-production flux, not external adult density.
    cots_external_source_recruits = parse(Float64, get(ENV, "COTS_EXTERNAL_SOURCE_DENSITY", "0.0"))
    cots_external_source_recruits >= 0.0 || error("COTS_EXTERNAL_SOURCE_DENSITY must be non-negative")
    cots_separated_recruitment && cots_external_source_recruits > 0.0 && error(
        "COTS_EXTERNAL_SOURCE_DENSITY is a legacy settled-supply control; use " *
        "COTS_EXTERNAL_OUTBREAK_PRODUCTION with separated recruitment"
    )
    cots_external_outbreak_production = parse(
        Float64, get(ENV, "COTS_EXTERNAL_OUTBREAK_PRODUCTION", "0.0")
    )
    cots_external_outbreak_production >= 0.0 ||
        error("COTS_EXTERNAL_OUTBREAK_PRODUCTION must be non-negative")
    # Labelled counterfactual: external outbreak supply multiplier for a calendar-year window.
    cots_supply_window = parse_cots_supply_window(get(ENV, "COTS_EXTERNAL_SUPPLY_WINDOW", ""))
    cots_external_outbreak_production > 0.0 && !cots_separated_recruitment && error(
        "Evidence-derived outbreak supply requires COTS_SEPARATED_RECRUITMENT=true"
    )
    cots_calendar_start_year = parse(
        Int, get(ENV, "COTS_CALENDAR_START_YEAR", "0")
    )
    raw_cots_timeframe = Int.(domain.env_layer_md.timeframe[1:tf])
    cots_calendar_years = if all(raw_cots_timeframe .>= 1000)
        raw_cots_timeframe
    elseif cots_calendar_start_year > 0
        collect(cots_calendar_start_year:(cots_calendar_start_year + tf - 1))
    else
        raw_cots_timeframe
    end
    # Opt-in forced grazing: prescribed adult density (e.g. observed CPUE / q) grazes
    # coral with the per-capita preferred-first response; COTS dynamics must be off.
    cots_forced_adults_csv = get(ENV, "COTS_FORCED_ADULTS_CSV", "")
    cots_forced_adults = if isempty(cots_forced_adults_csv)
        nothing
    else
        cots_enabled && error("COTS_FORCED_ADULTS_CSV requires ADRIA_COTS_ENABLED=false")
        load_cots_forced_adults(cots_forced_adults_csv, domain.loc_ids, cots_calendar_years)
    end
    cots_forced_habitable = isnothing(cots_forced_adults) ? nothing :
        Float64.(location_k(domain))
    # Opt-in observed culling: pre-spawning removals before the COTS step, post-spawning
    # removals after it (before dispersal). Unset = established path.
    # Opt-in larval-survival index (multiplier on each site's larval production at source,
    # applied before dispersal). Unset = established path.
    cots_larval_survival_csv = get(ENV, "COTS_LARVAL_SURVIVAL_CSV", "")
    cots_larval_survival = if isempty(cots_larval_survival_csv)
        nothing
    else
        cots_enabled || error("COTS_LARVAL_SURVIVAL_CSV requires ADRIA_COTS_ENABLED=true")
        cots_separated_recruitment || error("COTS_LARVAL_SURVIVAL_CSV requires COTS_SEPARATED_RECRUITMENT=true")
        load_cots_larval_survival(cots_larval_survival_csv, domain.loc_ids, cots_calendar_years)
    end
    # Opt-in reef-habitat settlement modifier (time-invariant multiplier on each site's
    # settlement from all sources, applied after the transition). Unset = established path.
    cots_settlement_modifier_csv = get(ENV, "COTS_SETTLEMENT_MODIFIER_CSV", "")
    cots_settlement_modifier = if isempty(cots_settlement_modifier_csv)
        nothing
    else
        cots_enabled || error("COTS_SETTLEMENT_MODIFIER_CSV requires ADRIA_COTS_ENABLED=true")
        cots_separated_recruitment || error("COTS_SETTLEMENT_MODIFIER_CSV requires COTS_SEPARATED_RECRUITMENT=true")
        load_cots_settlement_modifier(cots_settlement_modifier_csv, domain.loc_ids)
    end
    # Opt-in year-specific survival multiplier on external pelagic supply (regional chl proxy).
    # Unset = established path.
    cots_external_survival_csv = get(ENV, "COTS_EXTERNAL_SURVIVAL_CSV", "")
    cots_external_survival = isempty(cots_external_survival_csv) ? nothing :
        load_cots_external_survival(cots_external_survival_csv, cots_calendar_years)
    cots_cull_csv = get(ENV, "COTS_CULL_REMOVALS_CSV", "")
    cots_cull = if isempty(cots_cull_csv)
        nothing
    else
        cots_enabled || error("COTS_CULL_REMOVALS_CSV requires ADRIA_COTS_ENABLED=true")
        load_cots_cull_removals(cots_cull_csv, domain.loc_ids, cots_calendar_years)
    end
    # Cull allocation: "site" (default) applies the table's site requests; "reef_proportional"
    # pools each reef's requested count and allocates it to the reef's sites by COTS present.
    cots_cull_allocation = get(ENV, "COTS_CULL_ALLOCATION", "site")
    cots_cull_allocation in ("site", "reef_proportional") ||
        error("COTS_CULL_ALLOCATION must be site or reef_proportional")
    cots_cull_groups, cots_cull_areas = if cots_cull_allocation == "reef_proportional"
        isnothing(cots_cull) && error("COTS_CULL_ALLOCATION=reef_proportional requires COTS_CULL_REMOVALS_CSV")
        reef_ids = string.(domain.loc_data.GBRMPA_ID)
        [findall(==(r), reef_ids) for r in unique(reef_ids)], Float64.(domain.loc_data.area)
    else
        nothing, nothing
    end
    # Projection switches (crash_mechanism_protocol.md, 2026-10-21). Unset = established path.
    # M1 forward effort culling: reef diver-hours x multiplier ("fixed") or the domain total
    # reallocated to outbreak reefs ("targeted"), from COTS_CULL_EFFORT_START_YEAR.
    cots_cull_effort_csv = get(ENV, "COTS_CULL_EFFORT_CSV", "")
    cots_reset_csv = get(ENV, "COTS_STATE_RESET_CSV", "")
    # Reef groups (GBRMPA reef id) and site areas, read only when a switch needs them.
    cots_reef_labels, cots_reef_groups, cots_site_areas =
        isempty(cots_cull_effort_csv) && isempty(cots_reset_csv) ? (nothing, nothing, nothing) : let
            ids = string.(domain.loc_data.GBRMPA_ID)
            labels = unique(ids)
            # loc_data.area is m^2; effort culling needs ha (COTS ha^-1 x ha = COTS).
            labels, [findall(==(r), ids) for r in labels], Float64.(domain.loc_data.area) ./ 1e4
        end
    cots_cull_effort = if isempty(cots_cull_effort_csv)
        nothing
    else
        cots_enabled || error("COTS_CULL_EFFORT_CSV requires ADRIA_COTS_ENABLED=true")
        load_cots_cull_effort(cots_cull_effort_csv, cots_reef_labels)
    end
    cots_effort_start = parse(Int, get(ENV, "COTS_CULL_EFFORT_START_YEAR", "0"))
    # Optional last year of forward effort culling (2026-10-22); unset = no end.
    cots_effort_end = parse(Int, get(ENV, "COTS_CULL_EFFORT_END_YEAR", "9999"))
    cots_effort_end >= cots_effort_start || error("COTS_CULL_EFFORT_END_YEAR precedes the start year")
    cots_effort_hours_log = zeros(Float64, tf)  # diver-hours deployed per year (both phases)
    cots_effort_multiplier = parse(Float64, get(ENV, "COTS_CULL_EFFORT_MULTIPLIER", "1.0"))
    cots_effort_mode = get(ENV, "COTS_CULL_EFFORT_MODE", "fixed")
    cots_effort_mode in ("fixed", "targeted") || error("COTS_CULL_EFFORT_MODE must be fixed or targeted")
    cots_effort_threshold = parse(Float64, get(ENV, "COTS_CULL_EFFORT_THRESHOLD_HA", "3.0"))
    cots_effort_cpue_max = parse(Float64, get(ENV, "COTS_CULL_CPUE_MAX", "20.0"))
    isnothing(cots_cull_effort) || (cots_effort_start > 0 && cots_effort_multiplier >= 0.0) ||
        error("COTS_CULL_EFFORT_CSV needs COTS_CULL_EFFORT_START_YEAR > 0 and a non-negative multiplier")
    # M2 projection-start reset (data assimilation) after the table's year.
    cots_reset = if isempty(cots_reset_csv)
        nothing
    else
        cots_enabled || error("COTS_STATE_RESET_CSV requires ADRIA_COTS_ENABLED=true")
        load_cots_state_reset(cots_reset_csv, domain.loc_ids)
    end
    cots_reset_juveniles = get(ENV, "COTS_STATE_RESET_JUVENILES", "true") == "true"
    cots_reset_others = Symbol(get(ENV, "COTS_STATE_RESET_OTHERS", "median"))
    COTS_RESET_LOG[] = nothing
    # M3 coral-coupled outside supply: factor min(1, domain coral / its pre-start mean).
    cots_coral_coupling_start = parse(Int, get(ENV, "COTS_EXTERNAL_CORAL_COUPLING_START_YEAR", "0"))
    cots_coral_weights = cots_coral_coupling_start > 0 ? Float64.(location_k(domain)) .* Float64.(domain.loc_data.area) : nothing
    cots_coral_reference = Float64[]  # domain coral before the start year
    cots_coral_factor_log = fill(NaN, tf)
    cots_external_outbreak_production > 0.0 &&
        !any(in(cots_forcing.outbreak_years), cots_calendar_years) && error(
            "Evidence-derived outbreak supply has no calendar-year overlap; " *
            "set COTS_CALENDAR_START_YEAR for relative-time domains"
        )

    # COTS population log: [timesteps, 3 age classes, locations]
    Ycots = zeros(tf, 3, n_locs)
    # COTS body condition log: [timesteps, locations]
    Ycots_bc = zeros(tf, n_locs)
    # Flow channels 1:7 preserve the legacy output contract. Channels 8:12 add
    # local retention, pelagic survivors, settled recruits, potential external
    # production, and pelagic external supply.
    Ycots_flow = zeros(tf, 12, n_locs)
    # Opt-in adult-size diagnostics: small, large, effective breeders, growth.
    Ycots_size = zeros(tf, 4, n_locs)
    # Burden (COTS ha^-1), one-year cumulative hazard, extra adult deaths (COTS ha^-1 yr^-1).
    Ycots_hazard = zeros(tf, 3, n_locs)
    # Remembered preferred cover (habitable fraction; NaN when off), realized food
    # survival factor, senescent adults and excess senescent deaths (COTS ha^-1),
    # the share of grazing demand aimed at preferred prey, juvenile predation removals
    # (COTS ha^-1 yr^-1) and the juvenile maturation fraction.
    Ycots_mechanism = zeros(tf, 10, n_locs)
    Ycots_forcing_year = zeros(Int, tf)
    # Observed cull removals actually applied (COTS ha^-1): pre juveniles, pre adults,
    # post juveniles, post adults. Zero unless COTS_CULL_REMOVALS_CSV is set.
    Ycots_culled = zeros(Float64, tf, 4, n_locs)
    # Larval-survival multipliers applied (1 unless COTS_LARVAL_SURVIVAL_CSV is set).
    Ycots_larval_survival = ones(Float64, tf, n_locs)

    for tstep::Int64 in 2:tf
        # Convert cover to absolute values to use within CoralBlox model
        C_cover_t[:, :, habitable_locs] .=
            C_cover[tstep - 1, :, :, habitable_locs] .* habitable_loc_areas′

        # To prevent overgrowth when adding recruitment to C_cover_t, set a threshold above
        # which recruits are either set to zero or rescaled.
        constrain_recruitment!(
            recruitment, agg_cover_above_threshold_mask, C_cover_t, habitable_loc_areas
        )

        # Settlers from t-1 grow into observable sizes.
        C_cover_t[:, 1, habitable_locs] .+= recruitment

        sum!(_recruitment_col_sum_2d, recruitment)
        habitable_max_projected_cover .=
            cache_habitable_max_projected_cover .+ _recruitment_col_sum

        relative_habitable_cover_cache[habitable_locs] .=
            loc_coral_cover(C_cover_t[:, :, habitable_locs]) ./
            vec_abs_k[habitable_locs]

        # Growth constrains need to be calculated seperately for differen growth rates
        growth_threshold_mask_cache .=
            relative_habitable_cover_cache .>= cover_transition_ub
        for idx in 1:total_calib_groups
            cover_threshold_mask .=
                cb_calib_group_masks[:, idx] .&& growth_threshold_mask_cache

            # Only apply linear_extension_scale_factors to locations with high cover
            @views growth_constraints[cover_threshold_mask] .= linear_extension_scale_factors(
                C_cover_t[:, :, cover_threshold_mask],
                vec_abs_k[cover_threshold_mask],
                biogrp_lin_ext[:, :, idx],
                _bin_edges,
                habitable_max_projected_cover[cover_threshold_mask]
            )
        end

        # When cover is between cover_transition_lb and cover_transition_ub use a weighted mean
        @. growth_threshold_mask_cache = (
            cover_transition_lb <= relative_habitable_cover_cache <
            cover_transition_ub
        )
        if sum(growth_threshold_mask_cache) > 0
            for idx in 1:total_calib_groups
                cover_threshold_mask .=
                    growth_threshold_mask_cache .&& cb_calib_group_masks[:, idx]
                transition_scale =
                    (
                        relative_habitable_cover_cache[cover_threshold_mask] .-
                        cover_transition_lb
                    ) ./ (cover_transition_ub - cover_transition_lb)
                growth_constraints[cover_threshold_mask] .=
                    (
                        transition_scale .*
                        linear_extension_scale_factors(
                            C_cover_t[:, :, cover_threshold_mask],
                            vec_abs_k[cover_threshold_mask],
                            biogrp_lin_ext[:, :, idx],
                            _bin_edges,
                            habitable_max_projected_cover[cover_threshold_mask]
                        )
                    ) .+
                    (
                        (1 .- transition_scale) .*
                        growth_acceleration.(
                            growth_acc_height[idx],
                            growth_acc_midpoint[idx],
                            growth_acc_steepness[idx],
                            relative_habitable_cover_cache[cover_threshold_mask]
                        )
                    )
            end
        end
        growth_threshold_mask_cache .=
            relative_habitable_cover_cache .< cover_transition_lb
        for idx in 1:total_calib_groups
            cover_threshold_mask .=
                growth_threshold_mask_cache .&& cb_calib_group_masks[:, idx]
            growth_constraints[cover_threshold_mask] .= growth_acceleration.(
                growth_acc_height[idx],
                growth_acc_midpoint[idx],
                growth_acc_steepness[idx],
                relative_habitable_cover_cache[cover_threshold_mask]
            )
        end

        @inbounds for i in habitable_loc_idxs
            cb_idx::Int64 = loc_cb_calib_group_idxs[i]
            net_growth_rates[:, :, i] .=
                @view(biogrp_lin_ext[:, :, cb_idx]) .*
                growth_constraints[i]

            timestep!(
                functional_groups[i],
                @view(recruitment[:, i]),
                @view(net_growth_rates[:, :, i]),
                @view(biogrp_survival[:, :, cb_idx])
            )

            # Write to the cover matrix
            coral_cover(functional_groups[i], @view(C_cover_t[:, :, i]))
        end

        # Check if size classes are inappropriately out-growing habitable area
        if any(loc_coral_cover(C_cover_t)[habitable_locs] .> habitable_loc_areas)
            @warn "Cover outgrowing habitable area at tstep $tstep. Constraining."
            outgrowing_locs_mask =
                loc_coral_cover(C_cover_t)[habitable_locs] .> habitable_loc_areas
            C_cover_t[:, :, outgrowing_locs_mask] .*=
                reshape(
                    vec_abs_k[outgrowing_locs_mask] ./
                    loc_coral_cover(C_cover_t)[outgrowing_locs_mask],
                    (1, 1, count(outgrowing_locs_mask))
                ) .* 0.999
        end

        # Convert C_cover_t to relative values after CoralBlox was run
        C_cover_t[:, :, habitable_locs] .= (
            @view(C_cover_t[:, :, habitable_locs]) ./ habitable_loc_areas′
        )
        C_cover[tstep, :, :, habitable_locs] .= @view(C_cover_t[:, :, habitable_locs])

        # Natural adaptation (doesn't change C_cover_t)
        if tstep <= tf
            adjust_DHW_distribution!(
                @view(C_cover[tstep - 1, :, :, :]), c_mean_t, net_growth_rates,
                c_mean_tol_ceil
            )

            if log_dhw_tols
                # Log dhw tolerances if requested
                permutedims!(_permuted_buf, c_mean_t, (2, 1, 3))
                dhw_tol_mean_log[tstep, :, :] .= reshape(
                    _permuted_buf, size(dhw_tol_mean_log)[2:3]
                )
            end

            if log_cover
                permutedims!(_permuted_buf, @view(C_cover[tstep, :, :, :]), (2, 1, 3))
                cover_log[tstep, :, :] .= reshape(_permuted_buf, size(cover_log)[2:3])
            end
        end

        # M4 (opt-in): annual rise in mean DHW tolerance, under the HEAT_UB ceiling.
        if thermal_acclimation_rate > 0.0 && coral_calendar_years[tstep] >= thermal_acclimation_start
            c_mean_t .= min.(c_mean_t .+ thermal_acclimation_rate, c_mean_tol_ceil)
        end

        # Snapshot juvenile tolerance after natural adaptation, before coral aquaculture modifies it.
        # This is the correct baseline for a_adapt enhancement: using the post-coral-aquaculture state
        # here would cause the enhancement to compound each deployment year.
        c_mean_pre_CAq .= c_mean_t[:, 1, :]

        # Reproduction
        # Calculates scope for coral fedundity for each size class and at each location
        fecundity_scope!(
            fec_scope, fecundity_per_m², C_cover_t, habitable_areas;
            apply_allee_effect=apply_allee_effect
        )

        for l = 1:n_locs
            s = sum(fec_scope[:, l])
            if s > 0.0
                prop_fecundity[:, l] .= fec_scope[:, l] ./ s
            else
                prop_fecundity[:, l] .= 0.0
            end
        end

        _loc_coral_cover .= loc_coral_cover(C_cover_t)
        leftover_space_m² .= relative_leftover_space(_loc_coral_cover) .* vec_abs_k

        # Reset potential settlers to zero
        potential_settlers .= 0.0
        recruitment .= 0.0

        # habitable_areas_view = @view habitable_areas[:, habitable_locs]

        # M5 (opt-in): outside larvae from the start year at `fraction` x the earlier domain mean.
        ext_larvae_t = nothing
        if !isnothing(ext_coral_larvae) && coral_calendar_years[tstep] >= ext_coral_start
            ext_coral_n[] > 0 || error("ADRIA_EXTERNAL_CORAL_LARVAE_START_YEAR leaves no reference years")
            ext_coral_larvae[:, habitable_locs] .= ext_coral_fraction .* ext_coral_sum ./ ext_coral_n[]
            ext_larvae_t = ext_coral_larvae
        end

        # Recruitment represents additional cover, relative to total location area
        # Recruitment/settlement occurs after the full moon in October/November
        @views recruitment[:, habitable_loc_idxs] .=
            settler_cover(
                fec_scope,
                conn,
                leftover_space_m²,
                sim_params.max_settler_density,
                sim_params.max_larval_density,
                basal_area_per_settler,
                potential_settlers;
                rng=rng,
                external_larvae=ext_larvae_t
            )[
                :, habitable_loc_idxs
            ] ./ habitable_areas[:, habitable_loc_idxs]
        if !isnothing(ext_coral_larvae) && coral_calendar_years[tstep] < ext_coral_start
            # Reference: habitable-area-weighted domain-mean larval density per group.
            ext_coral_sum .+= vec(sum(@view(potential_settlers[:, habitable_locs]) .*
                habitable_loc_areas'; dims=2)) ./ sum(habitable_loc_areas)
            ext_coral_n[] += 1
        end

        # Reset fecundity scope before next run
        fec_scope .= 0.0

        sum!(_recruitment_col_sum_2d, recruitment)
        leftover_space_m² .-= _recruitment_col_sum .* vec_abs_k

        # Determine intervention locations whose deployment is assumed to occur
        # between November to February.
        #
        # Nominal sequence of events is conceptualised as:
        # - Corals spawn and are recruited
        # - Shading is applied
        # - Fogging is applied next
        # - Coral aquaculture interventions occur
        # - bleaching then occurs
        # - then cyclones hit
        #
        # Coral aquaculture occurs before bleaching as current tech only allows deployments shortly
        # after spawning. If bio-banking or similar tech comes up to speed then we could
        # potentially deploy at alternate times.

        # Shading
        # Apply regional cooling effect before selecting locations for coral aquaculture
        dhw_t .= dhw_scen[tstep, :]  # subset of DHW for given timestep
        if apply_Shd && Shd_decision_years[tstep]
            YShd[tstep, :] .= Shd

            # Apply reduction in DHW due to shading
            dhw_t[Shd_locs_mask] .= max.(0.0, dhw_t[Shd_locs_mask] .- Shd)
        end

        if is_guided
            # Use modified projected DHW (may have been affected by shading)
            dhw_p .= dhw_scen
            dhw_p[tstep, :] .= dhw_t
            dhw_projection .= weighted_projection(dhw_p, tstep, plan_horizon, decay, tf)

            if wave_idx > 0.0
                wave_projection .= weighted_projection(
                    wave_scen, tstep, plan_horizon, decay, tf
                )
            end

            # Determine connectivity strength weighting by area.
            # Accounts for strength of connectivity where there is low/no coral cover
            in_conn, out_conn, _ = connectivity_strength(
                area_weighted_conn, vec(_loc_coral_cover), conn_cache
            )

            # Calculate current location cover
            current_loc_cover .= dropdims(sum(C_cover_t; dims=(1, 2)); dims=(1, 2))
        end

        # Fogging
        # if is_guided
        #     if fog_decision_years[tstep] && (fogging .> 0.0)
        #         selected_Fog_ranks = select_locations(
        #             Fog_pref,
        #             decision_mat,
        #             MCDA_approach,
        #             min_iv_locs
        #         )

        #         if !isempty(selected_Fog_ranks)
        #             log_location_ranks[tstep, At(selected_Fog_ranks), At(:fog)] .=
        #                 1:length(selected_Fog_ranks)
        #         end
        #     end
        # elseif apply_fogging && fog_decision_years[tstep]
        #     selected_Fog_ranks = unguided_selection(
        #         domain.loc_ids,
        #         min_iv_locs,
        #         vec(leftover_space_m²)
        #         # depth_criteria
        #     )

        #     log_location_ranks[tstep, At(selected_Fog_ranks), At(:fog)] .= 1.0
        # end

        # Fog location selection
        if !isnothing(Fog_strategy)
            state = if is_guided
                build_state(
                    domain, Fog_strategy,
                    (
                        current_cover=current_loc_cover,
                        recent_cover_losses=recent_cover_losses,
                        last_deployment=iv_CAq_last_deployment
                    )
                )
            else
                nothing
            end

            # Get candidate locations from strategy
            candidate_locs = filter_candidate_locations(Fog_strategy, tstep, state)
            candidate_loc_indices = findall(
                in.(domain.loc_ids, Ref(candidate_locs))
            )

            if !isempty(candidate_locs)
                selected_Fog_ranks = []
                if is_guided
                    # Update decision matrix with current conditions
                    update_criteria_values!(
                        Fog_decision_mat[location = At(candidate_locs)];
                        heat_stress=dhw_projection[candidate_loc_indices],
                        wave_stress=wave_projection[candidate_loc_indices],
                        coral_cover=current_loc_cover[candidate_loc_indices],
                        in_connectivity=in_conn[candidate_loc_indices],
                        out_connectivity=out_conn[candidate_loc_indices]
                    )

                    # Build state for target locations only
                    selected_Fog_ranks = select_locations(
                        Fog_pref,
                        Fog_decision_mat[location = At(candidate_locs)],
                        MCDA_approach,
                        min_iv_locs
                    )
                else
                    selected_Fog_ranks = unguided_selection(
                        candidate_locs,
                        min_iv_locs,
                        vec(leftover_space_m²[candidate_loc_indices]);
                        rng=rng
                    )
                end
                if !isempty(selected_Fog_ranks)
                    log_val = is_guided ? (1:length(selected_Fog_ranks)) : 1.0
                    log_location_ranks[tstep, At(selected_Fog_ranks), At(:fog)] .=
                        log_val
                    selected_Fog_loc_idx = findall(
                        in.(domain.loc_ids, Ref(selected_Fog_ranks))
                    )
                    last_Fog_deployment[selected_Fog_loc_idx] .= tstep
                end
            end
        end

        has_Fog_locs::Bool = !isempty(selected_Fog_ranks)

        # Fog selected locations
        if has_Fog_locs  # fog_decision_years[tstep] &&
            Fog_locs = findall(log_location_ranks.locations .∈ [selected_Fog_ranks])
            Fog_locations!(@view(Yfog[tstep, :]), Fog_locs, dhw_t, fogging)

            # Empty selected_Fog_ranks before the next iteration
            selected_Fog_ranks = []
        end

        # Larval methods intervention
        if !isnothing(LvM_strategy)
            state = if is_guided
                build_state(
                    domain, LvM_strategy,
                    (
                        current_cover=current_loc_cover,
                        recent_cover_losses=recent_cover_losses,
                        last_deployment=last_LvM_deployment
                    )
                )
            else
                nothing
            end

            if is_decision_year(LvM_strategy, tstep)
                # Get candidate locations from strategy
                candidate_locs = filter_candidate_locations(LvM_strategy, tstep, state)

                # LvM_share is a @NamedTuple{weight::Float64, target_locs::Vector{String}}
                for LvM_share in domain.LvM_target_locations
                    share_candidate_locs = intersect(
                        candidate_locs, LvM_share.target_locs
                    )
                    if isempty(share_candidate_locs)
                        if ADRIA.decision.strategy_type(param_set, "iv_LvM") ==
                            PeriodicStrategy
                            @warn """
                                tstep $tstep: Deployment with PeriodicStrategy on
                                $(length(LvM_share.target_locs)) reefs with weight
                                $(LvM_share.weight) skipped because no candidate reefs
                                were found.
                            """
                        end
                        continue
                    end

                    share_candidate_loc_idx = findall(
                        in.(domain.loc_ids, Ref(share_candidate_locs))
                    )

                    if is_guided
                        # Update decision matrix with current conditions
                        update_criteria_values!(
                            LvM_decision_mat[location = At(share_candidate_locs)];
                            heat_stress=dhw_projection[share_candidate_loc_idx],
                            wave_stress=wave_projection[share_candidate_loc_idx],
                            coral_cover=current_loc_cover[share_candidate_loc_idx],
                            in_connectivity=in_conn[share_candidate_loc_idx],
                            out_connectivity=out_conn[share_candidate_loc_idx]
                        )

                        selected_LvM_ranks = select_locations(
                            LvM_pref,
                            LvM_decision_mat[location = At(share_candidate_locs)],
                            MCDA_approach,
                            LvM_min_iv_locs
                        )
                    else
                        # Unguided deployment, anywhere with available space > 0
                        selected_LvM_ranks = unguided_selection(
                            share_candidate_locs,
                            LvM_min_iv_locs,
                            vec(leftover_space_m²[share_candidate_loc_idx]),
                            depth_criteria[share_candidate_loc_idx];
                            rng=rng
                        )
                    end
                    if !isempty(selected_LvM_ranks)
                        log_val = is_guided ? (1:length(selected_LvM_ranks)) : 1.0
                        log_location_ranks[tstep, At(selected_LvM_ranks), At(:lvm)] .=
                            log_val
                        selected_LvM_loc_idx = findall(
                            in.(domain.loc_ids, Ref(selected_LvM_ranks))
                        )
                        last_LvM_deployment[selected_LvM_loc_idx] .= tstep
                    end

                    # Check if locations are selected
                    has_LvM_locs::Bool = !isempty(selected_LvM_ranks)

                    # Apply Larval Methods (assumed to occur after spawning)
                    if has_LvM_locs
                        LvM_loc_idx = findall(domain.loc_ids .∈ [selected_LvM_ranks])

                        @views available_space = leftover_space_m²[LvM_loc_idx]

                        locs_with_space = findall(available_space .> 0.0)

                        # If there are locations with space to select from, deploy what we can
                        if length(locs_with_space) > 0
                            LvM_loc_idx = LvM_loc_idx[locs_with_space]
                            available_space = available_space[locs_with_space]

                            @views LvM_proportional_increase, n_LvM_corals = distribute_LvM_corals(
                                vec_abs_k[LvM_loc_idx],
                                available_space,
                                iv_LvM_N_settlers * LvM_share.weight,
                                colony_areas[_CAq_size_groups],
                                prop_fecundity[:, LvM_loc_idx]
                            )

                            @views recruitment[:, LvM_loc_idx] .+=
                                LvM_proportional_increase
                            leftover_space_m²[LvM_loc_idx] .-=
                                vec(sum(Array(LvM_proportional_increase); dims=1)) .*
                                vec_abs_k[LvM_loc_idx]

                            # Log estimated number of corals moved
                            YLvM[tstep, :, LvM_loc_idx] .= n_LvM_corals
                        end

                        # Empty selected_LvM_ranks before the next iteration
                        selected_LvM_ranks = []
                    end
                end
            end
        end

        # Set next settlers DHW tolerance after reproducton and larval methods intervention
        settler_DHW_tolerance!(
            c_mean_t_1,
            c_mean_t,
            C_cover_t,
            vec_abs_k,
            TP_data,  # ! IMPORTANT: Pass in transition probability matrix, not connectivity!
            recruitment,
            fecundity_per_m²,
            param_set[At("heritability")],
            c_mean_tol_ceil
        )

        # Coral aquaculture
        # IDs of valid locations considering locations that have space for corals
        # locs_with_space = vec(leftover_space_m²) .> 0.0

        # if is_guided && CAq_decision_years[tstep] && (length(considered_locs) > 0)
        #     considered_locs = findall(_valid_locs .& locs_with_space)

        #     # Use modified projected DHW (may have been affected by fogging or shading)
        #     dhw_p = copy(dhw_scen)
        #     dhw_p[tstep, :] .= dhw_t

        #     dhw_projection = weighted_projection(dhw_p, tstep, plan_horizon, decay, tf)
        #     wave_projection = weighted_projection(wave_scen, tstep, plan_horizon, decay, tf)

        #     # Determine connectivity strength weighting by area.
        #     # Accounts for strength of connectivity where there is low/no coral cover
        #     in_conn, out_conn, _ = connectivity_strength(
        #         area_weighted_conn, vec(_loc_coral_cover), conn_cache
        #     )

        #     update_criteria_values!(
        #         decision_mat;
        #         heat_stress=dhw_projection[_valid_locs],
        #         wave_stress=wave_projection[_valid_locs],
        #         coral_cover=_loc_coral_cover[_valid_locs],  # Coral cover relative to `k`
        #         in_connectivity=in_conn[_valid_locs],  # area weighted connectivities for time `t`
        #         out_connectivity=out_conn[_valid_locs]
        #     )

        #     selected_CAq_ranks = select_locations(
        #         CAq_pref,
        #         decision_mat[location=locs_with_space[_valid_locs]],
        #         MCDA_approach,
        #         considered_locs,
        #         min_iv_locs
        #     )

        #     # Log rankings as appropriate
        #     if !isempty(selected_CAq_ranks)
        #         log_location_ranks[tstep, At(selected_CAq_ranks), At(:caq)] .=
        #             1:length(selected_CAq_ranks)
        #     end
        # elseif apply_CAq && CAq_decision_years[tstep]
        #     # Unguided deployment, coral-aquaculture/fog corals anywhere, so long as available space > 0
        #     selected_CAq_ranks = unguided_selection(
        #         domain.loc_ids,
        #         min_iv_locs,
        #         vec(leftover_space_m²),
        #         depth_criteria
        #     )

        #     log_location_ranks[tstep, At(selected_CAq_ranks), At(:caq)] .= 1.0

        #     # Estimate proportional change in cover to apply to cubes
        # end

        # Coral aquaculture location selection
        if !isnothing(CAq_strategy)
            state = if is_guided
                build_state(
                    domain,
                    CAq_strategy,
                    (
                        current_cover=current_loc_cover,
                        recent_cover_losses=recent_cover_losses,
                        last_deployment=iv_CAq_last_deployment
                    )
                )
            else
                nothing
            end

            if is_decision_year(CAq_strategy, tstep)
                # Get candidate locations from strategy
                candidate_locs = filter_candidate_locations(CAq_strategy, tstep, state)

                # CAq_share is a @NamedTuple{weight::Float64, target_locs::Vector{String}}
                for CAq_share in domain.CAq_target_locations
                    share_candidate_locs = intersect(
                        candidate_locs, CAq_share.target_locs
                    )
                    if isempty(share_candidate_locs)
                        if ADRIA.decision.strategy_type(param_set, "iv_CAq") ==
                            PeriodicStrategy
                            @warn """
                                tstep $tstep: Deployment with PeriodicStrategy on
                                $(length(CAq_share.target_locs)) reefs with weight
                                $(CAq_share.weight) skipped because no candidate reefs
                                were found.
                            """
                        end
                        continue
                    end

                    share_candidate_loc_idx = findall(
                        in.(domain.loc_ids, Ref(share_candidate_locs))
                    )

                    if is_guided
                        # Update decision matrix with current conditions
                        update_criteria_values!(
                            CAq_decision_mat[location = At(share_candidate_locs)];
                            heat_stress=dhw_projection[share_candidate_loc_idx],
                            wave_stress=wave_projection[share_candidate_loc_idx],
                            coral_cover=current_loc_cover[share_candidate_loc_idx],
                            in_connectivity=in_conn[share_candidate_loc_idx],
                            out_connectivity=out_conn[share_candidate_loc_idx]
                        )

                        # Build state for target locations only
                        selected_CAq_ranks = select_locations(
                            CAq_pref,
                            CAq_decision_mat[location = At(share_candidate_locs)],
                            MCDA_approach,
                            min_iv_locs
                        )
                    else
                        # Unguided deployment, coral-aquaculture/fog corals anywhere, so long as available space > 0
                        selected_CAq_ranks = unguided_selection(
                            share_candidate_locs,
                            min_iv_locs,
                            vec(leftover_space_m²[share_candidate_loc_idx]),
                            depth_criteria[share_candidate_loc_idx];
                            rng=rng
                        )
                    end
                    if !isempty(selected_CAq_ranks)
                        log_val = is_guided ? (1:length(selected_CAq_ranks)) : 1.0
                        log_location_ranks[tstep, At(selected_CAq_ranks), At(:caq)] .=
                            log_val
                        selected_CAq_loc_idx = findall(
                            in.(domain.loc_ids, Ref(selected_CAq_ranks))
                        )
                        iv_CAq_last_deployment[selected_CAq_loc_idx] .= tstep
                    end

                    # Check if locations are selected (can reuse previous selection)
                    has_CAq_locs::Bool = !isempty(selected_CAq_ranks)

                    # Apply coral aquaculture (assumed to occur after spawning)
                    if has_CAq_locs  # CAq_decision_years[tstep] &&
                        # Deploy coral aquaculture at selected locations
                        # Selected locations can fill up over time so avoid locations with no space´
                        CAq_loc_idx = findall(domain.loc_ids .∈ [selected_CAq_ranks])

                        available_space = leftover_space_m²[CAq_loc_idx]
                        locs_with_space = findall(available_space .> 0.0)

                        # If there are locations with space to select from, then deploy what we can
                        # Otherwise, do nothing.
                        if length(locs_with_space) > 0
                            # Calculate proportion for coral aquaculture based on current available space
                            CAq_loc_idx = CAq_loc_idx[locs_with_space]
                            available_space = available_space[locs_with_space]

                            proportional_increase, n_corals_CAq = distribute_CAq_corals(
                                vec_abs_k[CAq_loc_idx],
                                available_space,
                                CAq_volume .* CAq_share.weight,
                                colony_areas[_CAq_size_groups],
                                CAq_devices_per_m2
                            )

                            # Log estimated number of corals deployed via coral aquaculture
                            YCAq[tstep, :, CAq_loc_idx] .= n_corals_CAq'

                            # Add coral aquaculture to recruitment
                            recruitment[:, CAq_loc_idx] .+= proportional_increase
                            leftover_space_m²[CAq_loc_idx] .-=
                                vec(sum(Array(proportional_increase); dims=1)) .*
                                vec_abs_k[CAq_loc_idx]

                            update_tolerance_distribution!(
                                proportional_increase,
                                C_cover_t,
                                c_mean_t,
                                c_mean_reference[:, :, end],
                                c_std,
                                CAq_loc_idx,
                                _CAq_size_groups,
                                a_adapt,
                                c_mean_tol_ceil
                            )
                        end

                        # Empty selected_CAq_ranks before the next iteration
                        selected_CAq_ranks = []
                    end
                end
            end
        end

        # Apply disturbances

        # ΔC_cover_t should only hold changes in cover due to env. disturbances
        ΔC_cover_t .= C_cover_t

        # Store current c_mean before bleaching. This will be used in the next timestep to
        # update the settler's DHW tolerance dists.
        c_mean_t_1 .= c_mean_t

        # Bleaching mortality
        # Calculate and apply bleaching mortality
        # Bleaching typically occurs in the warmer months (November - February)
        # but mortality may occur within this timeframe or some time after it.
        # An earlier attempt to account for the cooling effect of storms / high wave
        # activity involves a linear reduction:
        #    `dhw_t .* (1.0 .- wave_scen[tstep, :])`
        # `wave_scen` is normalized to the maximum value found for the given wave scenario
        # so what causes 100% mortality can differ between runs.
        _c_cover_before_bleach = copy(C_cover_t)
        # Wave activity is said to also promote nutrient cycling, and improve growth.
        # Removed as I could not find any justification behind the implementation.
        eff_dhw_t .= effective_dhw_at_depth.(
            dhw_t, loc_data.depth_med, eff_dhw_base, eff_dhw_mix
        )
        bleaching_mortality!(
            C_cover_t,
            eff_dhw_t,
            c_std,
            c_mean_t,
            @view(bleach_dhw[(tstep - 1):tstep, :, :, :]),
            c_mean_tol_ceil;
            heat_lb=bleaching_dhw_floor,
            low_dhw_scale=bleaching_low_dhw_scale,
            low_dhw_threshold=bleaching_low_dhw_threshold
        )
        
        bleach_scalar = parse(Float64, get(ENV, "ADRIA_DEBUG_BLEACHING_SCALAR", "1.0"))
        if bleach_scalar != 1.0
            C_cover_t .= _c_cover_before_bleach .- (_c_cover_before_bleach .- C_cover_t) .* bleach_scalar
            clamp!(C_cover_t, 0.0, 1.0)
        end

        # Store current means to be used in future timesteps.
        # Use pre-coral-aquaculture snapshot so that a_adapt enhancement is not compounded
        # into the reference year-over-year.
        if a_adapt_ref > 0
            c_mean_reference[:, :, 2:end] .= c_mean_reference[:, :, 1:(end - 1)]
        end
        c_mean_reference[:, :, 1] .= c_mean_pre_CAq

        # Coral deaths due to selected cyclone scenario
        # Peak cyclone period is January to March
        # TODO: Update cyclone data to hold data for relevant functional groups
        _c_cover_before_cyclone = copy(C_cover_t)
        cyclone_mortality!(C_cover_t, cyclone_mortality_scen[tstep, :, :]')
        
        cyclone_scalar = parse(Float64, get(ENV, "ADRIA_DEBUG_CYCLONE_SCALAR", "1.0"))
        if cyclone_scalar != 1.0
            C_cover_t .= _c_cover_before_cyclone .- (_c_cover_before_cyclone .- C_cover_t) .* cyclone_scalar
            clamp!(C_cover_t, 0.0, 1.0)
        end

        # M2: rescale the state left by the reset year before this year's processes.
        if !isnothing(cots_reset) && cots_calendar_years[tstep] == cots_reset[1] + 1
            factors, other = apply_cots_state_reset!(cots_models, cots_reset[2], cots_reset[3],
                cots_site_areas; juveniles=cots_reset_juveniles, others=cots_reset_others)
            COTS_RESET_LOG[] = (year=cots_reset[1], groups=cots_reset[4], factors=factors, other=other,
                juveniles=cots_reset_juveniles)
        end
        # M1: this year's effort per reef (half before spawning, half after).
        cots_effort_on = !isnothing(cots_cull_effort) &&
            cots_effort_start <= cots_calendar_years[tstep] <= cots_effort_end
        cots_effort_hours(phase_models) = cots_effort_mode == "fixed" ?
            0.5 * cots_effort_multiplier .* cots_cull_effort :
            targeted_cull_hours(phase_models, 0.5 * cots_effort_multiplier * sum(cots_cull_effort),
                cots_reef_groups, cots_site_areas; threshold_ha=cots_effort_threshold)

        # COTS predation mortality (operates on relative cover)
        isnothing(cots_cull) || apply_cots_cull!(cots_models, cots_cull, Ycots_culled, tstep, :pre;
            groups=cots_cull_groups, areas=cots_cull_areas)
        if cots_effort_on
            hours_pre = cots_effort_hours(cots_models)
            cots_effort_hours_log[tstep] += sum(hours_pre)
            apply_cots_effort_cull!(cots_models, hours_pre, cots_reef_groups, cots_site_areas, Ycots_culled, tstep, :pre;
                cpue_max=cots_effort_cpue_max)
        end
        apply_cots_predation_diagnostic!(C_cover_t, cots_models, cots_prey_map,
            cots_calendar_years[tstep], cots_grazing_off_years;
            habitable_fraction=cots_habitable_fraction)
        isnothing(cots_cull) || apply_cots_cull!(cots_models, cots_cull, Ycots_culled, tstep, :post;
            groups=cots_cull_groups, areas=cots_cull_areas)
        if cots_effort_on
            hours_post = cots_effort_hours(cots_models)
            cots_effort_hours_log[tstep] += sum(hours_post)
            apply_cots_effort_cull!(cots_models, hours_post, cots_reef_groups, cots_site_areas, Ycots_culled, tstep, :post;
                cpue_max=cots_effort_cpue_max)
        end
        if !isnothing(cots_larval_survival)
            apply_cots_larval_survival!(cots_models, @view(cots_larval_survival[tstep, :]))
            Ycots_larval_survival[tstep, :] .= @view(cots_larval_survival[tstep, :])
        end
        isnothing(cots_settlement_modifier) ||
            apply_cots_settlement_modifier!(cots_models, cots_settlement_modifier)
        if !isnothing(cots_forced_adults)
            COTSMod.forced_grazing!(C_cover_t, @view(cots_forced_adults[tstep, :]),
                cots_prey_map; habitable_fraction=cots_forced_habitable,
                consumption_m2_per_adult_year=cots_params.consumption_m2_per_adult_year)
        end

        # Disperse COTS larvae between locations via connectivity
        _cots_imm_scalar = cots_internal_settlement_scalar(
            parse(Float64, get(ENV, "COTS_IMMIGRATION_SCALAR", "1.0")),
            cots_calendar_years[tstep],
            cots_internal_blackout
        )
        forcing_idx = 0
        if cots_temporal_mode == "mean" && !cots_separated_recruitment
            disperse_larvae!(cots_models, cots_conn; scalar=_cots_imm_scalar)
        elseif cots_temporal_mode == "mean"
            source_survival = cots_apply_larval_survival ?
                cots_mean_source_survival : nothing
            disperse_larvae_by_reef!(
                cots_models,
                cots_mean_reef_conn,
                cots_forcing.site_to_reef;
                scalar=_cots_imm_scalar,
                source_survival=source_survival,
                reef_areas=cots_forcing.reef_habitable_area_ha
            )
        else
            forcing_idx = cots_forcing_index(
                cots_forcing, tstep - 1, cots_temporal_mode, cots_connectivity_seed
            )
            source_survival = cots_apply_larval_survival ?
                @view(cots_forcing.source_survival[:, forcing_idx]) : nothing
            disperse_larvae_by_reef!(
                cots_models,
                cots_forcing.reef_connectivity[forcing_idx],
                cots_forcing.site_to_reef;
                scalar=_cots_imm_scalar,
                source_survival=source_survival,
                reef_areas=cots_separated_recruitment ?
                    cots_forcing.reef_habitable_area_ha : nothing
            )
            Ycots_forcing_year[tstep] = cots_forcing.years[forcing_idx]
            if cots_external_source_recruits > 0.0
                external_supply = cots_external_supply_by_site(
                    cots_forcing, forcing_idx, cots_external_source_recruits
                )
                apply_external_supply!(cots_models, external_supply)
            end
        end
        if cots_external_outbreak_production > 0.0
            calendar_year = cots_calendar_years[tstep]
            potential_supply, pelagic_supply = cots_external_outbreak_supply_by_site(
                cots_forcing,
                forcing_idx,
                calendar_year,
                cots_external_outbreak_production
            )
            supply_factor = cots_supply_window_factor(cots_supply_window, calendar_year)
            if cots_coral_coupling_start > 0
                # M3: domain coral at the end of the previous step (relative to habitable area).
                prev = vec(sum(@view(C_cover[tstep - 1, :, :, :]); dims=(1, 2)))
                c_dom = sum(prev .* cots_coral_weights) / sum(cots_coral_weights)
                if calendar_year < cots_coral_coupling_start
                    push!(cots_coral_reference, c_dom)
                else
                    isempty(cots_coral_reference) &&
                        error("COTS_EXTERNAL_CORAL_COUPLING_START_YEAR leaves no reference years")
                    g = cots_external_coral_factor(c_dom, mean(cots_coral_reference))
                    cots_coral_factor_log[tstep] = g
                    supply_factor *= g
                end
            end
            if supply_factor != 1.0
                potential_supply .*= supply_factor
                pelagic_supply .*= supply_factor
            end
            !cots_apply_larval_survival && (pelagic_supply = potential_supply)
            # Not in place: pelagic_supply may alias potential_supply.
            isnothing(cots_external_survival) ||
                (pelagic_supply = pelagic_supply .* cots_external_survival[tstep])
            apply_external_larval_supply!(
                cots_models, potential_supply, pelagic_supply
            )
        end

        # Mimic incoming larvae/recruitment from upstream outbreak (Cairns Initiation Box).
        # External pulses are off by default. When enabled, they add a relative
        # pulse scaled to each location's running maximum internal larval supply.
        for loc in 1:n_locs
            cots_max_larval_supply[loc] = max(cots_max_larval_supply[loc], cots_models[loc].N[1])
        end

        _cots_pulse_enabled = cots_enabled && get(ENV, "COTS_EXTERNAL_PULSE", "false") == "true"
        _pulse_elapsed = tstep - _cots_pulse_start
        _pulse_in_window = _pulse_elapsed >= 0 && (
            (_cots_pulse_repeat_interval <= 0 && _pulse_elapsed < _cots_pulse_duration) ||
            (_cots_pulse_repeat_interval > 0 && mod(_pulse_elapsed, _cots_pulse_repeat_interval) < _cots_pulse_duration)
        )
        if _cots_pulse_enabled && !isnothing(_cots_seed_locs) && _cots_pulse_relative_magnitude > 0.0 && _pulse_in_window
            pulse_vals = _cots_pulse_relative_magnitude .* cots_max_larval_supply
            apply_external_supply!(cots_models, _cots_seed_locs, pulse_vals)
        end

        # Log COTS populations and decomposed demographic/connectivity flows.
        cots_flows = cots_flow_diagnostics(cots_models)
        cots_hazard = cots_hazard_diagnostics(cots_models)
        cots_mechanism = cots_mechanism_diagnostics(cots_models)
        for loc in 1:n_locs
            Ycots[tstep, :, loc] .= cots_models[loc].N
            Ycots_bc[tstep, loc] = cots_models[loc].body_condition
            Ycots_flow[tstep, 1, loc] = cots_flows.local_fecundity[loc]
            Ycots_flow[tstep, 2, loc] = cots_flows.background_immigration[loc]
            Ycots_flow[tstep, 3, loc] = cots_flows.internal_immigration[loc]
            Ycots_flow[tstep, 4, loc] = cots_flows.external_immigration[loc]
            Ycots_flow[tstep, 5, loc] = cots_flows.maturation[loc]
            Ycots_flow[tstep, 6, loc] = cots_flows.retained_juveniles[loc]
            Ycots_flow[tstep, 7, loc] = cots_flows.settlement_gate[loc]
            Ycots_flow[tstep, 8, loc] = cots_flows.local_retention[loc]
            Ycots_flow[tstep, 9, loc] = cots_flows.pelagic_survivors[loc]
            Ycots_flow[tstep, 10, loc] = cots_flows.settled_recruits[loc]
            Ycots_flow[tstep, 11, loc] = cots_flows.external_potential[loc]
            Ycots_flow[tstep, 12, loc] = cots_flows.external_pelagic[loc]
            Ycots_size[tstep, 1, loc] = cots_flows.small_adults[loc]
            Ycots_size[tstep, 2, loc] = cots_flows.large_adults[loc]
            Ycots_size[tstep, 3, loc] = cots_flows.effective_breeders[loc]
            Ycots_size[tstep, 4, loc] = cots_flows.small_adult_growth[loc]
            Ycots_hazard[tstep, 1, loc] = cots_hazard.burden_ha[loc]
            Ycots_hazard[tstep, 2, loc] = cots_hazard.hazard[loc]
            Ycots_hazard[tstep, 3, loc] = cots_hazard.adult_deaths_ha[loc]
            Ycots_mechanism[tstep, 1, loc] = cots_mechanism.preferred_food_memory[loc]
            Ycots_mechanism[tstep, 2, loc] = cots_mechanism.food_survival[loc]
            Ycots_mechanism[tstep, 3, loc] = cots_mechanism.senescent_adults[loc]
            Ycots_mechanism[tstep, 4, loc] = cots_mechanism.senescent_deaths_ha[loc]
            Ycots_mechanism[tstep, 5, loc] = cots_mechanism.preferred_diet_share[loc]
            Ycots_mechanism[tstep, 6, loc] = cots_mechanism.juvenile_predation_ha[loc]
            Ycots_mechanism[tstep, 7, loc] = cots_mechanism.maturation_fraction[loc]
            Ycots_mechanism[tstep, 8, loc] = cots_mechanism.outbreak_clock[loc]
            Ycots_mechanism[tstep, 9, loc] = cots_mechanism.collapse_removed_ha[loc]
            Ycots_mechanism[tstep, 10, loc] = cots_mechanism.low_food_years[loc]
        end

        # Calculate survival_rate due to env. disturbances
        @inbounds for i in eachindex(survival_rate_cache)
            d = ΔC_cover_t[i]
            survival_rate_cache[i] = d == 0.0 ? 1.0 : C_cover_t[i] / d
        end
        @assert !any(>(1.0), survival_rate_cache) "Survival rate should be <= 1"

        for loc = 1:n_locs
            apply_mortality!(
                functional_groups[loc], @view(survival_rate_cache[:, :, loc])
            )
        end

        recruitment .*= (view(survival_rate_cache,:,1,:) .* habitable_areas)

        C_cover[tstep, :, :, :] .= C_cover_t

        # Track cover loss for reactive strategies
        if is_guided && _is_reactive && (tstep > 1)
            # Calculate proportional cover loss at each location
            # due to disturbances this timestep
            Δcover_loss_proportion .= 0.0
            for loc = 1:n_locs
                cover_before = sum(@view(ΔC_cover_t[:, :, loc]))
                cover_after = sum(@view(C_cover_t[:, :, loc]))

                if cover_before > 0.0
                    Δcover_loss_proportion[loc] =
                        (cover_before - cover_after) / cover_before
                end
            end

            # Store in rolling buffer
            push!(recent_cover_losses, Δcover_loss_proportion)
        end
    end

    # Could collate critical DHW threshold log for corals to reduce disk space...
    # dhw_tol_mean = dropdims(mean(dhw_tol_mean_log, dims=3), dims=3)
    # dhw_tol_mean_std = dropdims(mean(dhw_tol_std_log, dims=3), dims=3)
    # collated_dhw_tol_log = DataCube(cat(dhw_tol_mean, dhw_tol_mean_std, dims=3);
    #     timesteps=1:tf, species=corals.coral_id, stat=[:mean, :stdev])
    if log_dhw_tols
        # copy() prevents aliasing when the underlying buffer is reused across scenarios
        collated_dhw_tol_log = DataCube(
            copy(dhw_tol_mean_log); timesteps=1:tf, species=corals.coral_id, sites=1:n_locs
        )
    else
        collated_dhw_tol_log = false
    end

    if log_cover
        # copy() prevents aliasing when the underlying buffer is reused across scenarios
        collated_cover_log = DataCube(
            copy(cover_log); timesteps=1:tf, species=corals.coral_id, sites=1:n_locs
        )
    else
        collated_cover_log = false
    end

    # Set variables to nothing so garbage collector clears them
    # Leads to memory leak issues in multiprocessing contexts without these.
    wave_scen = nothing
    dhw_tol_mean_log = nothing
    cover_log = nothing

    return (
        raw=C_cover,
        CAq_log=YCAq,
        LvM_log=YLvM,
        Fog_log=Yfog,
        Shd_log=YShd,
        site_ranks=log_location_ranks,
        bleaching_mortality=bleach_dhw,
        coral_dhw_log=collated_dhw_tol_log,
        coral_cover_log=collated_cover_log,
        cots_log=Ycots,
        cots_condition_log=Ycots_bc,
        cots_flow_log=Ycots_flow,
        cots_size_log=Ycots_size,
        cots_hazard_log=Ycots_hazard,
        cots_mechanism_log=Ycots_mechanism,
        cots_forcing_year=Ycots_forcing_year,
        cots_calendar_year=cots_calendar_years,
        cots_culled_log=Ycots_culled,
        cots_larval_survival_log=Ycots_larval_survival,
        cots_external_coral_factor_log=cots_coral_factor_log,
        cots_reset_log=COTS_RESET_LOG[],
        cots_effort_hours_log=cots_effort_hours_log
    )
end
