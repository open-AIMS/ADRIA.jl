using CoralBlox

import CoralBlox:
    SizeClass,
    FunctionalGroup,
    reuse_buffers!,
    apply_mortality!,
    timestep!,
    coral_cover,
    max_projected_cover,
    linear_extension_scale_factors,
    LinearExtensionCache

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

    # Fail before running anything if a COTS scenario needs inputs the domain lacks
    validate_cots_settings(dom, scens)
    _validate_submodels(dom, scens)
    # `prepare` once per run; the caches are shared read-only by all scenarios and threads
    submodel_caches = _prepare_submodels(dom)

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

    # Cap the chunk size so that the per-thread batch buffers stay within a working budget
    # and every Zarr chunk file stays below the 2 GB limit enforced by BloscCompressor
    # (see `_max_chunk_by_size`).
    sub_specs = _submodel_log_specs(dom, scens)
    _max_chunk_by_size::Int = _max_chunk_by_size_for(dom, sub_specs)

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
    _bin_edges::Matrix{Float64} = bin_edges(; unit=:m).data

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
                            sim_cache=sc, submodel_caches=submodel_caches
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
                        dfx -> _collect_scenario_results(
                            dfx[1], dfx[3], functional_groups; submodel_caches
                        ),
                        batch
                    )
                    _write_batch!(data_store, first(batch)[2], batch_results)
                end
            else
                for batch in all_batches
                    batch_results = map(
                        dfx -> _collect_scenario_results(
                            dfx[1], dfx[3], functional_groups; submodel_caches
                        ),
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
- `coral_dhw_final` : Array, Mean DHW tolerance per group/size class/location (`n_groups × n_sizes × n_locs`) at the end of the run, i.e. after the last timestep's bleaching update. Unlike `coral_dhw_log`, whose row `t` is written before timestep `t`'s bleaching update, this includes the last bleaching update.
"""
function run_model(
    domain::Domain,
    param_set::Union{DataFrameRow,YAXArray};
    apply_allee_effect::Bool=true,
    submodel_caches=nothing
)
    n_locs::Int64 = domain.coral_growth.n_locs
    n_sizes::Int64 = domain.coral_growth.n_sizes
    n_groups::Int64 = domain.coral_growth.n_groups

    edges = bin_edges().data
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
        domain, param_set, functional_groups;
        apply_allee_effect=apply_allee_effect, submodel_caches=submodel_caches
    )
end
function run_model(
    domain::Domain,
    param_set::DataFrameRow,
    functional_groups::Vector{Vector{FunctionalGroup}};
    sim_cache=nothing,
    apply_allee_effect::Bool=true,
    submodel_caches=nothing
)::NamedTuple
    setup()
    ps = DataCube(Vector(param_set); factors=names(param_set))
    return run_model(
        domain,
        ps,
        functional_groups;
        sim_cache=sim_cache,
        apply_allee_effect=apply_allee_effect,
        submodel_caches=submodel_caches
    )
end
function run_model(
    domain::Domain,
    param_set::YAXArray,
    functional_groups::Vector{Vector{FunctionalGroup}};
    sim_cache=nothing,
    apply_allee_effect::Bool=true,
    submodel_caches=nothing
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

    # Submodels attached to the domain (none by default). Called directly (not through
    # `run_scenarios`), validate and prepare in-line, like `sim_cache`.
    if isnothing(submodel_caches) && !isempty(attached_submodels(domain))
        _validate_submodels(
            domain, DataFrame(permutedims(Vector(param_set.data)), factor_names)
        )
        submodel_caches = _prepare_submodels(domain)
    end

    rng = set_random_seed(param_set; excluded=_submodel_seed_excluded(domain))

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
    else
        MCDA_approach = nothing
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

    CAq_t_locs = vcat(getproperty.(domain.CAq_target_locations, :target_locs)...)
    LvM_t_locs = vcat(getproperty.(domain.LvM_target_locations, :target_locs)...)
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
        # Only used by guided selection; defined so they can be passed to the phases
        CAq_pref = CAq_decision_mat = nothing
        Fog_pref = Fog_decision_mat = nothing
        LvM_pref = LvM_decision_mat = nothing
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

    # Cache matrix to store potential settlers
    potential_settlers = zeros(size(fec_scope)...)
    _linear_extensions = _to_group_size(domain.coral_growth, corals.linear_extension)
    _bin_edges = bin_edges().data
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
    n_cb_calib_groups::Int64 = length(unique_cb_calib_groups)
    cb_calib_group_masks::BitMatrix = falses(n_locs, n_cb_calib_groups)
    for (idx, cb_calib_group) in enumerate(unique_cb_calib_groups)
        cb_calib_group_masks[:, idx] .= domain.loc_data.CB_CALIB_GROUPS .== cb_calib_group
    end

    # Index into unique_biogroups for each location
    loc_cb_calib_group_idxs::Vector{Int64} = [
        findfirst(x -> x == cb_group, unique_cb_calib_groups)
        for cb_group in domain.loc_data.CB_CALIB_GROUPS
    ]

    is_growth_acc_mask = occursin.("growth_acceleration", factor_names)
    growth_acc_steepness::Vector{Float64} =
        param_set[is_growth_acc_mask .&& occursin.("steepness", factor_names)].data
    growth_acc_height::Vector{Float64} =
        param_set[is_growth_acc_mask .&& occursin.("height", factor_names)].data
    growth_acc_midpoint::Vector{Float64} =
        param_set[is_growth_acc_mask .&& occursin.("midpoint", factor_names)].data

    # linear_extension scale factors with dimensions (cb_calib_groups ⋅ functional_groups)
    _linear_extension_scale_factors = reshape(
        param_set[occursin.("linear_extension_scale", factor_names)].data,
        n_cb_calib_groups,
        n_groups
    )

    # mb_rate scale factors with dimensions (cb_calib_groups ⋅ functional_groups)
    _mb_rate_scale_factors = reshape(
        param_set[occursin.("mb_rate_scale", factor_names)].data,
        n_cb_calib_groups,
        n_groups
    )

    biogrp_lin_ext::Array{Float64,3} = repeat(_linear_extensions, 1, 1, n_cb_calib_groups)
    biogrp_survival::Array{Float64,3} = repeat(survival_rate, 1, 1, n_cb_calib_groups)
    for i = 1:n_cb_calib_groups
        biogrp_lin_ext[:, :, i] .*= _linear_extension_scale_factors[i, :]
        biogrp_survival[:, :, i] .= apply_survival_scaling(
            biogrp_survival[:, :, i], _mb_rate_scale_factors[i, :]
        )
    end

    # Pre-compute Δd matrices once for the entire simulation — _bin_edges is constant,
    # so this avoids repeated allocation inside linear_extension_scale_factors each timestep.
    _le_cache = LinearExtensionCache(_bin_edges)

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

    # Optional reduction of bleaching mortality at low effective DHW (see `LowDHWBleaching`)
    low_dhw_bleach_reduction::Float64 = param_set[At("low_dhw_bleach_reduction")]
    low_dhw_bleach_threshold::Float64 = param_set[At("low_dhw_bleach_threshold")]

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
    c_mean_reference .= c_mean_t[:, JUVENILE_SIZE_CLASS, :]

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

    # Optional submodels attached to the domain (ADRIAinterfaces); `()` when none
    submodel_slots = _init_submodels(domain, param_set, factor_names, submodel_caches, rng)
    has_submodels = !isempty(submodel_slots)

    # Per-cell cover before the submodel phase, to enforce decrease-only; empty (no
    # allocation of note, same type) when nothing is attached
    submodel_cover_before = zeros(has_submodels ? size(C_cover_t) : (0, 0, 0))

    # COTS submodel (COTSMod.jl), run only when the scenario sets `cots_on` = 1
    cots = param_set[At("cots_on")] == 1 ? CotsRun(domain, param_set, rng) : nothing
    # COTS densities (timesteps x 3 age classes x locations, COTS ha^-1) and body condition
    Ycots = zeros(tf, 3, n_locs)
    Ycots_bc = zeros(tf, n_locs)

    # State used by the growth phase (arrays are mutated in place by `_grow_corals!`)
    growth_state = (;
        # Habitable locations and their areas
        habitable_locs, habitable_loc_idxs, habitable_loc_areas, habitable_loc_areas′,
        vec_abs_k,

        # Carrying capacity (maximum projected cover)
        habitable_max_projected_cover, cache_habitable_max_projected_cover,

        # Recruitment scratch space
        _recruitment_col_sum_2d, _recruitment_col_sum,

        # Growth-constraint masks and caches
        agg_cover_above_threshold_mask, growth_threshold_mask_cache, cover_threshold_mask,
        relative_habitable_cover_cache, growth_constraints, net_growth_rates,

        # Cover thresholds for the growth transition
        cover_transition_lb, cover_transition_ub,

        # Calibration groups (per-location masks and indices) and their scaled rates
        cb_calib_group_masks, n_cb_calib_groups, loc_cb_calib_group_idxs, biogrp_lin_ext,
        biogrp_survival, _le_cache,

        # Growth acceleration curve parameters
        growth_acc_height, growth_acc_midpoint, growth_acc_steepness,

        # Functional group definitions
        functional_groups
    )

    # State used by the reproduction phase (arrays are mutated in place by `_reproduce!`)
    reproduction_state = (;
        # Locations and areas
        n_locs, habitable_loc_idxs, habitable_areas, vec_abs_k, _loc_coral_cover,

        # Fecundity and settlement
        fec_scope, fecundity_per_m², prop_fecundity, potential_settlers,
        basal_area_per_settler, conn, sim_params,

        # Space and recruitment scratch space
        leftover_space_m², _recruitment_col_sum_2d, _recruitment_col_sum,

        # Heat tolerance snapshot taken before coral aquaculture
        c_mean_pre_CAq
    )

    # State shared by the intervention phases and per-intervention constants; arrays are
    # mutated in place by the phases
    intervention_state = (;
        # Domain, run mode and random numbers
        domain, is_guided, MCDA_approach, rng,

        # Projections and current state
        current_loc_cover, dhw_projection, wave_projection, leftover_space_m²,

        # Intervention tracking
        recent_cover_losses, log_location_ranks, iv_CAq_last_deployment
    )
    lvm_state = (;
        # Strategy and location selection
        LvM_strategy, LvM_pref, LvM_decision_mat, LvM_min_iv_locs, depth_criteria,

        # Deployment settings and tracking
        iv_LvM_N_settlers, last_LvM_deployment, YLvM,

        # Colony sizes and areas
        colony_areas, _CAq_size_groups, vec_abs_k
    )
    caq_state = (;
        # Strategy and location selection
        CAq_strategy, CAq_pref, CAq_decision_mat, min_iv_locs, depth_criteria,

        # Deployment settings and tracking
        CAq_volume, CAq_devices_per_m2, YCAq,

        # Colony sizes and areas
        colony_areas, _CAq_size_groups, vec_abs_k,

        # Heat tolerance of the outplanted corals
        c_mean_reference, c_std, a_adapt, c_mean_tol_ceil
    )
    fog_state = (;
        # Strategy and location selection
        Fog_strategy, Fog_pref, Fog_decision_mat, min_iv_locs,

        # Deployment settings and tracking
        fogging, last_Fog_deployment, Yfog
    )

    # State used by the disturbance phase
    disturbance_state = (;
        # Location data and cyclone mortality
        loc_data, cyclone_mortality_scen,

        # DHW and bleaching
        eff_dhw_base, eff_dhw_mix, bleach_dhw, low_dhw_bleach_reduction,
        low_dhw_bleach_threshold,

        # Heat tolerance
        c_std, c_mean_tol_ceil, c_mean_reference, c_mean_pre_CAq, a_adapt_ref
    )

    for tstep::Int64 = 2:tf
        _grow_corals!(C_cover, C_cover_t, recruitment, growth_state, tstep)

        _adapt_and_log!(
            C_cover, c_mean_t, net_growth_rates, c_mean_tol_ceil, dhw_tol_mean_log,
            cover_log, _permuted_buf, log_dhw_tols, log_cover, tstep, tf
        )

        _reproduce!(
            recruitment, C_cover_t, c_mean_t, reproduction_state, rng, apply_allee_effect
        )

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

        _apply_shading!(
            dhw_t, YShd, dhw_scen, Shd_locs_mask, Shd_decision_years, apply_Shd, Shd, tstep
        )

        in_conn, out_conn = nothing, nothing
        if is_guided
            in_conn, out_conn = _update_guided_state!(
                dhw_projection, wave_projection, current_loc_cover, dhw_p, dhw_scen, dhw_t,
                wave_scen, wave_idx, area_weighted_conn, _loc_coral_cover, conn_cache,
                C_cover_t, plan_horizon, decay, tf, tstep
            )
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

        _apply_fogging!(fog_state, intervention_state, in_conn, out_conn, dhw_t, tstep)

        _apply_larval_methods!(
            lvm_state, intervention_state, in_conn, out_conn, recruitment, prop_fecundity,
            param_set, tstep
        )

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

        _apply_coral_aquaculture!(
            caq_state, intervention_state, in_conn, out_conn, recruitment, C_cover_t,
            c_mean_t, param_set, tstep
        )

        _apply_disturbances!(
            disturbance_state, C_cover_t, ΔC_cover_t, c_mean_t, c_mean_t_1, dhw_t,
            eff_dhw_t, tstep
        )

        _step_cots!(
            cots, Ycots, Ycots_bc, C_cover_t, tstep,
            Int64(domain.env_layer_md.timeframe[tstep])
        )

        if has_submodels
            _step_submodel_hooks!(
                submodel_slots, submodel_cover_before, C_cover_t, dhw_t, vec_abs_k, tstep,
                Int64(domain.env_layer_md.timeframe[tstep])
            )
        end

        _propagate_survival!(
            survival_rate_cache, ΔC_cover_t, C_cover_t, functional_groups, recruitment,
            habitable_areas, C_cover, n_locs, tstep
        )

        # Track cover loss for reactive strategies
        if is_guided && _is_reactive && (tstep > 1)
            _track_cover_loss!(
                recent_cover_losses, Δcover_loss_proportion, ΔC_cover_t, C_cover_t, n_locs
            )
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
        coral_dhw_final=copy(c_mean_t),
        cots_log=Ycots,
        cots_condition_log=Ycots_bc,
        submodel_logs=_submodel_logs(submodel_slots)
    )
end
