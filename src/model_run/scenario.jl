"""Scenario running functions"""

using CoralBlox

import CoralBlox: SizeClass, FunctionalGroup
import CoralBlox: coral_cover, max_projected_cover, linear_extension_scale_factors
import CoralBlox: reuse_buffers!, apply_mortality!, timestep!

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

using .decision

"""
    _reshape_init_cover(data::AbstractMatrix{<:Union{Float32, Float64}})

Reshape initial coral cover of shape [groups_and_sizes ⋅ locations] to shape
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
    _to_group_size(growth_spec::CoralDetails, data::AbstractVector{<:Union{Float32, Float64}})::Matrix{<:Union{Float32, Float64}}

Reshape vector to shape [functional_groups ⋅ sizes]
"""
function _to_group_size(
    growth_spec::CoralDetails, data::AbstractVector{T}
)::Matrix{T} where {T<:Union{AbstractFloat,Bool}}
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
```julia-repl
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
    dom::Domain, scens::DataFrame, RCP::String; show_progress=true, remove_workers=true
)::ResultSet
    return run_scenarios(dom, scens, [RCP]; show_progress, remove_workers)
end
function run_scenarios(
    dom::Domain,
    scens::DataFrame,
    RCP::Vector{String};
    show_progress=true,
    remove_workers=true
)::ResultSet
    # Initialize ADRIA configuration options
    setup()

    # Sort RCPs so the dataframe order match the output filepath
    RCP = sort(RCP)

    @info "Running $(nrow(scens)) scenarios over $(length(RCP)) RCPs: $RCP"

    # Cross product between rcps and param_df to have every row of param_df for each rcp
    rcps_df = DataFrame(; RCP=parse.(Int64, RCP))
    scenarios_df = crossjoin(scens, rcps_df)
    sort!(scenarios_df, :RCP)

    @info "Setting up Result Set"
    dom, data_store = ADRIA.setup_result_store!(dom, scenarios_df)

    # Convert DataFrame to named matrix for faster iteration
    scenarios_matrix::YAXArray = DataCube(
        Matrix(scenarios_df); scenarios=1:nrow(scenarios_df), factors=names(scenarios_df)
    )

    # The idea to initialise the functional groups here is so that each scenario run can
    # reuse the memory. After a few scenarios runs there will be very few new
    # allocations as buffers grow larger and are not exceeded
    #
    # The problem is that even if a subsequent run only has to reallocate once, the
    # buffer is so large it takes very long regardless.
    #
    # This is also forced me to added an argument to run_model which breaks
    # encapsulation, so its quite ugly
    n_locs::Int64 = n_locations(dom)
    n_sizes::Int64 = dom.coral_details.n_sizes
    n_groups::Int64 = dom.coral_details.n_groups
    _bin_edges::Matrix{Float64} = bin_edges()
    functional_groups = [
        FunctionalGroup.(
            eachrow(_bin_edges[:, 1:(end - 1)]),
            eachrow(_bin_edges[:, 2:end]),
            eachrow(zeros(n_groups, n_sizes))
        ) for _ in 1:n_locs
    ]

    para_threshold::Int64 =
        ((typeof(dom) == RMEDomain) || (typeof(dom) == ReefModDomain)) ? 8 : 256
    active_cores::Int64 = parse(Int64, ENV["ADRIA_NUM_CORES"])
    parallel::Bool =
        (parse(Bool, ENV["ADRIA_DEBUG"]) == false) && (active_cores > 1) &&
        (nrow(scens) >= para_threshold)
    if parallel && nworkers() == 1
        @info "Setting up parallel processing..."
        spinup_time = @elapsed begin
            _setup_workers()

            # Load ADRIA on workers and define helper function
            # Note: Workers do not share the same cache in parallel workloads.
            #       Previously, cache would be reserialized so each worker has access to
            #       a separate cache.
            #       Using CachingPool() resolves the the repeated reserialization but it
            #       seems each worker was then attempting to use the same cache, causing the
            #       Julia kernel to crash in multi-processing contexts.
            #       Getting each worker to create its own cache reduces serialization time
            #       (at the cost of increased run time) but resolves the kernel crash issue.
            @sync @async @everywhere @eval begin
                using ADRIA
                func = (dfx) -> run_scenario(dfx..., functional_groups, data_store)
            end
        end

        @info "Time taken to spin up workers: $(round(spinup_time; digits=2)) seconds"
    end

    if parallel
        # Define local helper
        func = (dfx) -> run_scenario(dfx..., functional_groups, data_store)

        try
            for rcp in RCP
                run_msg = "Running $(nrow(scens)) scenarios for RCP $rcp"

                # Switch RCPs so correct data is loaded
                dom = switch_RCPs!(dom, rcp)
                target_rows = findall(
                    scenarios_matrix[factors=At("RCP")] .== parse(Float64, rcp)
                )
                scen_args = _scenario_args(dom, scenarios_matrix, rcp, length(target_rows))

                if show_progress
                    @showprogress desc = run_msg dt = 4 pmap(
                        func, CachingPool(workers()), scen_args
                    )
                else
                    pmap(func, CachingPool(workers()), scen_args)
                end
            end
        catch err
            # Remove hanging workers if any error occurs
            _remove_workers()
            rethrow(err)
        end
    else
        # Define local helper
        func = dfx -> run_scenario(dfx..., functional_groups, data_store)

        for rcp in RCP
            run_msg = "Running $(nrow(scens)) scenarios for RCP $rcp"

            # Switch RCPs so correct data is loaded
            dom = switch_RCPs!(dom, rcp)
            scen_args = _scenario_args(
                dom, scenarios_matrix, rcp, size(scenarios_matrix, 1)
            )

            if show_progress
                @showprogress desc = run_msg dt = 4 map(func, scen_args)
            else
                map(func, scen_args)
            end
        end
    end

    if parallel && remove_workers
        _remove_workers()
    end

    return load_results(_result_location(dom, RCP))
end

function _scenario_args(dom, scenarios_matrix::YAXArray, rcp::String, n::Int)
    target_rows = findall(scenarios_matrix[factors=At("RCP")] .== parse(Float64, rcp))
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
- `domain` : Domain
- `idx` : Scenario index

# Returns
Nothing
"""
function run_scenario(
    domain::Domain,
    idx::Int64,
    scenario::Union{AbstractVector,DataFrameRow},
    functional_groups::Vector{Vector{FunctionalGroup}}, # additional argument for reusable buffer
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

    result_set = run_model(domain, scenario, functional_groups)

    # Capture results to disk
    # Set values below threshold to 0 to save space
    threshold = parse(Float32, ENV["ADRIA_THRESHOLD"])

    rs_raw::Array{Float64} = result_set.raw
    vals = relative_cover(rs_raw)
    vals[vals .< threshold] .= 0.0
    data_store.relative_cover[:, :, idx] .= vals

    vals = absolute_shelter_volume(rs_raw, loc_k_area(domain), scenario)
    vals[vals .< threshold] .= 0.0
    data_store.absolute_shelter_volume[:, :, idx] .= vals
    vals = relative_shelter_volume(rs_raw, loc_k_area(domain), scenario)
    vals[vals .< threshold] .= 0.0
    data_store.relative_shelter_volume[:, :, idx] .= vals

    coral_spec::DataFrame = to_coral_spec(scenario)
    vals = relative_juveniles(rs_raw, coral_spec)
    vals[vals .< threshold] .= 0.0
    data_store.relative_juveniles[:, :, idx] .= vals

    vals = juvenile_indicator(rs_raw, coral_spec, loc_k_area(domain))
    vals[vals .< threshold] .= 0.0
    data_store.juvenile_indicator[:, :, idx] .= vals

    vals = relative_taxa_cover(rs_raw, loc_k_area(domain), domain.coral_details.n_groups)
    vals[vals .< threshold] .= 0.0
    data_store.relative_taxa_cover[:, :, idx] .= vals

    vals = relative_loc_taxa_cover(
        rs_raw, loc_k_area(domain), domain.coral_details.n_groups
    )

    vals = coral_evenness(vals.data)
    vals[vals .< threshold] .= 0.0
    data_store.coral_evenness[:, :, idx] .= vals

    # Store logs
    c_dim = Base.ndims(result_set.raw) + 1
    log_stores = (:site_ranks, :seed_log, :fog_log, :shade_log, :coral_dhw_log)
    for k in log_stores
        if k == :seed_log || k == :site_ranks
            concat_dim = c_dim
        else
            concat_dim = c_dim - 1
        end

        vals = getfield(result_set, k)

        try
            vals[vals .< threshold] .= Float32(0.0)
        catch err
            err isa MethodError ? nothing : rethrow(err)
        end

        if k == :seed_log
            getfield(data_store, k)[:, :, :, idx] .= vals
        elseif k == :site_ranks
            if !isnothing(data_store.site_ranks)
                # Squash site ranks down to average rankings over environmental repeats
                data_store.site_ranks[:, :, :, idx] .= vals
            end
        elseif k == :coral_dhw_log
            # Only log coral DHW tolerances if in debug mode
            if parse(Bool, ENV["ADRIA_DEBUG"]) == true
                getfield(data_store, k)[:, :, :, idx] .= vals
            end
        else
            getfield(data_store, k)[:, :, idx] .= vals
        end
    end

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

function _loc_coral_cover(C_cover_t::Array{Float64,3})
    return dropdims(sum(C_cover_t; dims=(1, 2)); dims=(1, 2))
end
