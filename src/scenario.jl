"""Scenario running functions"""

using DynamicCoralCoverModel
import DynamicCoralCoverModel.blocks_model: CoverBlock
import DynamicCoralCoverModel.blocks_model: SizeClass

using ADRIA.metrics:
    relative_cover,
    relative_loc_taxa_cover,
    total_absolute_cover,
    absolute_shelter_volume,
    relative_shelter_volume
using ADRIA.metrics: relative_juveniles, relative_taxa_cover, juvenile_indicator
using ADRIA.metrics: coral_evenness
using ADRIA.decision

"""
    setup_cache(domain::Domain)::NamedTuple

Establish tuple of matrices/vectors for use as reusable data stores to avoid repeated memory allocations.
"""
function setup_cache(domain::Domain)::NamedTuple

    # Simulation constants
    n_locs::Int64 = domain.coral_growth.n_locs
    n_group_and_size::Int64 = domain.coral_growth.n_group_and_size
    n_sizes::Int64 = domain.coral_growth.n_sizes
    n_groups::Int64 = domain.coral_growth.n_groups
    tf = length(timesteps(domain))

    cache = (
        # sf=zeros(n_groups, n_sites),  # stressed fecundity, commented out as it is disabled
        fec_all=zeros(n_groups, n_sizes, n_locs),  # all fecundity
        fec_scope=zeros(n_groups, n_locs),  # fecundity scope
        recruitment=zeros(n_groups, n_locs),  # coral recruitment
        dhw_step=zeros(n_locs),  # DHW for each time step
        cov_tmp=zeros(n_groups, n_sizes, n_locs),  # Cover for previous timestep
        depth_coeff=zeros(n_locs),  # store for depth coefficient
        site_area=Matrix{Float64}(site_area(domain)'),  # area of locations
        site_k_area=Matrix{Float64}(site_k_area(domain)'),  # location carrying capacity
        wave_damage=zeros(tf, n_group_and_size, n_locs),  # damage coefficient for each size class
        dhw_tol_mean_log=zeros(tf, n_group_and_size, n_locs),  # tmp log for mean dhw tolerances
    )

    return cache
end

"""
    _cover_to_groups(growth_spec::CoralGrowth, data::AbstractVector{<:Union{Float32, Float64}})::Matrix{<:Union{Float32, Float64}}

Reshape coral cover to shape [sizes ⋅ locations ⋅ n_groups]
"""
function _cover_to_groups(growth_spec::CoralGrowth, data::AbstractMatrix{<:Union{Float32, Float64}})
    return Array(reshape(
        data,
        (growth_spec.n_sizes, growth_spec.n_locs, growth_spec.n_groups)
    ))
end

"""
    _group_cover_locs(growth_spec::CoralGrowth, data::AbstractMatrix{<:Union{Float32, Float64}})

Reshape coral cover of shape [groups_and_sizes ⋅ locations] to shape [groups ⋅ sizes ⋅ locations]
"""
function _group_cover_locs(growth_spec::CoralGrowth, data::AbstractMatrix{<:Union{Float32, Float64}})
    return permutedims(reshape(
        data,
        (growth_spec.n_sizes, growth_spec.n_groups, growth_spec.n_locs)
    ), (2, 1, 3))
end

"""
    _flatten_cover(growth_spec::CoralGrowth, data::AbstractArray{<:Union{Float32, Float64}})::Matrix{<:Union{Float32, Float64}}

Reshape coral cover of shape [groups ⋅ sizes ⋅ locations] to shape [groups_and_sizes ⋅ locations]
"""
function _flatten_cover(growth_spec::CoralGrowth, data::AbstractArray{<:Union{Float32, Float64}})::Matrix{<:Union{Float32, Float64}}
    return reshape(
        permutedims(data, [2, 1, 3]),
        (growth_spec.n_group_and_size, growth_spec.n_locs)
    )
end

"""
    _to_group_size(growth_spec::CoralGrowth, data::AbstractVector{<:Union{Float32, Float64}})::Matrix{<:Union{Float32, Float64}}

Reshape vector to shape [groups ⋅ sizes]
"""
function _to_group_size(growth_spec::CoralGrowth, data::AbstractVector{T})::Matrix{T} where {T<:Union{Float32, Float64, Bool}}
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
    remove_workers=true,
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

    para_threshold = ((typeof(dom) == RMEDomain) || (typeof(dom) == ReefModDomain)) ? 8 : 256
    active_cores::Int64 = parse(Int64, ENV["ADRIA_NUM_CORES"])
    parallel = (parse(Bool, ENV["ADRIA_DEBUG"]) == false)  && (active_cores > 1) && (nrow(scens) >= para_threshold)
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
                func = (dfx) -> run_scenario(dfx..., data_store)
            end
        end

        @info "Time taken to spin up workers: $(round(spinup_time; digits=2)) seconds"
    end

    if parallel
        # Define local helper
        func = (dfx) -> run_scenario(dfx..., data_store)

        try
            for rcp in RCP
                run_msg = "Running $(nrow(scens)) scenarios for RCP $rcp"

                # Switch RCPs so correct data is loaded
                dom = switch_RCPs!(dom, rcp)
                target_rows = findall(scenarios_matrix[factors=At("RCP")] .== parse(Float64, rcp))
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
        func = dfx -> run_scenario(dfx..., data_store)

        for rcp in RCP
            run_msg = "Running $(nrow(scens)) scenarios for RCP $rcp"

            # Switch RCPs so correct data is loaded
            dom = switch_RCPs!(dom, rcp)
            scen_args = _scenario_args(dom, scenarios_matrix, rcp, size(scenarios_matrix, 1))

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
    run_scenario(domain::Domain, idx::Int64, scenario::Union{AbstractVector, DataFrameRow}, data_store::NamedTuple)::Nothing
    run_scenario(domain::Domain, idx::Int64, scenario::Union{AbstractVector, DataFrameRow}, domain::Domain, data_store::NamedTuple)::Nothing
    run_scenario(domain::Domain, scenario::Union{AbstractVector, DataFrameRow})::NamedTuple
    run_scenario(domain::Domain, scenario::NamedTuple)::NamedTuple

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
    data_store::NamedTuple,
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

    result_set = run_model(domain, scenario)

    # Capture results to disk
    # Set values below threshold to 0 to save space
    threshold = parse(Float32, ENV["ADRIA_THRESHOLD"])

    rs_raw::Array{Float64} = result_set.raw
    vals = relative_cover(rs_raw)
    vals[vals.<threshold] .= 0.0
    data_store.relative_cover[:, :, idx] .= vals

    vals .= absolute_shelter_volume(rs_raw, site_k_area(domain), scenario)
    vals[vals.<threshold] .= 0.0
    data_store.absolute_shelter_volume[:, :, idx] .= vals

    vals .= relative_shelter_volume(rs_raw, site_k_area(domain), scenario)
    vals[vals.<threshold] .= 0.0
    data_store.relative_shelter_volume[:, :, idx] .= vals

    coral_spec::DataFrame = to_coral_spec(scenario)
    vals .= relative_juveniles(rs_raw, coral_spec)
    vals[vals.<threshold] .= 0.0
    data_store.relative_juveniles[:, :, idx] .= vals

    vals .= juvenile_indicator(rs_raw, coral_spec, site_k_area(domain))
    vals[vals.<threshold] .= 0.0
    data_store.juvenile_indicator[:, :, idx] .= vals

    vals = relative_taxa_cover(rs_raw, site_k_area(domain), domain.coral_growth.n_groups)
    vals[vals.<threshold] .= 0.0
    data_store.relative_taxa_cover[:, :, idx] .= vals

    vals = relative_loc_taxa_cover(rs_raw, site_k_area(domain), domain.coral_growth.n_groups)
    vals = coral_evenness(vals.data)
    vals[vals.<threshold] .= 0.0
    data_store.coral_evenness[:, :, idx] .= vals

    # Store raw results if no metrics specified
    # if length(metrics) == 0
    #     data_store.raw[:, :, :, idx] .= r.raw
    # end

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
            vals[vals.<threshold] .= Float32(0.0)
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

"""
    run_model(domain::Domain, param_set::Union{NamedTuple,DataFrameRow})::NamedTuple

Core scenario running function.

# Returns
NamedTuple of collated results
- `raw` : Array, Coral cover relative to k area
- `seed_log` : Array, Log of seeded locations
- `fog_log` : Array, Log of fogged locations
- `shade_log` : Array, Log of shaded locations
- `site_ranks` : Array, Log of location rankings
- `bleaching_mortality` : Array, Log of mortalities caused by bleaching
- `coral_dhw_log` : Array, Log of DHW tolerances / adaptation over time (only logged in debug mode)
"""
function run_model(domain::Domain, param_set::DataFrameRow)::NamedTuple
    ps = DataCube(Vector(param_set); factors=names(param_set))
    return run_model(domain, ps)
end
function run_model(domain::Domain, param_set::YAXArray)::NamedTuple
    p = domain.coral_growth.ode_p
    corals = to_coral_spec(param_set)
    cache = setup_cache(domain)

    # Set random seed using intervention values
    # TODO: More robust way of getting intervention/criteria values
    rnd_seed_val::Int64 = floor(Int64, sum(param_set[Where(x -> x != "RCP")]))  # select everything except RCP
    Random.seed!(rnd_seed_val)

    dhw_idx::Int64 = Int64(param_set[At("dhw_scenario")])
    wave_idx::Int64 = Int64(param_set[At("wave_scenario")])
    cyclone_mortality_idx::Int64 = Int64(param_set[At("cyclone_mortality_scenario")])

    # Extract environmental data
    dhw_scen = @view(domain.dhw_scens[:, :, dhw_idx])

    # TODO: Better conversion of Ub to wave mortality
    #       Currently scaling significant wave height by its max to non-dimensionalize values
    wave_scen = copy_datacube(domain.wave_scens[:, :, wave_idx])
    wave_scen .= wave_scen ./ maximum(wave_scen)
    replace!(wave_scen, Inf => 0.0, NaN => 0.0)

    cyclone_mortality_scen = @view(
        domain.cyclone_mortality_scens[:, :, :, cyclone_mortality_idx]
    )

    # Environment variables are stored as strings, so convert to bool for use
    in_debug_mode = parse(Bool, ENV["ADRIA_DEBUG"]) == true

    # Sim constants
    sim_params = domain.sim_constants
    tf::Int64 = size(dhw_scen, 1)
    n_locs::Int64 = domain.coral_growth.n_locs
    n_groups::Int64 = domain.coral_growth.n_groups
    n_sizes::Int64 = domain.coral_growth.n_sizes
    n_group_and_size::Int64 = domain.coral_growth.n_group_and_size

    # Locations to intervene
    min_iv_locs::Int64 = param_set[At("min_iv_locations")]
    max_members::Int64 = param_set[At("cluster_max_member")]

    # Years to start seeding/shading/fogging
    seed_start_year::Int64 = param_set[At("seed_year_start")]
    shade_start_year::Int64 = param_set[At("shade_year_start")]
    fog_start_year::Int64 = param_set[At("fog_year_start")]

    fogging::Real = param_set[At("fogging")]  # proportion of bleaching mortality reduction through fogging
    srm::Real = param_set[At("SRM")]  # DHW equivalents reduced by some shading mechanism
    seed_years::Int64 = param_set[At("seed_years")]  # number of years to seed
    shade_years::Int64 = param_set[At("shade_years")]  # number of years to shade
    fog_years::Int64 = param_set[At("fog_years")]  # number of years to fog

    loc_k_area::Matrix{Float64} = cache.site_k_area
    fec_params_per_m²::Matrix{Float64} = _to_group_size(domain.coral_growth, corals.fecundity) # number of larvae produced per m²

    # Caches
    conn = domain.conn

    # Determine contribution of each source to a sink location
    # i.e., columns should sum to 1!
    TP_data = conn ./ sum(conn; dims=1)

    # sf = cache.sf  # unused as it is currently deactivated
    fec_all = cache.fec_all
    fec_scope = cache.fec_scope
    recruitment = cache.recruitment
    dhw_t = cache.dhw_step
    C_t = cache.cov_tmp
    depth_coeff = cache.depth_coeff

    site_data = domain.site_data
    depth_coeff .= depth_coefficient.(site_data.depth_med)

    # Coral cover relative to available area (i.e., 1.0 == site is filled to max capacity)
    C_cover::Array{Float64,4} = zeros(tf, n_groups, n_sizes, n_locs)
    C_cover[1, :, :, :] .= _group_cover_locs(domain.coral_growth, domain.init_coral_cover)
    cover_tmp = zeros(n_locs)

    # Locations that can support corals
    valid_locs::BitVector = location_k(domain) .> 0.0

    # Avoid placing importance on sites that were not considered
    # Lower values are higher importance/ranks.
    # Values of n_locs+1 indicate locations that were not considered in rankings.
    log_location_ranks = ZeroDataCube(;     # log seeding/fogging ranks
        T=Float64,
        timesteps=1:tf,
        locations=domain.site_ids,
        intervention=interventions()
    )

    Yshade = SparseArray(spzeros(tf, n_locs))
    Yfog = SparseArray(spzeros(tf, n_locs))
    Yseed = SparseArray(zeros(tf, 3, n_locs))  # 3 = the number of seeded coral types

    # Prep scenario-specific flags/values
    # Intervention strategy: < 0 is no intervention, 0 is random location selection, > 0 is guided
    is_guided = param_set[At("guided")] > 0
    if is_guided
        MCDA_approach = mcda_methods()[Int64(param_set[At("guided")])]
    end

    # Decisions should place more weight on environmental conditions
    # closer to the decision point
    α = 0.99
    decay = α .^ (1:Int64(param_set[At("plan_horizon")])+1) .^ 2

    # Years at which intervention locations are re-evaluated and deployed
    seed_decision_years = decision_frequency(seed_start_year, tf, seed_years, param_set[At("seed_deployment_freq")])
    fog_decision_years = decision_frequency(fog_start_year, tf, fog_years, param_set[At("fog_deployment_freq")])
    shade_decision_years = decision_frequency(shade_start_year, tf, shade_years, param_set[At("shade_deployment_freq")])
    # Define taxa and size class to seed, and identify their factor names
    taxa_to_seed = [2, 3, 5]
    target_class_id::BitArray = corals.class_id .== 2  # seed second smallest size class
    taxa_names::Vector{String} = collect(param_set.factors[occursin.("N_seed_", param_set.factors)])

    # Identify taxa and size class to be seeded
    seed_sc = _to_group_size(domain.coral_growth, (corals.taxa_id .∈ [taxa_to_seed]) .& target_class_id)

    # Extract colony areas for sites selected in m^2 and add adaptation values
    colony_areas = _to_group_size(
        domain.coral_growth, colony_mean_area(corals.mean_colony_diameter_m)
    )
    seeded_area = colony_areas[seed_sc] .* param_set[At(taxa_names)]

    # Set up assisted adaptation values
    a_adapt = zeros(n_groups, n_sizes)
    a_adapt[seed_sc] .= param_set[At("a_adapt")]

    # Flag indicating whether to seed or not to seed when unguided
    is_unguided = param_set[At("guided")] == 0.0
    seeding = any(param_set[At(taxa_names)] .> 0.0)
    apply_seeding = is_unguided && seeding

    # Flag indicating whether to fog or not fog
    apply_fogging = is_unguided && (fogging > 0.0)
    # Flag indicating whether to apply shading
    apply_shading = srm > 0.0

    # Calculate total area to seed respecting tolerance for minimum available space to still
    # seed at a site
    area_to_seed = sum(seeded_area)

    depth_criteria = identify_within_depth_bounds(
        site_data.depth_med, param_set[At("depth_min")], param_set[At("depth_offset")]
    )

    coral_habitable_locs = site_data.k .> 0.0
    if is_guided
        seed_pref = SeedPreferences(domain, param_set)
        fog_pref = FogPreferences(domain, param_set)

        # Create shared decision matrix, setting criteria values that do not change
        # between time steps
        decision_mat = decision_matrix(
            domain.site_ids,
            seed_pref.names;
            depth=site_data.depth_med
        )

        # Unsure what to do with this because it is usually empty
        # decision_mat[criteria=At("seed_zone")]

        # Remove locations that cannot support corals or are out of depth bounds
        # from consideration
        _valid_locs = coral_habitable_locs .& depth_criteria
        decision_mat = decision_mat[_valid_locs, :]

        # Number of time steps in environmental layers to look ahead when making decisions
        plan_horizon::Int64 = Int64(param_set[At("plan_horizon")])
    end

    # Set up distributions for natural adaptation/heritability
    c_mean_t_1::Array{Float64, 3} = repeat(_to_group_size(
        domain.coral_growth, corals.dist_mean
    ), 1, 1, n_locs)
    c_std::Array{Float64, 2} = _to_group_size(
        domain.coral_growth, corals.dist_std
    )

    c_mean_t = copy(c_mean_t_1)

    # Log of distributions
    dhw_tol_mean_log = cache.dhw_tol_mean_log  # tmp log for mean dhw tolerances

    # Cache for proportional mortality and coral population increases
    bleaching_mort = zeros(tf, n_groups, n_sizes, n_locs)

    #### End coral constants

    ## Update ecological parameters based on intervention option

    # Treat as enhancement from mean of "natural" DHW tolerance
    a_adapt[a_adapt.>0.0] .+= _to_group_size(
        domain.coral_growth, corals.dist_mean
    )[a_adapt.>0.0]

    # Pre-calculate proportion of survivers from wave stress
    # Sw_t = wave_damage!(cache.wave_damage, wave_scen, corals.wavemort90, n_species)

    p.r .= _to_group_size(domain.coral_growth, corals.growth_rate)
    p.mb .= _to_group_size(domain.coral_growth, corals.mb_rate)

    area_weighted_conn = conn .* site_k_area(domain)
    conn_cache = similar(area_weighted_conn.data)

    # basal_area_per_settler is the area in m^2 of a size class one coral
    basal_area_per_settler = colony_mean_area(
        corals.mean_colony_diameter_m[corals.class_id.==1]
    )

    # Dummy vars to fill/replace with ranks of selected locations
    selected_seed_ranks = []
    selected_fog_ranks = []

    # Cache matrix to store potential settlers
    potential_settlers = zeros(size(fec_scope)...)
    C_bins::Matrix{Float64} = hcat(
        zeros(n_groups),
        _to_group_size(domain.coral_growth, corals.bin_ub)
    )

    cover_blocks::Vector{Matrix{CoverBlock}} = [
        DynamicCoralCoverModel.blocks_model.CoverBlock.(
            C_cover[1, :, :, loc] .* (site_data.area[loc] .* site_data.k[loc]),
            C_bins[:, 1:end-1],
            C_bins[:, 2:end]
        ) for loc in 1:n_locs
    ]

    linear_extension = _to_group_size(domain.coral_growth, corals.linear_extension)
    survival_rate = 1.0 .- _to_group_size(domain.coral_growth, corals.mb_rate)

    size_classes::Vector{Matrix{SizeClass}} = [
        SizeClass.(
            cover_blocks[loc],
            repeat(1:n_sizes, 1, n_groups)',
            linear_extension,
            survival_rate
        ) for loc in 1:n_locs
    ]

    # Preallocate memory for temporaries
    temp_change = ones(n_groups, n_sizes, n_locs)
    C_t = zeros(n_groups, n_sizes, n_locs)
    cover_copy = zeros(n_groups, n_sizes, n_locs)

    for tstep::Int64 in 2:tf
        change_view = [@view temp_change[:, :, loc] for loc in 1:n_locs]
        DynamicCoralCoverModel.blocks_model.apply_changes!.(size_classes, change_view)
        recruitment .*= reshape(temp_change[:, 1, :], (n_groups, n_locs))

        C_t .= C_cover[tstep-1, :, :, :] .* reshape(
            site_data.area .* site_data.k, (1, 1, n_locs)
        )

        for i in 1:n_locs
            C_t[:, :, i] .= DynamicCoralCoverModel.blocks_model.timestep(
                C_t[:, :, i],
                recruitment[:, i],
                size_classes[i],
                site_data.k[i] * site_data.area[i],
                tstep,
                false
            )
        end

        C_t ./= reshape(site_data.area .* site_data.k, (1, 1, n_locs))
        replace!(C_t, NaN=>0.0)
        cover_copy .= copy(C_t)

        # Check if size classes are inappropriately out-growing available space
        proportional_adjustment!(
            @view(C_t[:, :, valid_locs]),
            cover_tmp[valid_locs]
        )

        # Update initial condition
        C_cover[tstep, :, :, valid_locs] .= C_t[:, :, valid_locs]

        if tstep <= tf
            # Natural adaptation
            adjust_DHW_distribution!(
                @view(C_cover[(tstep-1), :, :, :]), c_mean_t, p.r
            )

            # Set values for t to t-1
            c_mean_t_1 .= c_mean_t

            if in_debug_mode
                # Log dhw tolerances if in debug mode
                dhw_tol_mean_log[tstep, :, :] .= mean.(c_mean_t)
            end
        end

        # Calculates scope for coral fedundity for each size class and at each location
        fecundity_scope!(fec_scope, fec_all, fec_params_per_m², C_t, loc_k_area)

        loc_coral_cover = dropdims(sum(C_t; dims=(1, 2)), dims=1)  # dims: 1 * nsites
        leftover_space_m² = relative_leftover_space(loc_coral_cover) .* loc_k_area

        # Reset potential settlers to zero
        potential_settlers .= 0.0
        recruitment .= 0.0

        # Recruitment represents additional cover, relative to total site area
        # Recruitment/settlement occurs after the full moon in October/November
        recruitment[:, valid_locs] .=
            settler_cover(
                fec_scope,
                conn,
                leftover_space_m²,
                sim_params.max_settler_density,
                sim_params.max_larval_density,
                basal_area_per_settler,
                potential_settlers,
            )[
                :, valid_locs
            ] ./ loc_k_area[:, valid_locs]

        settler_DHW_tolerance!(
            c_mean_t_1,
            c_mean_t,
            site_k_area(domain),
            TP_data,  # ! IMPORTANT: Pass in transition probability matrix, not connectivity!
            recruitment,
            fec_params_per_m²,
            param_set[At("heritability")],
        )

        # Add recruits to current cover
        C_t[:, 1, :] .+= recruitment

        # Cover copy needs to include recruits so overall mortality can be calculated to
        # apply to cover blocks
        cover_copy[:, 1, :] .= C_t[:, 1, :]

        # Update available space
        loc_coral_cover = dropdims(sum(C_t, dims=(1, 2)), dims=1)  # dims: 1 * nsites
        leftover_space_m² = relative_leftover_space(loc_coral_cover) .* loc_k_area

        # Determine intervention locations whose deployment is assumed to occur
        # between November to February.
        #
        # Nominal sequence of events is conceptualised as:
        # - Corals spawn and are recruited
        # - SRM is applied
        # - Fogging is applied next
        # - Seeding interventions occur
        # - bleaching then occurs
        # - then cyclones hit
        #
        # Seeding occurs before bleaching as current tech only allows deployments shortly
        # after spawning. If bio-banking or similar tech comes up to speed then we could
        # potentially deploy at alternate times.

        # Apply regional cooling effect before selecting locations to seed
        dhw_t .= dhw_scen[tstep, :]  # subset of DHW for given timestep
        if apply_shading && shade_decision_years[tstep]
            Yshade[tstep, :] .= srm

            # Apply reduction in DHW due to SRM
            dhw_t .= max.(0.0, dhw_t .- srm)
        end

        if is_guided
            if fog_decision_years[tstep] && (fogging .> 0.0)
                selected_fog_ranks = select_locations(
                    fog_pref,
                    decision_mat,
                    MCDA_approach,
                    min_iv_locs
                )

                if !isempty(selected_fog_ranks)
                    log_location_ranks[tstep, At(selected_fog_ranks), At(:fog)] .= 1:length(selected_fog_ranks)
                end
            end
        elseif apply_fogging && fog_decision_years[tstep]
            selected_fog_ranks = unguided_selection(
                domain.site_ids,
                min_iv_locs,
                vec(leftover_space_m²),
                depth_criteria
            )

            log_location_ranks[tstep, At(selected_fog_ranks), At(:fog)] .= 1.0
        end

        has_fog_locs::Bool = !isempty(selected_fog_ranks)

        # Fog selected locations
        if fog_decision_years[tstep] && has_fog_locs
            fog_locs = findall(log_location_ranks.locations .∈ [selected_fog_ranks])
            fog_locations!(@view(Yfog[tstep, :]), fog_locs, dhw_t, fogging)
        end

        if is_guided && seed_decision_years[tstep]
            # Use modified projected DHW (may have been affected by fogging or shading)
            dhw_p = copy_datacube(dhw_scen)
            dhw_p[tstep, :] .= dhw_t

            dhw_projection = weighted_projection(dhw_p, tstep, plan_horizon, decay, tf)
            wave_projection = weighted_projection(wave_scen, tstep, plan_horizon, decay, tf)

            # Determine connectivity strength weighting by area.
            # Accounts for strength of connectivity where there is low/no coral cover
            in_conn, out_conn, _ = connectivity_strength(area_weighted_conn, vec(loc_coral_cover), conn_cache)

            update_criteria_values!(
                decision_mat;
                heat_stress=dhw_projection[_valid_locs],
                wave_stress=wave_projection[_valid_locs],
                coral_cover=loc_coral_cover[_valid_locs],  # Coral cover relative to `k`
                in_connectivity=in_conn[_valid_locs],  # area weighted connectivities for time `t`
                out_connectivity=out_conn[_valid_locs]
            )

            # IDs of valid locations considering locations that have space for corals
            locs_with_space = vec(leftover_space_m²) .> 0.0
            considered_locs = findall(_valid_locs .& locs_with_space)

            # IDs of valid locations considering locations that have space for corals
            locs_with_space = vec(leftover_space_m²) .> 0.0
            considered_locs = findall(_valid_locs .& locs_with_space)

            selected_seed_ranks = select_locations(
                seed_pref,
                decision_mat[location=locs_with_space[_valid_locs]],
                MCDA_approach,
                site_data.cluster_id,
                area_to_seed,
                considered_locs,
                vec(leftover_space_m²),
                min_iv_locs,
                max_members
            )

            # Log rankings as appropriate
            if !isempty(selected_seed_ranks)
                log_location_ranks[tstep, At(selected_seed_ranks), At(:seed)] .= 1:length(selected_seed_ranks)
            end
        elseif apply_seeding && seed_decision_years[tstep]
            # Unguided deployment, seed/fog corals anywhere, so long as available space > 0
            selected_seed_ranks = unguided_selection(
                domain.site_ids,
                min_iv_locs,
                vec(leftover_space_m²),
                depth_criteria
            )

            log_location_ranks[tstep, At(selected_seed_ranks), At(:seed)] .= 1.0

            # Estimate proportional change in cover to apply to cubes
        end

        # Check if locations are selected (can reuse previous selection)
        has_seed_locs::Bool = !isempty(selected_seed_ranks)

        # Apply seeding (assumed to occur after spawning)
        if seed_decision_years[tstep] && has_seed_locs
            # Seed selected locations
            seed_locs = findall(log_location_ranks.locations .∈ [selected_seed_ranks])
            seed_corals!(
                C_t,
                vec(loc_k_area),
                vec(leftover_space_m²),
                seed_locs,  # use location indices
                seeded_area,
                seed_sc,
                a_adapt,
                @view(Yseed[tstep, :, :]),
                c_std,
                c_mean_t,
            )

            # Add coral seeding to cover copy and recruitment
            recruitment += C_t[:, 1, :] .- cover_copy[:, 1, :]
            cover_copy[:, 1, :] .= C_t[:, 1, :]
        end

        # Calculate and apply bleaching mortality
        # Bleaching typically occurs in the warmer months (November - February)
        #    This: `dhw_t .* (1.0 .- wave_scen[tstep, :])`
        #    attempts to account for the cooling effect of storms / high wave activity
        # `wave_scen` is normalized to the maximum value found for the given wave scenario
        # so what causes 100% mortality can differ between runs.
        bleaching_mortality!(
            C_t,
            dhw_t,  # collect(dhw_t .* (1.0 .- @view(wave_scen[tstep, :]))),
            depth_coeff,
            c_std,
            c_mean_t_1,
            c_mean_t,
            @view(bleaching_mort[(tstep-1):tstep, :, :, :])
        )

        # Coral deaths due to selected cyclone scenario
        # Peak cyclone period is January to March
        # TODO: Update cyclone data to hold data for relevant functional groups
        cyclone_mortality!(@views(C_t), p, cyclone_mortality_scen[tstep, :, :]')

        # Update record
        C_cover[tstep, :, :] .= C_t
        cover_copy[cover_copy .== 0] .= 1.0
        temp_change = _group_cover_locs(
            domain.coral_growth,
            C_cover[tstep, :, :]
        ) ./ cover_copy
    end

    # Could collate critical DHW threshold log for corals to reduce disk space...
    # dhw_tol_mean = dropdims(mean(dhw_tol_mean_log, dims=3), dims=3)
    # dhw_tol_mean_std = dropdims(mean(dhw_tol_std_log, dims=3), dims=3)
    # collated_dhw_tol_log = DataCube(cat(dhw_tol_mean, dhw_tol_mean_std, dims=3);
    #     timesteps=1:tf, species=corals.coral_id, stat=[:mean, :stdev])
    if in_debug_mode
        collated_dhw_tol_log = DataCube(
            dhw_tol_mean_log; timesteps=1:tf, species=corals.coral_id, sites=1:n_locs
        )
    else
        collated_dhw_tol_log = false
    end

    # Set variables to nothing so garbage collector clears them
    # Leads to memory leak issues in multiprocessing contexts without these.
    wave_scen = nothing
    dhw_tol_mean_log = nothing

    return (
        raw=C_cover,
        seed_log=Yseed,
        fog_log=Yfog,
        shade_log=Yshade,
        site_ranks=log_location_ranks,
        bleaching_mortality=bleaching_mort,
        coral_dhw_log=collated_dhw_tol_log,
    )
end


"""
    cyclone_mortality!(coral_cover, coral_params, cyclone_mortality)::Nothing

Apply cyclone mortalities.

# Arguments
- `coral_cover` : Coral cover for current time step
- `coral_params` : Coral parameters indicating indices of small/mid/large size classes
- `cyclone_mortality` : Mortalities for each functional group and size class
"""
function cyclone_mortality!(coral_cover, coral_params, cyclone_mortality)::Nothing
    # TODO: Move to own file.

    # Small class coral mortality
    coral_deaths_small = coral_cover[coral_params.small, :] .* cyclone_mortality
    coral_cover[coral_params.small, :] -= coral_deaths_small

    # Mid class coral mortality
    coral_mid = hcat(collect(Iterators.partition(coral_params.mid, length(coral_params.small)))...)
    for i in size(coral_mid, 1)
        coral_deaths_mid = coral_cover[coral_mid[i, :], :] .* cyclone_mortality
        coral_cover[coral_mid[i, :], :] -= coral_deaths_mid
    end

    # Large class coral mortality
    coral_deaths_large = coral_cover[coral_params.large, :] .* cyclone_mortality
    coral_cover[coral_params.large, :] -= coral_deaths_large

    # Ensure no negative values
    clamp!(coral_cover, 0.0, 1.0)

    return nothing
end
