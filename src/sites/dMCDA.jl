"""Objects and methods for Dynamic Multi-Criteria Decision Analysis/Making"""

using StatsBase
using Distances
using Combinatorics
using JMcDM
using ADRIA: order_ranking, mcda_vikor, mcda_topsis


global mcda_methods = [
    order_ranking,
    mcda_vikor,
    mcda_topsis,
    [ArasMethod(), true],
    [CocosoMethod(), true],
    [CodasMethod(), true],
    [CoprasMethod(), false],
    [EdasMethod(), true],
    [GreyMethod(), true],
    [MabacMethod(), true],
    [MaircaMethod(), false],
    [MarcosMethod(), true],
    [MooraMethod(), true],
    #[MoosraMethod(), true],
    [PIVMethod(), true],
    [PSIMethod(), true],
    [ROVMethod(), true],
    [SawMethod(), true],
    [TopsisMethod(), true],
    [VikorMethod(), false],
    [WPMMethod(), true],
    [WaspasMethod(), true]
]

"""
    create_criteria_store(site_ids::AbstractArray, criteria::NamedTuple) 

    Constructs the criteria NamedDimsArray for performing site selection. 
    This is used to construct decision matrices for the mcda methods.

# Arguments
- `site_ids` : site ids as integers.
- `criteria` : NamedTuple of vectors of length nsites containing criteria values to 
be used to construct mcda matrices. Keys should correspond to weight names and begin with "iv__".
Can also be entered as named varargs.
E.g. the criteria for heat stress will be iv__heat_stress and it's weight will be iv__heat_stress__seed_shade,
indicating it is used for seeding and shading.
"""
function create_criteria_store(site_ids::AbstractArray; criteria...)
    return create_criteria_store(site_ids, criteria)
end
function create_criteria_store(site_ids::AbstractArray, criteria::NamedTuple)
    criteria_matrix = zeros(length(site_ids), length(criteria))

    for (ind, crit_key) in enumerate(keys(criteria))
        criteria_matrix[:, ind] .= criteria[crit_key][site_ids]
    end
    return NamedDimsArray(criteria_matrix, locations=site_ids, criteria=collect(keys(criteria)))

end

"""
    create_tolerances_store(tolerances::NamedTuple)  

    Constructs the tolerances NamedDimsArray which is used to filter the decision matrix.

# Arguments
- `tolerances` : NamedTuple with format (criteria_name1=(operation,value),criteria_name2=(operation,value),...)
where operation is < or > and value is the tolerance value.

"""
function create_tolerances_store(tolerances::NamedTuple)
    tol_store = [x -> tolerances[tol_key][1](x, tolerances[tol_key][2]) for tol_key in keys(tolerances)]
    return NamedDimsArray(tol_store, criteria=collect(keys(tolerances)))
end

"""
    mcda_normalize(x::Vector)::Vector

Normalize a Vector (wse/wsh) for MCDA.
"""
function mcda_normalize(x::Vector)::Vector
    return x ./ sum(x)
end

"""
    mcda_normalize(x::Matrix)::Matrix

Normalize a Matrix (SE/SH) for MCDA.
"""
function mcda_normalize(x::Matrix)::Matrix
    return x ./ sqrt.(sum(x .^ 2, dims=1))
end


"""
    align_rankings!(rankings::Array, s_order::Matrix, col::Int64)::Nothing

Align a vector of site rankings to match the indicated order in `s_order`.
"""
function align_rankings!(rankings::Array, s_order::Matrix, col::Int64)::Nothing
    # Fill target ranking column
    for (i, site_id) in enumerate(s_order[:, 1])
        rankings[rankings[:, 1].==site_id, col] .= s_order[i, 3]
    end

    return
end

"""
    rank_sites!(S, weights, rankings, n_site_int, rank_col)
    rank_seed_sites!(S, weights, rankings, n_site_int)
    rank_shade_sites!(S, weights, rankings, n_site_int)

# Arguments
- `S` : Matrix, Site preference values
- `weights` : weights to apply
- `rankings` : vector of site ranks to update
- `n_site_int` : number of sites to select for interventions
- `rank_col` : column to fill with rankings (2 for seed, 3 for shade)

# Returns
- `prefsites` : sites in order of their rankings
"""
function rank_sites!(S, weights, rankings, n_site_int, site_ids, mcda_func, rank_col)::Tuple{Vector{Int64},Matrix{Union{Float64,Int64}}}
    # Filter out all non-preferred sites
    selector = vec(.!all(S .== 0, dims=1))

    # weights in order of: in_conn, out_conn, wave, heat, predecessors, low cover
    weights = weights[selector]
    S = S[:, selector]

    s_order = retrieve_ranks(S, weights, mcda_func, site_ids)

    last_idx = min(n_site_int, size(s_order, 1))
    prefsites = Int.(s_order[1:last_idx, 1])

    # Match by site_id and assign rankings to log
    align_rankings!(rankings, s_order, rank_col)

    return prefsites, s_order
end

function rank_seed_sites!(S, weights, rankings, n_site_int, site_ids, mcda_func)::Tuple{Vector{Int64},Matrix{Union{Float64,Int64}}}
    rank_sites!(S, weights, rankings, n_site_int, site_ids, mcda_func, 2)
end
function rank_shade_sites!(S, weights, rankings, n_site_int, site_ids, mcda_func)::Tuple{Vector{Int64},Matrix{Union{Float64,Int64}}}
    rank_sites!(S, weights, rankings, n_site_int, site_ids, mcda_func, 3)
end

"""
    filter_decision_matrix(criteria_store::NamedDimsArray, tolerances::NamedDimsArray)  

# Arguments
- `criteria_store` : contains criteria in each column for sites in each row.
- `tolerances` : contains thresholds for specified criteria, with keys matching those in criteria_store.
                Each key specifies the operation function >tolerance or <tolerance.

# Returns
- `criteria_store` : filtered version of criteria_store input.
"""
function filter_decision_matrix(criteria_store::NamedDimsArray, tolerances::NamedDimsArray)

    rule = sum(map.(tolerances.data, criteria_store(tolerances.criteria)'), dims=:criteria) .== length(tolerances.criteria)
    criteria_store = criteria_store[rule[criteria=1], :]

    return criteria_store
end


"""
    create_intervention_matrix(criteria_store::NamedDimsArray, params::NamedDimsArray, int_type::String)
  
# Arguments
- `criteria_store` : pre-filtered criteria to be used in the decision matrix for mcda.
- `params` : scenario parameters, must include weights for criteria for the intervention of interest.
- `int_type` : intervention type indicated as a string. Used to find relevant weights in `params`.

# Returns
- `S` : decision matrix for use in mcda technique of choice.
- 'ws': weights for each of the criteria in the decision matrix.

"""
function create_intervention_matrix(criteria_store::NamedDimsArray, params::NamedDimsArray, int_type::String)
    # Define intervention decision matrix
    int_params = params.factors[occursin.("iv__", params.factors).&occursin.(int_type, params.factors)]
    crit_inds = [findall.([occursin.(String(crit_name), int_params) for crit_name in criteria_store.criteria])...;]

    ws = mcda_normalize(Array(params(int_params)))
    S = Matrix(criteria_store[criteria=crit_inds])
    return S, ws
end


"""
    guided_site_selection(d_vars::DMCDA_vars, criteria_df::DataFrame, alg_ind::Int64, log_seed::Bool, log_shade::Bool, prefseedsites::AbstractArray{Int64}, prefshadesites::AbstractArray{Int64}, rankingsin::Matrix{Int64})

# Arguments
- `d_vars` : DMCDA_vars type struct containing weightings and criteria values for site selection.
- `alg_ind` : integer indicating MCDA aggregation method to use (0: none, 1: order ranking, 2:topsis, 3: vikor)
- `log_seed` : boolean indicating whether seeding sites are being re-assesed at current time
- `log_shade` : boolean indicating whether shading/fogging sites are being re-assesed at current time
- `prefshadesites` : previous time step's selection of sites for shading
- `prefseedsites` : previous time step's selection of sites for seeding
- `rankingsin` : pre-allocated store for site rankings
- `in_conn` : in-degree centrality
- `out_conn` : out-degree centrality
- `strong_pred` : strongest predecessors

# Returns
Tuple :
    - `prefseedsites` : n_site_int highest ranked seeding sites
    - `prefshadesites` : n_site_int highest ranked shading/fogging sites
    - `rankings` : n_sites ⋅ 3 matrix holding [site_id, seeding_rank, shading_rank],
        Values of 0 indicate sites that were not considered
"""
function guided_site_selection(criteria_store::NamedDimsArray,
    params::NamedDimsArray, thresholds::NamedDimsArray, n_site_int::Int64,
    distances::Matrix, minimum_distance::Float64, log_seed::B, log_shade::B,
    prefseedsites::IA, prefshadesites::IA,
    rankingsin::Matrix{T}
)::Tuple where {T<:Int64,IA<:AbstractArray{<:Int64},IB<:AbstractArray{<:Int64},B<:Bool}

    site_ids::Array{Int64} = criteria_store.reefs
    use_dist::Int64 = d_vars.use_dist
    min_dist::Float64 = d_vars.min_dist
    site_ids = copy(d_vars.site_ids)
    n_sites::Int64 = length(site_ids)

    # if no sites are available, abort
    if n_sites == 0
        return zeros(Int64, length(prefseedsites)), zeros(Int64, length(prefshadesites)), rankingsin
    end

    # site_id, seeding rank, shading rank
    rankings = Int64[site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]
    mcda_func = mcda_methods[Int(params("guided"))]

    criteria_store = filter_decision_matrix(criteria_store, thresholds)

    if isempty(criteria_store)
        # if all rows have nans and A is empty, abort mission
        return zeros(Int64, length(prefseedsites)), zeros(Int64, length(prefshadesites)), rankingsin
    end

    # cap to number of sites left after risk filtration
    n_site_int = min(n_site_int, length(criteria_store.locations))

    # if seeding, create seeding specific decision matrix
    if log_seed
        SE, wse = create_intervention_matrix(criteria_store, params, "seed")
    end

    # if shading, create shading specific decision matrix
    if log_shade
        SH, wsh = create_intervention_matrix(criteria_store, params, "shade")
    end

    if log_seed && isempty(SE)
        prefseedsites = zeros(Int64, n_site_int)
    elseif log_seed
        prefseedsites, s_order_seed = rank_seed_sites!(SE, wse, rankings, n_site_int, criteria_store.locations, mcda_func)
        if use_dist != 0
            prefseedsites, rankings = distance_sorting(prefseedsites, s_order_seed, distances, minimum_distance, rankings, 2)
        end
    end

    if log_shade && isempty(SH)
        prefshadesites = zeros(Int64, n_site_int)
    elseif log_shade
        prefshadesites, s_order_shade = rank_shade_sites!(SH, wsh, rankings, n_site_int, criteria_store.locations, mcda_func)

        if use_dist != 0
            prefshadesites, rankings = distance_sorting(prefshadesites, s_order_shade, distances, minimum_distance, rankings, 3)
        end
    end

    # Replace with input rankings if seeding or shading rankings have not been filled
    if sum(prefseedsites) == 0
        rankings[:, 2] .= rankingsin[:, 2]
    end

    if sum(prefshadesites) == 0
        rankings[:, 3] .= rankingsin[:, 3]
    end

    return prefseedsites, prefshadesites, rankings
end


"""
    distance_sorting(pref_locations::AbstractArray{Int}, location_order::AbstractVector, dist::Array{Float64}, dist_thresh::Float64)::AbstractArray{Int}

Find selected locations with distances between each other < dist_thresh*(median distance).
Replaces these locations with locations in the top_n ranks if the distance between these locations is greater.

# Arguments
- `pref_locations` : original n highest ranked locations selected for seeding or shading.
- `location_order` : current order of ranked locations in terms of numerical site ID.
- `dist` : Matrix of unique distances between locations.
- `min_dist` : minimum distance between locations for selected locations.

# Returns
- `rep_locations` : new set of selected locations for seeding or shading.
"""
function distance_sorting(pref_locations::AbstractArray{Int}, l_order::Matrix{Union{Float64,Int64}}, dist::Array{Float64},
    min_dist::Float64, rankings::Matrix{Int64}, rank_col::Int64)::Tuple{Vector{Union{Float64,Int64}},Matrix{Int64}}
    # set-up
    n_locations = length(pref_locations)
    location_order = l_order[:, 1]

    # locations to select alternatives from
    alt_locations = setdiff(location_order, pref_locations)[1:length(location_order)-n_locations]

    # find all selected locations closer than the min distance
    pref_dists = findall(dist[pref_locations, pref_locations] .< min_dist)
    # indices to replace
    inds_rep = sort(unique(reinterpret(Int64, pref_dists)))
    # number of locations to replace
    select_n = length(inds_rep)
    # indices to keep
    inds_keep = collect(1:length(pref_locations))
    inds_keep = setdiff(inds_keep, inds_rep)

    # storage for new set of locations
    rep_locations = pref_locations

    while (length(alt_locations) .>= select_n)
        rep_locations = [rep_locations[inds_keep[:]]; alt_locations[1:select_n]]

        # Find all locations within these highly ranked but unselected locations which are further apart
        alt_dists = dist[rep_locations, rep_locations] .> min_dist

        # Select from these locations those far enough away from all locations
        inds_keep = sum(alt_dists, dims=2) .== n_locations - 1

        # Keep locations that were far enough away last iteration
        inds_keep[1:end-select_n] .= true
        if length(inds_keep) == n_locations
            select_n = 0
            break
        else
            # remove checked alt_locations
            alt_locations = setdiff(alt_locations, alt_locations[1:select_n])
            select_n = sum(.!inds_keep)
        end
    end

    # If not all locations could be replaced, just use highest ranked remaining pref_locations
    if (select_n != 0) && !isempty(setdiff(pref_locations, rep_locations))
        rem_pref_locations = setdiff(pref_locations, rep_locations)
        rep_locations[end-select_n+1:end] .= rem_pref_locations[1:select_n]
    end

    new_location_order = setdiff(location_order, rep_locations)
    new_location_order = [rep_locations; new_location_order]
    l_order[:, 1] .= new_location_order
    # Match by site_id and assign rankings to log

    align_rankings!(rankings, l_order, rank_col)

    return rep_locations, rankings
end


"""
    retrieve_ranks(S::Matrix, weights::Array{Float64}, mcda_func::Function, site_ids::Array{Int64})

Get location ranks using mcda technique specified in mcda_func, weights and a decision matrix S.

# Arguments
- `S` : decision matrix containing criteria values for each location (n locations)*(m criteria)
- `weights` : importance weights for each criteria. 
- `mcda_func` : function to use for mcda, specified as an element from mcda_methods.
- `site_ids` : array of integers indicating site ids still remaining after filtering.

# Returns
- `s_order` : [site_ids, criteria values, ranks]
"""
function retrieve_ranks(S::Matrix, weights::Array{Float64}, mcda_func::Function, site_ids::Array{Int64})
    S = mcda_normalize(S)
    S .= S .* repeat(weights', size(S, 1), 1)
    scores = mcda_func(S)

    return retrieve_ranks(S, scores, true, site_ids)
end
function retrieve_ranks(S::Matrix, weights::Array, mcda_func::Vector, site_ids::Array{Int64})
    fns = repeat([maximum], length(weights))
    results = mcdm(MCDMSetting(S, weights, fns), mcda_func[1])

    return retrieve_ranks(S, results.scores, mcda_func[2], site_ids)
end
function retrieve_ranks(S::Matrix, scores::Array, rev_val::Bool, site_ids::Array{Int64})
    s_order = Union{Float64,Int64}[Int.(site_ids) scores 1:size(S, 1)]
    s_order .= sortslices(s_order, dims=1, by=x -> x[2], rev=rev_val)
    @views s_order[:, 3] .= Int.(1:size(S, 1))

    return s_order
end


"""
    site_selection(criteria_store::NamedDimsArray, scenario::NamedDimsArray, tolerances::NamedTuple,
        site_ids::AbstractArray, site_distances::Matrix, med_site_distance::Float64, n_site_int::Int64)

Perform site selection using a set of criteria, tolerances, locations and location distances.
# Arguments
- `criteria_store` : contains criteria for a single location selection instance.
- `scenario` : contains parameters for a single location selection instance, including tolerance values 
and parameters for distance sorting.
- `tolerances` : specifies criteria tolerances and has keys with names corresponding to criteria in criteria_store. 
For example, (iv__heat_stress=(<,0.5)), implies the criteria "iv__heat_stress" must not have values greater than 0.5.
- `location_ids` : array of length nlocations containing indices of locations to be selected from.
- `site_distances` : Matrix of distances between locations.
- `med_site_distance` : Median distance between locations in the site_distances matrix.
- `n_site_int` : number of locations to select to perform intervention.

# Returns
- `ranks` : n_reps * sites * 3 (last dimension indicates: site_id, seeding rank, shading rank)
    containing ranks for single scenario.
"""
function site_selection(criteria_store::NamedDimsArray, scenario::NamedDimsArray, tolerances::NamedTuple,
    site_ids::AbstractArray, site_distances::Matrix, med_site_distance::Float64, n_site_int::Int64)

    tolerances_store = create_tolerances_store(tolerances)
    min_distance = med_site_distance .* scenario("dist_thresh")
    n_sites = length(site_ids)

    # site_id, seeding rank, shading rank
    rankingsin = [site_ids zeros(Int64, (n_sites, 1)) zeros(Int64, (n_sites, 1))]
    prefseedsites = zeros(Int64, (1, n_site_int))
    prefshadesites = zeros(Int64, (1, n_site_int))

    (_, _, ranks) = guided_site_selection(criteria_store, scenario, tolerances_store, n_site_int,
        site_distances, min_distance, true, true, prefseedsites, prefshadesites, rankingsin)
    return ranks
end

"""
    run_site_selection(domain::Domain, scenarios::DataFrame, tolerances::NamedTuple, coral_covers::NamedDimsArray;
        aggregation_method=nothing, target_seed_sites=nothing, target_shade_sites=nothing)

Perform site selection for a given domain for multiple scenarios defined in a dataframe.

# Arguments
- `domain` : ADRIA Domain type, indicating geographical domain to perform site selection over.
- `scenarios` : DataFrame of criteria weightings and thresholds for each scenario.
- `tolerances` : NamedTuple specifying tolerances for pre-selection filtering. E.g. `tolerances = (iv__coral_cover=(>, x -> f_coral_cover(x)),
    iv__heat_stress=(<, x -> x), iv__wave_stress=(<, x -> x))`, where the keys correspond to criteria names in the Domain,
    the first element difines the filtering operation (> or <) and the second element defines the operation on the tolerance parameter
    (x->x is directly using parameter value). Tolerance values in `scenarios` DataFrame will have names "iv__"+(criteria name)+"__tol".
- `coral_covers` : contains coral covers for each selection scenario, size (N locations, M scenarios).
- `target_seed_sites` : optional additional set of sites to only consider during seeding (must be a subset of the location ids in site data in Domain).
- `target_shade_sites` : optional additional set of sites to only consider during shading (must be a subset of the location ids in site data in Domain).

# Returns
- `ranks_store` : number of scenarios * sites * 3 (last dimension indicates: site_id, seed rank, shade rank)
    containing ranks for each scenario run.
- `aggregated_ranks_store` : if aggregation method is selected, the aggregated ranks_store output.
"""
function run_site_selection(domain::Domain, scenarios::DataFrame, tolerances::NamedTuple, coral_covers::NamedDimsArray;
    target_seed_sites=nothing, target_shade_sites=nothing)
    ranks_store = NamedDimsArray(
        zeros(nrow(scenarios), length(domain.site_ids), 3),
        scenarios=1:nrow(scenarios),
        locations=domain.site_ids,
        ranks=["site_id", "seed_rank", "shade_rank"],
    )

    site_data = domain.site_data

    # Pre-calculate maximum depth to consider
    scenarios[:, "max_depth"] .= scenarios.depth_min .+ scenarios.depth_offset
    target_dhw_scens = unique(scenarios[:, "dhw_scenario"])
    target_wave_scens = unique(scenarios[:, "wave_scenario"])

    target_site_ids = Int64[]
    dhw_scens = domain.dhw_scens
    wave_scens = domain.wave_scens

    if !isnothing(target_seed_sites)
        append!(target_site_ids, target_seed_sites)
    end

    if !isnothing(target_shade_sites)
        append!(target_site_ids, target_shade_sites)
    end

    # Pre-calculate maximum depth to consider
    scenarios[:, "max_depth"] .= scenarios.depth_min .+ scenarios.depth_offset
    criteria_store = create_criteria_store(collect(1:length(domain.site_ids)), domain.mcda_criteria)

    coral_cover, coral_space = coral_cover_criteria(domain.site_data, coral_covers)

    in_connectivity = connectivity_criteria(domain.in_conn, coral_covers, domain.site_data.area)
    out_connectivity = connectivity_criteria(domain.out_conn, coral_covers, domain.site_data.area)

    criteria_store(:iv__wave_stress) .= env_stress_criteria(Array(dropdims(mean(wave_scens, dims=(:timesteps, :scenarios)) .+ var(wave_scens, dims=(:timesteps, :scenarios)), dims=:timesteps)))
    criteria_store(:iv__heat_stress) .= env_stress_criteria(Array(dropdims(mean(dhw_scens, dims=(:timesteps, :scenarios)) .+ var(dhw_scens, dims=(:timesteps, :scenarios)), dims=:timesteps)))

    for (cover_ind, scen) in enumerate(eachrow(scenarios))

        tol_temp = ()
        for tol = keys(tolerances)
            tol_temp = (; tol_temp..., tol => (tolerances[tol][1], map(tolerances[tol][2], scen[string(String(tol), "__tol")])))
        end

        depth_criteria = (domain.site_data.depth_med .<= scen.max_depth) .& (domain.site_data.depth_med .>= scen.depth_min)
        depth_priority = findall(depth_criteria)

        criteria_store(:iv__coral_cover) .= coral_cover[scenarios=cover_ind]
        criteria_store(:iv__coral_space) .= coral_space[scenarios=cover_ind]
        criteria_store(:iv__in_connectivity) .= in_connectivity[scenarios=cover_ind]
        criteria_store(:iv__out_connectivity) .= out_connectivity[scenarios=cover_ind]

        considered_sites = target_site_ids[findall(in(depth_priority), target_site_ids)]
        scen_set = NamedDimsArray(Vector(scen), factors=names(scen))

        ranks_store(scenarios=cover_ind, sites=domain.site_ids[depth_priority]) .= site_selection(
            criteria_store[locations=depth_priority],
            scen_set,
            tol_temp,
            depth_priority,
            domain.site_distances,
            domain.median_site_distance,
            domain.sim_constants.n_site_int
        )

    end

    return ranks_store
end

"""
    site_selection(domain::Domain, scenario::DataFrameRow{DataFrame,DataFrames.Index}, w_scens::NamedDimsArray, dhw_scens::NamedDimsArray, sum_cover::AbstractArray, area_to_seed::Float64)

Perform site selection using a chosen mcda aggregation method, domain, initial cover, criteria weightings and thresholds.

# Arguments
- `scenario` : contains criteria weightings and thresholds for a single scenario.
- `mcda_vars` : site selection criteria and weightings structure
- `w_scens` : array of length nsites containing wave scenario.
- `dhw_scens` : array of length nsites containing dhw scenario.
- `sum_cover` : summed cover (over species) for single scenario being run, for each site.
- `area_to_seed` : area of coral to be seeded at each time step in km^2

# Returns
- `ranks` : n_reps * sites * 3 (last dimension indicates: site_id, seeding rank, shading rank)
    containing ranks for single scenario.
"""
function site_selection(domain::Domain, scenario::DataFrameRow, w_scens::NamedDimsArray, dhw_scens::NamedDimsArray,
    site_ids::AbstractArray, sum_cover::AbstractArray, area_to_seed::Float64)::Matrix{Int64}

    mcda_vars = DMCDA_vars(domain, scenario, site_ids, sum_cover, area_to_seed, w_scens, dhw_scens)
    n_sites = length(mcda_vars.site_ids)

    # site_id, seeding rank, shading rank
    rankingsin = [mcda_vars.site_ids zeros(Int64, (n_sites, 1)) zeros(Int64, (n_sites, 1))]

    prefseedsites::Matrix{Int64} = zeros(Int64, (1, mcda_vars.n_site_int))
    prefshadesites::Matrix{Int64} = zeros(Int64, (1, mcda_vars.n_site_int))

    (_, _, ranks) = guided_site_selection(mcda_vars, scenario.guided, true, true, prefseedsites, prefshadesites, rankingsin)

    return ranks
end

"""
    unguided_site_selection(prefseedsites, prefshadesites, seed_years, shade_years, n_site_int, max_cover)

Randomly select seed/shade site locations for the given year, constraining to sites with max. carrying capacity > 0.
Here, `max_cover` represents the max. carrying capacity for each site (the `k` value).

# Arguments
- `prefseedsites` : Previously selected sites
- `seed_years` : bool, indicating whether to seed this year or not
- `shade_years` : bool, indicating whether to shade this year or not
- `n_site_int` : int, number of sites to intervene on
- `available_space` : vector/matrix : space available at each site (`k` value)
- `depth` : vector of site ids found to be within desired depth range
"""
function unguided_site_selection(prefseedsites, prefshadesites, seed_years, shade_years, n_site_int, available_space, depth)
    # Unguided deployment, seed/shade corals anywhere so long as available_space > 0.0
    # Only sites that have available space are considered, otherwise a zero-division error may occur later on.

    # Select sites (without replacement to avoid duplicate sites)
    candidate_sites = depth[(available_space.>0.0)[depth]]  # Filter down to site ids to be considered
    num_sites = length(candidate_sites)
    s_n_site_int = num_sites < n_site_int ? num_sites : n_site_int

    if seed_years
        prefseedsites = zeros(Int64, n_site_int)
        prefseedsites[1:s_n_site_int] .= StatsBase.sample(candidate_sites, s_n_site_int; replace=false)
    end

    if shade_years
        prefshadesites = zeros(Int64, n_site_int)
        prefshadesites[1:s_n_site_int] .= StatsBase.sample(candidate_sites, s_n_site_int; replace=false)
    end

    return prefseedsites[prefseedsites.>0], prefshadesites[prefshadesites.>0]
end
