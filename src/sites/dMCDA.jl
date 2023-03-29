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
    [MoosraMethod(), true],
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
    create_criteria_df(site_ids::AbstractArray, coral_cover::AbstractArray,
                coral_space::AbstractArray, connectivity_in::AbstractArray, 
                connectivity_out::AbstractArray, heat_stress::AbstractArray, 
                wave_stress::AbstractArray, criteria...)    

    Constructs the criteria dataframe for performing site selection.

# Arguments
- `site_ids` : site ids as integers.
- `coral_cover` : array containing coral cover (m^2) summed across species for each site.
- `coral_space` : array containing space available for coral (m^2) for each site.
- `connectivity_in` : array containing in-coming connectivity for each site.
- `connectivity_out` : array containing out-going connectivity for each site.
- `heat_stress` : array containing heat stress for each site.
- `wave_stress` : array containing wave stress for each site.
- `criteria...` : any number of additional criteria to use in site selection, input as Tuples
                    ("criteria name", [criteria array]).

"""
function create_criteria_store(site_ids::AbstractArray; criteria...)

    criteria_matrix = site_ids
    for crit_key in keys(criteria)
        criteria_matrix = hcat(criteria_matrix, criteria[crit_key][site_ids])
    end
    return NamedDimsArray(criteria_matrix[:, 2:end], reefs=site_ids, criteria=collect(keys(criteria)))

end

function create_criteria_store(site_ids::AbstractArray, domain::Domain, criteria::DataFrame)
    criteria_names = collect(criteria.columns)

    criteria_names = replace.(replace.(criteria_names, "_seed" => ""), "_shade" => "")
    criteria_names = Symbol.(unique(criteria_names[!occursin.("tol", criteria_names)]))
    criteria_values = Tuple(getproperty(domain, criteria_n) for criteria_n in criteria_names)
    criteria = (; zip(criteria_names, criteria_values)...)

    return create_criteria_store(site_ids; criteria)
end


function create_tolerances_store(; tolerances...)
    tol_store = [x -> tolerances[tol_key][2](x, tolerances[tol_key][1]) for tol_key in keys(tolerances)]
    return NamedDimsArray(tol_store, tols=collect(keys(tolerances)))
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
function rank_sites!(S, weights, rankings, n_site_int, mcda_func, rank_col)::Tuple{Vector{Int64},Matrix{Union{Float64,Int64}}}
    # Filter out all non-preferred sites
    selector = vec(.!all(S .== 0, dims=1))

    # weights in order of: in_conn, out_conn, wave, heat, predecessors, low cover
    weights = weights[selector]
    S = S[:, Bool[1, selector...]]

    s_order = retrieve_ranks(S, weights, mcda_func)

    last_idx = min(n_site_int, size(s_order, 1))
    prefsites = Int.(s_order[1:last_idx, 1])

    # Match by site_id and assign rankings to log
    align_rankings!(rankings, s_order, rank_col)

    return prefsites, s_order
end

function rank_seed_sites!(S, weights, rankings, n_site_int, mcda_func)::Tuple{Vector{Int64},Matrix{Union{Float64,Int64}}}
    rank_sites!(S, weights, rankings, n_site_int, mcda_func, 2)
end
function rank_shade_sites!(S, weights, rankings, n_site_int, mcda_func)::Tuple{Vector{Int64},Matrix{Union{Float64,Int64}}}
    rank_sites!(S, weights, rankings, n_site_int, mcda_func, 3)
end

"""
    create_decision_matrix(criteria_df::DataFrame, tolerances::DataFrame)

# Arguments
- `criteria_df` : contains criteria in each column for sites in each row.
- `tolerances` : contains thresholds for specified criteria, with names matching those in criteria_df.
                First row is the threshold value, second is 'gt' if criteria should be greater
                than threshold and 'lt' if criteria should be less than.

# Returns
- `A` : Decision matrix
- 'filtered': indices for sites not filtered due to threshold specifications.
"""
function create_decision_matrix(criteria_store::NamedDimsArray, tolerances::NamedDimsArray)

    for tol_key in tolerances.tols
        rule = map(tolerances(tol_key), criteria_store(tol_key))
        criteria_store = criteria_store[rule, :]
    end

    return criteria_store
end

"""
    create_intervention_matrix(A::Matrix, weights::DataFrame, criteria_df::DataFrame, int_crit_names::Vector{String})

# Arguments
- `A` : Criteria_df as a Matrix and filtered according to criteria thresholds set in tolerances.
- `weights` : contains weights for all criteria in criteria_df (not including site_ids).
- 

# Returns
- `A` : Decision matrix
- 'filtered': indices for sites not filtered due to threshold specifications.
- `criteria_df` : contains criteria for site selection in each column for sites in each row.
- 'int_crit_names': specifies criteria to be use in this decision instance (seeding, shading etc).
                    Must be a subset of the columns of criteria_df.

"""
function create_intervention_matrix(criteria_store::NamedDimsArray, params::NamedDimsArray, int_type::String)
    # Define intervention decision matrix
    weights_ind = occursin.("iv__", params.factors) .& occursin.(int_type, params.factors)
    temp_names = split.(params.factors[weights_ind], "__")
    crit_names = [findall(criteria_store.criteria .== Symbol(temp_name[2]))[1] for temp_name in temp_names]

    ws = mcda_normalize(Array(params[factors=findall(weights_ind)]))
    S = Matrix(criteria_store[criteria=crit_names])
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

    A = create_decision_matrix(criteria_store, thresholds)
    if isempty(A)
        # if all rows have nans and A is empty, abort mission
        return zeros(Int64, length(prefseedsites)), zeros(Int64, length(prefshadesites)), rankingsin
    end

    # cap to number of sites left after risk filtration
    n_site_int = min(n_site_int, length(A[:, 1]))

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
        prefseedsites, s_order_seed = rank_seed_sites!(SE, wse, rankings, n_site_int, mcda_func)
        if use_dist != 0
            prefseedsites, rankings = distance_sorting(prefseedsites, s_order_seed, distances, minimum_distance, rankings, 2)
        end
    end

    if log_shade && isempty(SH)
        prefshadesites = zeros(Int64, n_site_int)
    elseif log_shade
        prefshadesites, s_order_shade = rank_shade_sites!(SH, wsh, rankings, n_site_int, mcda_func)
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
    distance_sorting(pref_sites::AbstractArray{Int}, site_order::AbstractVector, dist::Array{Float64}, dist_thresh::Float64, top_n::Int64)::AbstractArray{Int}

Find selected sites with distances between each other < median distance-dist_thresh*(median distance).
Replaces these sites with sites in the top_n ranks if the distance between these sites is greater.

# Arguments
- `pref_sites` : original n highest ranked sites selected for seeding or shading.
- `site_order` : current order of ranked sites in terms of numerical site ID.
- `dist` : Matrix of unique distances between sites.
- `min_dist` : minimum distance between sites for selected sites.
- `top_n` : number of top ranked sites to re-select from.

# Returns
- `prefsites` : new set of selected sites for seeding or shading.
"""
function distance_sorting(pref_sites::AbstractArray{Int}, s_order::Matrix{Union{Float64,Int64}}, dist::Array{Float64},
    min_dist::Float64, rankings::Matrix{Int64}, rank_col::Int64)::Tuple{Vector{Union{Float64,Int64}},Matrix{Int64}}
    # set-up
    n_sites = length(pref_sites)
    site_order = s_order[:, 1]

    # sites to select alternatives from
    alt_sites = setdiff(site_order, pref_sites)[1:length(site_order)-n_sites]

    # find all selected sites closer than the min distance
    pref_dists = findall(dist[pref_sites, pref_sites] .< min_dist)
    # indices to replace
    inds_rep = sort(unique(reinterpret(Int64, pref_dists)))
    # number of sites to replace
    select_n = length(inds_rep)
    # indices to keep
    inds_keep = collect(1:length(pref_sites))
    inds_keep = setdiff(inds_keep, inds_rep)

    # storage for new set of sites
    rep_sites = pref_sites

    while (length(alt_sites) .>= select_n)
        rep_sites = [rep_sites[inds_keep[:]]; alt_sites[1:select_n]]

        # Find all sites within these highly ranked but unselected sites which are further apart
        alt_dists = dist[rep_sites, rep_sites] .> min_dist

        # Select from these sites those far enough away from all sites
        inds_keep = sum(alt_dists, dims=2) .== n_sites - 1

        # Keep sites that were far enough away last iteration
        inds_keep[1:end-select_n] .= true
        if length(inds_keep) == n_sites
            select_n = 0
            break
        else
            # remove checked alt_sites
            alt_sites = setdiff(alt_sites, alt_sites[1:select_n])
            select_n = sum(.!inds_keep)
        end
    end

    # If not all sites could be replaced, just use highest ranked remaining pref_sites
    if (select_n != 0) && !isempty(setdiff(pref_sites, rep_sites))
        rem_pref_sites = setdiff(pref_sites, rep_sites)
        rep_sites[end-select_n+1:end] .= rem_pref_sites[1:select_n]
    end

    new_site_order = setdiff(site_order, rep_sites)
    new_site_order = [rep_sites; new_site_order]
    s_order[:, 1] .= new_site_order
    # Match by site_id and assign rankings to log
    align_rankings!(rankings, s_order, rank_col)

    return rep_sites, rankings
end

function retrieve_ranks(S::Matrix, weights::Array, mcda_func::Function, site_ids::Array{Int64})
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
function retrieve_ranks(S::Matrix, scores::Array{Float64}, rev_val::Bool, site_ids::Array{Int64})
    s_order = Union{Float64,Int64}[Int.(site_ids) scores 1:size(S, 1)]
    s_order .= sortslices(s_order, dims=1, by=x -> x[2], rev=rev_val)
    @views s_order[:, 3] .= Int.(1:size(S, 1))

    return s_order
end


"""
    site_selection(domain::Domain, criteria::DataFrameRow{DataFrame,DataFrames.Index}, w_scens::NamedDimsArray, dhw_scens::NamedDimsArray, sum_cover::AbstractArray, area_to_seed::Float64)

Perform site selection using a chosen mcda aggregation method, domain, initial cover, criteria weightings and thresholds.

# Arguments
- `criteria` : contains criteria weightings and thresholds for a single scenario.
- `mcda_vars` : site selection criteria and weightings structure
- `w_scens` : array of length nsites containing wave scenario.
- `dhw_scens` : array of length nsites containing dhw scenario.
- `sum_cover` : summed cover (over species) for single scenario being run, for each site.
- `area_to_seed` : area of coral to be seeded at each time step in km^2

# Returns
- `ranks` : n_reps * sites * 3 (last dimension indicates: site_id, seeding rank, shading rank)
    containing ranks for single scenario.
"""
function site_selection(domain::Domain, criteria::DataFrameRow, tolerances::Tuple, site_ids::AbstractArray)::Matrix{Int64}

    criteria_store = create_criteria_store(domain, criteria, site_ids)
    tolerances_store = create_tolerance_store(tolerances)

    min_distance = domain.median_site_distance .* criteria.min_dist
    n_sites = length(site_ids)
    # site_id, seeding rank, shading rank
    rankingsin = [site_ids zeros(Int64, (n_sites, 1)) zeros(Int64, (n_sites, 1))]
    prefseedsites = zeros(Int64, (1, domain.sim_constants.n_site_int))
    prefshadesites = zeros(Int64, (1, domain.sim_constants.n_site_int))

    (_, _, ranks) = guided_site_selection(criteria_store, criteria, tolerances_store, domain.site_distances,
        min_distance, true, true, prefseedsites, prefshadesites, rankingsin)
    return ranks
end


"""
    run_site_selection(domain::Domain, criteria::DataFrame, sum_cover::AbstractArray, area_to_seed::Float64, time_step::Int64)

Perform site selection for a given domain for multiple scenarios defined in a dataframe.

# Arguments
- `domain` : ADRIA Domain type, indicating geographical domain to perform site selection over.
- `scenarios` : DataFrame of criteria weightings and thresholds for each scenario.
- `sum_cover` : array of size (number of scenarios * number of sites) containing the summed coral cover for each site selection scenario.
- `area_to_seed` : area of coral to be seeded at each time step in km^2
- `target_seed_sites` : list of candidate locations for seeding (indices)
- `target_shade_sites` : list of candidate location to shade (indices)

# Returns
- `ranks_store` : number of scenarios * sites * 3 (last dimension indicates: site_id, seed rank, shade rank)
    containing ranks for each scenario run.
"""
function run_site_selection(domain::Domain, criteria::DataFrame, tolerances::Tuple)
    ranks_store = NamedDimsArray(
        zeros(nrow(scenarios), length(dom.site_ids), 3),
        scenarios=1:nrow(scenarios),
        sites=dom.site_ids,
        ranks=["site_id", "seed_rank", "shade_rank"],
    )

    site_data = domain.site_data

    # Pre-calculate maximum depth to consider
    scenarios[:, "max_depth"] .= scenarios.depth_min .+ scenarios.depth_offset
    target_dhw_scens = unique(scenarios[:, "dhw_scenario"])
    target_wave_scens = unique(scenarios[:, "wave_scenario"])

    target_site_ids = Int64[]
    if !isnothing(target_seed_sites)
        append!(target_site_ids, target_seed_sites)
    end

    if !isnothing(target_shade_sites)
        append!(target_site_ids, target_shade_sites)
    end

        ranks_temp = site_selection(
            domain,
            scen_criteria,
            tolerances,
            depth_priority,
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
