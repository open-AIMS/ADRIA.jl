"""Objects and methods for Dynamic Multi-Criteria Decision Analysis/Making"""

using StatsBase
using Distances
using Combinatorics
using JMcDM
using InteractiveUtils: subtypes
using ADRIA: order_ranking, adria_vikor, adria_topsis


struct DMCDA_vars  # {V, I, F, M} where V <: Vector
    site_ids  # ::V
    n_site_int  # ::I
    priority_sites  # ::V
    priority_zones # ::V
    zones # ::V
    strong_pred  # ::V
    conn
    dam_prob  # ::A
    heat_stress_prob  # ::A
    site_depth #::V
    sum_cover  # ::F
    max_cover  # ::V
    area  # ::M
    min_area # ::F
    risk_tol  # ::F
    dist # ::M
    use_dist # ::Int64
    min_dist # ::Float64
    top_n # ::Int64
    wt_in_conn_seed  # ::F
    wt_out_conn_seed  # ::F
    wt_conn_shade  # ::F
    wt_waves # ::F
    wt_heat  # ::F
    wt_hi_cover  # ::F
    wt_lo_cover  # ::F
    wt_predec_seed  # ::F
    wt_predec_shade  # ::F
    wt_zones_seed # ::F
    wt_zones_shade # ::F
end

const jmcdm_methods = subtypes(MCDMMethod)

jmcdm_ignore = [
    JMcDM.CRITIC.CriticMethod,
    JMcDM.COPRAS.CoprasMethod,
    JMcDM.MOOSRA.MoosraMethod,
    JMcDM.MEREC.MERECMethod,
    JMcDM.ELECTRE.ElectreMethod,
    JMcDM.PROMETHEE.PrometheeMethod
]

const methods_mcda = [
    order_ranking,
    adria_vikor,
    adria_topsis,
    setdiff(jmcdm_methods, jmcdm_ignore)...
]

"""
    DMCDA_vars(domain::Domain, criteria::NamedDimsArray,
               site_ids::AbstractArray, sum_cover::AbstractArray, area_to_seed::Float64,
               waves::AbstractArray, dhws::AbstractArray)::DMCDA_vars
    DMCDA_vars(domain::Domain, criteria::NamedDimsArray, site_ids::AbstractArray,
               sum_cover::AbstractArray, area_to_seed::Float64)::DMCDA_vars
    DMCDA_vars(domain::Domain, criteria::DataFrameRow, site_ids::AbstractArray,
               sum_cover::AbstractArray, area_to_seed::Float64)::DMCDA_vars
    DMCDA_vars(domain::Domain, criteria::DataFrameRow, site_ids::AbstractArray,
               sum_cover::AbstractArray, area_to_seed::Float64,
               waves::AbstractArray, dhw::AbstractArray)::DMCDA_vars

Constuctors for DMCDA variables.
"""
function DMCDA_vars(domain::Domain, criteria::NamedDimsArray,
    site_ids::AbstractArray, sum_cover::AbstractArray, area_to_seed::Float64,
    waves::AbstractArray, dhws::AbstractArray)::DMCDA_vars

    # Site Data
    site_d = domain.site_data
    n_sites = n_locations(domain)
    area = site_area(domain)

    mcda_vars = DMCDA_vars(
        site_ids,
        domain.sim_constants.n_site_int,
        domain.sim_constants.priority_sites,
        domain.sim_constants.priority_zones,
        site_d.zone_type,
        domain.strong_pred,
        domain.TP_data .* site_k_area(domain),
        waves,
        dhws,
        site_d.depth_med,
        sum_cover,
        site_k(domain),
        area,
        criteria("coral_cover_tol") .* area_to_seed,
        criteria("deployed_coral_risk_tol"),
        domain.site_distances,
        criteria("use_dist"),
        domain.median_site_distance - domain.median_site_distance * criteria("dist_thresh"),
        criteria("top_n"),
        criteria("in_seed_connectivity"),
        criteria("out_seed_connectivity"),
        criteria("shade_connectivity"),
        criteria("wave_stress"),
        criteria("heat_stress"),
        criteria("coral_cover_high"),
        criteria("coral_cover_low"),
        criteria("seed_priority"),
        criteria("shade_priority"),
        criteria("zone_seed"),
        criteria("zone_shade")
    )

    return mcda_vars
end
function DMCDA_vars(domain::Domain, criteria::NamedDimsArray, site_ids::AbstractArray, sum_cover::AbstractArray, area_to_seed::Float64)::DMCDA_vars
    num_sites = n_locations(domain)
    return DMCDA_vars(domain, criteria, site_ids, sum_cover, area_to_seed, zeros(num_sites, 1), zeros(num_sites, 1))
end
function DMCDA_vars(domain::Domain, criteria::DataFrameRow, site_ids::AbstractArray,
    sum_cover::AbstractArray, area_to_seed::Float64, waves::AbstractArray, dhw::AbstractArray)::DMCDA_vars

    criteria_vec::NamedDimsArray = NamedDimsArray(collect(criteria), rows=names(criteria))
    return DMCDA_vars(domain, criteria_vec, site_ids, sum_cover, area_to_seed, waves, dhw)
end
function DMCDA_vars(domain::Domain, criteria::DataFrameRow, site_ids::AbstractArray,
    sum_cover::AbstractArray, area_to_seed::Float64)::DMCDA_vars

    criteria_vec::NamedDimsArray = NamedDimsArray(collect(criteria), rows=names(criteria))
    return DMCDA_vars(domain, criteria_vec, site_ids, sum_cover, area_to_seed)
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
    rank_sites!(S, weights, rankings, n_site_int, mcda_func, rank_col)

# Arguments
- `S` : Matrix, Site preference values
- `weights` : weights to apply
- `rankings` : vector of site ranks to update
- `n_site_int` : number of sites to select for interventions
- `mcda_func` : function or JMcDM DataType, designates mcda method to use
- `rank_col` : column to fill with rankings (2 for seed, 3 for shade)

# Returns
- `prefsites` : sites in order of their rankings
"""
function rank_sites!(S, weights, rankings, n_site_int, mcda_func, rank_col)::Tuple{Vector{Int64},Matrix{Union{Float64,Int64}}}
    # Filter out all non-preferred sites
    selector = vec(.!all(S[:, 2:end] .== 0, dims=1))

    # weights in order of: in_conn, out_conn, wave, heat, predecessors, low cover
    weights = weights[selector]
    S = S[:, Bool[1, selector...]]

    s_order = retrieve_ranks(S[:, 2:end], S[:, 1], weights, mcda_func)

    last_idx = min(n_site_int, size(s_order, 1))
    prefsites = Int64.(s_order[1:last_idx, 1])

    # Match by site_id and assign rankings to log
    align_rankings!(rankings, s_order, rank_col)

    return prefsites, s_order
end

"""
    retrieve_ranks(S::Matrix, site_ids::Vector, weights::Vector{Float64}, mcda_func::Function)
    retrieve_ranks(S::Matrix, site_ids::Vector, weights::Vector{Float64}, mcda_func::Type{<:MCDMMethod})
    retrieve_ranks(S::Matrix, site_ids::Vector, scores::Vector, maximize::Bool)

Get location ranks using mcda technique specified in mcda_func, weights and a decision matrix S.

# Arguments
- `S` : decision matrix containing criteria values for each location (n locations)*(m criteria)
- `site_ids` : array of site ids still remaining after filtering.
- `weights` : importance weights for each criteria. 
- `mcda_func` : function/JMcDM DataType to use for mcda, specified as an element from methods_mcda.
- `scores` : set of scores derived from applying an mcda ranking method.
- `maximize` : Boolean indicating whether a mcda method is maximizing score (true), or minimizing (false). 

# Returns
- `s_order` : [site_ids, criteria values, ranks]
"""
function retrieve_ranks(S::Matrix, site_ids::Vector, weights::Vector{Float64}, mcda_func::Function)
    S = mcda_normalize(S) .* weights'
    scores = mcda_func(S)

    return retrieve_ranks(S, site_ids, vec(scores), true)
end
function retrieve_ranks(S::Matrix, site_ids::Vector, weights::Vector{Float64}, mcda_func::Type{<:MCDMMethod})
    fns = fill(maximum, length(weights))
    results = mcdm(MCDMSetting(S, weights, fns), mcda_func())
    maximize = results.bestIndex == argmax(results.scores)

    return retrieve_ranks(S, site_ids, results.scores, maximize)
end
function retrieve_ranks(S::Matrix, site_ids::Vector, scores::Vector, maximize::Bool)
    s_order = Union{Float64,Int64}[Int64.(site_ids) scores Int64.(1:size(S, 1))]
    s_order .= sortslices(s_order, dims=1, by=x -> x[2], rev=maximize)

    return s_order
end


"""
    create_decision_matrix(site_ids, in_conn, out_conn, sum_cover, max_cover, area, wave_stress, heat_stress, predec, risk_tol)

Creates criteria matrix `A`, where each column is a selection criterium and each row is a site.
Sites are then filtered based on heat and wave stress risk.

Where no sites are filtered, size of ``A := n_sites × 7 criteria``.

Columns indicate:
1. Site ID
2. Incoming Node Connectivity Centrality
3. Outgoing Node Connectivity Centrality
4. Wave stress
5. Heat stress
6. Priority Predecessors
7. Available Area (relative to max cover)

# Arguments
- `site_ids` : vector of site ids
- `in_conn` : site incoming centrality (relative strength of connectivity) (0 <= c <= 1.0)
- `out_conn` : site outgoing centrality (relative strength of connectivity) (0 <= c <= 1.0)
- `sum_cover` : vector, sum of coral cover (across species) for each site (i.e., [x₁, x₂, ..., xₙ] where x_{1:n} <= 1.0)
- `max_cover` : maximum possible proportional coral cover (k) for each site, relative to total site area (k <= 1.0)
- `area` : total absolute area (in m²) for each site
- `wave_stress` : Probability of wave damage
- `heat_stress` : Probability of site being affected by heat stress
- `predec` : list of priority predecessors (sites strongly connected to priority sites)
- `risk_tol` : tolerance for wave and heat risk (∈ [0,1]). Sites with heat or wave risk> risk_tol are filtered out.
"""
function create_decision_matrix(site_ids, in_conn, out_conn, sum_cover, max_cover, area, wave_stress, heat_stress, site_depth, predec, zones_criteria, risk_tol)
    A = zeros(length(site_ids), 9)
    A[:, 1] .= site_ids  # Column of site ids

    # node connectivity centrality, need to instead work out strongest predecessors to priority sites
    A[:, 2] .= maximum(in_conn) != 0.0 ? in_conn / maximum(in_conn) : 0.0
    A[:, 3] .= maximum(out_conn) != 0.0 ? out_conn / maximum(out_conn) : 0.0

    # Wave damage, account for cases where no chance of damage or heat stress
    # if max > 0 then use damage probability from wave exposure
    A[:, 4] .= maximum(wave_stress) != 0.0 ? (wave_stress .- minimum(wave_stress)) ./ (maximum(wave_stress) - minimum(wave_stress)) : 0.0
    #A[:, 4] .= wave_stress ./ maximum(wave_stress)

    # risk from heat exposure
    # A[:, 5] .= maximum(heat_stress) != 0.0 ? (heat_stress .- minimum(heat_stress)) ./ (maximum(heat_stress) - minimum(heat_stress)) : 0.0
    A[:, 5] .= min.(heat_stress ./ 20.0, 1.0)

    # priority predecessors
    A[:, 6] .= predec[:, 3]

    # priority zone predecessors and sites
    A[:, 7] .= zones_criteria

    # Proportion of empty space (no coral) compared to max possible cover
    A[:, 8] = max.((max_cover - sum_cover), 0.0) .* area

    A[:, 9] = site_depth

    # Filter out sites that have high risk of wave damage, specifically
    # exceeding the risk tolerance
    A[A[:, 4].>risk_tol, 4] .= NaN
    rule = (A[:, 4] .<= risk_tol) .& (A[:, 5] .> risk_tol)
    A[rule, 5] .= NaN

    filtered = vec(.!any(isnan.(A), dims=2))
    # remove rows with NaNs
    A = A[filtered, :]
    return A, filtered
end


"""
    create_seed_matrix(A, min_area, in_conn_seed, out_conn_seed, waves, heat, predec, low_cover)

Create seeding specific decision matrix from criteria matrix. The weight criteria and filter.

# Arguments
- `A` : Criteria matrix
- `min_area` : Minimum available area for a site to be considered
- `wt_in_conn_seed` : Seed connectivity weight for seeding
- `wt_out_conn_seed` : Seed connectivity weight for seeding
- `wt_waves` : Wave stress weight
- `wt_heat` : Heat stress weight
- `wt_predec_seed` : Priority predecessor weight
- `wt_predec_zones_seed` : Priority zones weight for seeding
- `wt_low_cover` : Weighting for low coral cover (coral real estate), when seeding

# Returns
Tuple (SE, wse)
- `SE` : Matrix of shape [n sites considered, 7]
    1. Site index ID
    2. Incoming Centrality
    3. Outgoing Centrality
    4. Wave risk (higher values = less risk)
    5. Damage risk (higher values = less risk)
    6. Priority predecessors relating to coral real estate relative to max capacity
    7. Available space
- `wse` : 5-element vector of criteria weights
    1. incoming connectivity
    2. outgoing connectivity
    3. wave
    4. heat
    5. seed predecessors (weights importance of sites highly connected to priority sites for seeding)
    6. seed zones (weights importance of sites highly connected to or within priority zones for seeding)
    7. low cover (weights importance of sites with low cover/high available real estate to plant corals)
"""
function create_seed_matrix(A, min_area, wt_in_conn_seed, wt_out_conn_seed, wt_waves, wt_heat, wt_predec_seed, wt_predec_zones_seed, wt_lo_cover)
    # Define seeding decision matrix, based on copy of A
    SE = copy(A)

    wse = [wt_in_conn_seed, wt_out_conn_seed, wt_waves, wt_heat, wt_predec_seed, wt_predec_zones_seed, wt_lo_cover, wt_heat]
    # wse .= mcda_normalize(wse)

    SE[:, 4] = (1 .- SE[:, 4]) # compliment of wave risk
    SE[:, 5] = (1 .- SE[:, 5]) # compliment of heat risk

    # Coral real estate as total area, sites with ≤ min_area to be seeded available filtered out
    SE[vec(A[:, 8] .<= min_area), 8] .= NaN

    # Mark "hot" locations for filter
    SE[vec(A[:, 5] .>= 0.75), 5] .= NaN

    # Filter out identified locations
    SE = SE[vec(.!any(isnan.(SE), dims=2)), :]

    return SE, wse
end


"""
    create_shade_matrix(A, wt_conn_shade , wt_waves, wt_heat, wt_predec_shade, wt_hi_cover)

Create shading specific decision matrix and apply weightings.

# Arguments
- `A` : Criteria  matrix
- `wt_conn_shade` : Shading connectivity weight
- `wt_waves` : Wave stress weight
- `wt_heat` : Heat stress weight
- `wt_predec_zones_shade` : Priority zones weight for shading
- `wt_predec_shade` : Priority predecessor weight for shading
- `wt_hi_cover` : Weighting for high coral cover when shading

# Returns
Tuple (SH, wsh)
- `SH` : Matrix of shape [n sites considered, 7]
    1. Site index ID
    2. Incoming Centrality
    3. Outgoing Centrality
    4. Wave risk (higher values = less risk)
    5. Damage risk (higher values = less risk)
    6. Priority predecessors relating to coral real estate relative to max capacity
    7. Available space
- `wsh` : 5-element vector of criteria weights
    1. shade connectivity
    2. wave
    3. heat
    4. shade predecessors (weights importance of sites highly connected to priority sites for shading)
    4. shade zones (weights importance of sites highly connected to or within priority zones)
    5. high cover (weights importance of sites with high cover of coral to shade)
"""
function create_shade_matrix(A, max_area, wt_conn_shade, wt_waves, wt_heat, wt_predec_shade, wt_predec_zones_shade, wt_hi_cover)
    # Set up decision matrix to be same size as A
    SH = copy(A)

    # Remove consideration of site depth as shading not accounted for in bleaching model yet
    SH = SH[:, 1:end-1]

    wsh = [wt_conn_shade, wt_conn_shade, wt_waves, wt_heat, wt_predec_shade, wt_predec_zones_shade, wt_hi_cover]
    wsh .= mcda_normalize(wsh)

    SH[:, 4] = (1.0 .- A[:, 4]) # complimentary of wave damage risk
    SH[:, 8] = (max_area .- A[:, 8]) # total area of coral cover

    SH[SH[:, 8].<0, 8] .= 0  # if any negative, scale back to zero
    return SH, wsh
end


"""
    guided_site_selection(d_vars::DMCDA_vars, alg_ind::Int64, log_seed::Bool, log_shade::Bool, prefseedsites::AbstractArray{Int64}, prefshadesites::AbstractArray{Int64}, rankingsin::Matrix{Int64})

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
function guided_site_selection(
    d_vars::DMCDA_vars, alg_ind::T,
    log_seed::B, log_shade::B,
    prefseedsites::IA, prefshadesites::IB,
    rankingsin::Matrix{T},
    in_conn::Vector{Float64},
    out_conn::Vector{Float64},
    strong_pred::Vector{Int64};
    methods_mcda=methods_mcda
)::Tuple where {T<:Int64,IA<:AbstractArray{<:Int64},IB<:AbstractArray{<:Int64},B<:Bool}

    use_dist::Int64 = d_vars.use_dist
    min_dist::Float64 = d_vars.min_dist
    site_ids = copy(d_vars.site_ids)
    n_sites::Int64 = length(site_ids)

    # if no sites are available, abort
    if n_sites == 0
        return zeros(Int64, length(prefseedsites)), zeros(Int64, length(prefshadesites)), rankingsin
    end

    n_site_int::Int64 = d_vars.n_site_int
    priority_sites::Array{Int64} = d_vars.priority_sites[in.(d_vars.priority_sites, [site_ids])]
    priority_zones::Array{String} = d_vars.priority_zones

    zones = d_vars.zones[site_ids]
    wave_stress = d_vars.dam_prob[site_ids]
    heat_stress = d_vars.heat_stress_prob[site_ids]
    site_depth = d_vars.site_depth[site_ids]
    sum_cover = d_vars.sum_cover[site_ids]
    max_cover = d_vars.max_cover[site_ids]
    area = d_vars.area[site_ids]

    risk_tol = d_vars.risk_tol
    w_in_conn = d_vars.wt_in_conn_seed
    w_out_conn = d_vars.wt_out_conn_seed
    w_shade_conn = d_vars.wt_conn_shade
    w_waves = d_vars.wt_waves
    w_heat = d_vars.wt_heat
    w_high_cover = d_vars.wt_hi_cover
    w_low_cover = d_vars.wt_lo_cover
    w_predec_seed = d_vars.wt_predec_seed
    w_predec_shade = d_vars.wt_predec_shade
    w_predec_zones_seed = d_vars.wt_zones_seed
    w_predec_zones_shade = d_vars.wt_zones_shade

    # site_id, seeding rank, shading rank
    rankings = Int64[site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

    # work out which priority predecessors are connected to priority sites
    predec::Matrix{Float64} = zeros(n_sites, 3)
    predec[:, 1:2] .= strong_pred
    predprior = predec[in.(predec[:, 1], [priority_sites']), 2]
    predprior = Int64[x for x in predprior if !isnan(x)]

    predec[predprior, 3] .= 1.0

    # for zones, find sites which are zones and strongest predecessors of sites in zones
    zone_ids = intersect(priority_zones, unique(zones))
    zone_weights = mcda_normalize(collect(length(zone_ids):-1:1))
    zone_preds = zeros(n_sites)
    zone_sites = zeros(n_sites)

    for (k::Int64, z_name::String) in enumerate(zone_ids)
        # find sites which are strongest predecessors of sites in the zone
        zone_preds_temp::Vector{Int64} = strong_pred[zones.==z_name]
        for s::Int64 in unique(zone_preds_temp)
            # for each predecessor site, add zone_weights * (no. of zone sites the site is a strongest predecessor for)
            zone_preds[site_ids.==s] .= zone_preds[site_ids.==s] .+ (zone_weights[k] .* sum(zone_preds_temp .== s))
        end
        # add zone_weights for sites in the zone (whether a strongest predecessor of a zone or not)
        zone_sites[zones.==z_name] .= zone_weights[k]
    end

    # add weights for strongest predecessors and zones to get zone criteria
    zones_criteria = zone_preds .+ zone_sites
    mcda_func = methods_mcda[alg_ind]

    A, filtered_sites = create_decision_matrix(site_ids, in_conn, out_conn, sum_cover, max_cover, area, wave_stress, heat_stress, site_depth, predec, zones_criteria, risk_tol)
    if isempty(A)
        # if all rows have nans and A is empty, abort mission
        return zeros(Int64, length(prefseedsites)), zeros(Int64, length(prefshadesites)), rankingsin
    end

    # cap to number of sites left after risk filtration
    n_site_int = min(n_site_int, length(A[:, 1]))

    # if seeding, create seeding specific decision matrix
    if log_seed
        SE, wse = create_seed_matrix(A, d_vars.min_area, w_in_conn, w_out_conn, w_waves, w_heat, w_predec_seed, w_predec_zones_seed, w_low_cover)
    end

    # if shading, create shading specific decision matrix
    if log_shade
        max_area = (area.*max_cover)[filtered_sites]
        SH, wsh = create_shade_matrix(A, max_area, w_shade_conn, w_waves, w_heat, w_predec_shade, w_predec_zones_shade, w_high_cover)
    end

    if log_seed && isempty(SE)
        prefseedsites = zeros(Int64, n_site_int)
    elseif log_seed
        prefseedsites, s_order_seed = rank_sites!(SE, wse, rankings, n_site_int, mcda_func, 2)
        if use_dist != 0
            prefseedsites, rankings = distance_sorting(prefseedsites, s_order_seed, d_vars.dist, min_dist, Int64(d_vars.top_n), rankings, 2)
        end
    end

    if log_shade && isempty(SH)
        prefshadesites = zeros(Int64, n_site_int)
    elseif log_shade
        prefshadesites, s_order_shade = rank_sites!(SH, wsh, rankings, n_site_int, mcda_func, 3)
        if use_dist != 0
            prefshadesites, rankings = distance_sorting(prefshadesites, s_order_shade, d_vars.dist, min_dist, Int64(d_vars.top_n), rankings, 3)
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
    min_dist::Float64, top_n::Int64, rankings::Matrix{Int64}, rank_col::Int64)::Tuple{Vector{Union{Float64,Int64}},Matrix{Int64}}
    # set-up
    n_sites = length(pref_sites)
    site_order = s_order[:, 1]

    # sites to select alternatives from
    alt_sites = setdiff(site_order, pref_sites)[1:min(top_n, length(site_order) - n_sites)]

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


"""
    run_site_selection(domain::Domain, scenarios::DataFrame, sum_cover::AbstractArray, area_to_seed::Float64, time_step::Int64)

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
function run_site_selection(dom::Domain, scenarios::DataFrame, sum_cover::AbstractArray, area_to_seed::Float64;
    target_seed_sites=nothing, target_shade_sites=nothing)
    ranks_store = NamedDimsArray(
        zeros(nrow(scenarios), length(dom.site_ids), 3),
        scenarios=1:nrow(scenarios),
        sites=dom.site_ids,
        ranks=["site_id", "seed_rank", "shade_rank"],
    )

    dhw_scens = dom.dhw_scens
    wave_scens = dom.wave_scens

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

    n_sites = length(dom.site_ids)
    for (scen_idx, scen) in enumerate(eachrow(scenarios))
        depth_criteria = (dom.site_data.depth_med .<= scen.max_depth) .& (dom.site_data.depth_med .>= scen.depth_min)
        depth_priority = findall(depth_criteria)

        considered_sites = target_site_ids[findall(in(depth_priority), target_site_ids)]

        ranks_store(scenarios=scen_idx, sites=dom.site_ids[considered_sites]) .= site_selection(
            dom,
            scen,
            (mean(wave_scens[:, :, target_wave_scens], dims=(:timesteps, :scenarios)) .+ std(wave_scens[:, :, target_wave_scens], dims=(:timesteps, :scenarios))) .* 0.5,
            (mean(dhw_scens[:, :, target_dhw_scens], dims=(:timesteps, :scenarios)) .+ std(dhw_scens[:, :, target_dhw_scens], dims=(:timesteps, :scenarios))) .* 0.5,
            considered_sites,
            sum_cover[scen_idx, :],
            area_to_seed
        )
    end

    return ranks_store
end

"""
    site_selection(domain::Domain, scenario::DataFrameRow{DataFrame,DataFrames.Index}, w_scens::NamedDimsArray, dhw_scens::NamedDimsArray, site_ids::Vector{Int64}, sum_cover::AbstractArray, area_to_seed::Float64)

Perform site selection using a chosen mcda aggregation method, domain, initial cover, criteria weightings and thresholds.

# Arguments
- `scenario` : contains criteria weightings and thresholds for a single scenario.
- `mcda_vars` : site selection criteria and weightings structure
- `w_scens` : array of length nsites containing wave scenario.
- `dhw_scens` : array of length nsites containing dhw scenario.
- `site_ids` : locations to consider
- `sum_cover` : summed cover (over species) for single scenario being run, for each site.
- `area_to_seed` : area of coral to be seeded at each time step in km^2

# Returns
- `ranks` : n_reps * sites * 3 (last dimension indicates: site_id, seeding rank, shading rank)
    containing ranks for single scenario.
"""
function site_selection(domain::Domain, scenario::DataFrameRow, w_scens::NamedDimsArray, dhw_scens::NamedDimsArray,
    site_ids::Vector{Int64}, sum_cover::Vector{Float64}, area_to_seed::Float64)::Matrix{Int64}

    mcda_vars = DMCDA_vars(domain, scenario, site_ids, sum_cover, area_to_seed, w_scens, dhw_scens)
    n_sites = length(site_ids)

    # site_id, seeding rank, shading rank
    rankingsin = [mcda_vars.site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

    prefseedsites::Vector{Int64} = zeros(Int64, mcda_vars.n_site_int)
    prefshadesites::Vector{Int64} = zeros(Int64, mcda_vars.n_site_int)

    # Determine connectivity strength
    # Account for cases where no coral cover
    in_conn, out_conn, strong_pred = connectivity_strength(domain.TP_data .* site_k_area(domain), sum_cover)
    in_conn = in_conn[site_ids]
    out_conn = out_conn[site_ids]
    strong_pred = strong_pred[site_ids]

    (_, _, ranks) = guided_site_selection(mcda_vars, scenario.guided, true, true, prefseedsites, prefshadesites, rankingsin, in_conn, out_conn, strong_pred)

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
