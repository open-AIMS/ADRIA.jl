"""Objects and methods for Dynamic Multi-Criteria Decision Analysis/Making"""

using StatsBase
using Distances
using Combinatorics

struct DMCDA_vars  # {V, I, F, M} where V <: Vector
    site_ids  # ::V
    nsiteint  # ::I
    prioritysites  # ::V
    priorityzones # ::V
    zones # ::V
    strongpred  # ::V
    in_conn  # ::v
    out_conn  # ::v
    damprob  # ::A
    heatstressprob  # ::A
    sitedepth #::V
    sumcover  # ::F
    maxcover  # ::V
    area  # ::M
    min_area # ::F
    risktol  # ::F
    dist # ::M
    dist_thresh # ::Float64
    top_n # ::Int64
    wtinconnseed  # ::F
    wtoutconnseed  # ::F
    wtconshade  # ::F
    wtwaves  # ::F
    wtheat  # ::F
    wthicover  # ::F
    wtlocover  # ::F
    wtpredecseed  # ::F
    wtpredecshade  # ::F
    wtzonesseed # ::F
    wtzonesshade # ::F
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
        rankings[findall(rankings[:, 1] .== site_id), col] .= s_order[i, 3]
    end

    return
end

"""
    rank_sites!(S, weights, rankings, nsiteint, rank_col)
    rank_seed_sites!(S, weights, rankings, nsiteint)
    rank_shade_sites!(S, weights, rankings, nsiteint)

# Arguments
- S : Matrix, Site preference values
- weights : weights to apply
- rankings : vector of site ranks to update
- nsiteint : number of sites to select for interventions
- rank_col : column to fill with rankings (2 for seed, 3 for shade)

# Returns
prefsites : sites in order of their rankings
"""
function rank_sites!(S, weights, rankings, nsiteint, mcda_func, rank_col)::Tuple{Vector{Int64},Matrix{Union{Float64,Int64}}}
    # Filter out all non-preferred sites
    selector = vec(.!all(S[:, 2:end] .== 0, dims=1))

    # weights in order of: in_conn, out_conn, wave, heat, predecessors, low cover
    weights = weights[selector]
    S = S[:, Bool[1, selector...]]

    # Skip first column as this holds site index ids
    S[:, 2:end] = mcda_normalize(S[:, 2:end])
    S[:, 2:end] .= S[:, 2:end] .* repeat(weights', size(S[:, 2:end], 1), 1)
    s_order = mcda_func(S)

    last_idx = min(nsiteint, size(s_order, 1))
    prefsites = Int.(s_order[1:last_idx, 1])

    # Match by site_id and assign rankings to log
    align_rankings!(rankings, s_order, rank_col)

    return prefsites, s_order
end

function rank_seed_sites!(S, weights, rankings, nsiteint, mcda_func)::Tuple{Vector{Int64},Matrix{Union{Float64,Int64}}}
    rank_sites!(S, weights, rankings, nsiteint, mcda_func, 2)
end
function rank_shade_sites!(S, weights, rankings, nsiteint, mcda_func)::Tuple{Vector{Int64},Matrix{Union{Float64,Int64}}}
    rank_sites!(S, weights, rankings, nsiteint, mcda_func, 3)
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
- site_ids : vector of site ids
- in_conn : site incoming centrality (relative strength of connectivity) (0 <= c <= 1.0)
- out_conn : site outgoing centrality (relative strength of connectivity) (0 <= c <= 1.0)
- sum_cover : vector, sum of coral cover (across species) for each site (i.e., [x₁, x₂, ..., xₙ] where x_{1:n} <= 1.0)
- max_cover : maximum possible proportional coral cover (k) for each site, relative to total site area (k <= 1.0)
- area : total absolute area (in m²) for each site
- wave_stress : Probability of wave damage
- heat_stress : Probability of site being affected by heat stress
- predec : list of priority predecessors (sites strongly connected to priority sites)
- risk_tol : tolerance for wave and heat risk (∈ [0,1]). Sites with heat or wave risk> risktol are filtered out.
"""
function create_decision_matrix(site_ids, in_conn, out_conn, sum_cover, max_cover, area, wave_stress, heat_stress, site_depth, predec, zones_criteria, risk_tol)
    A = zeros(length(site_ids), 9)

    A[:, 1] .= site_ids  # Column of site ids

    # Account for cases where no coral cover
    c_cov_area = in_conn .* sum_cover .* area
    o_cov_area = out_conn .* sum_cover .* area

    # node connectivity centrality, need to instead work out strongest predecessors to priority sites
    A[:, 2] .= maximum(c_cov_area) != 0.0 ? c_cov_area / maximum(c_cov_area) : 0.0
    A[:, 3] .= maximum(o_cov_area) != 0.0 ? o_cov_area / maximum(o_cov_area) : 0.0

    # Wave damage, account for cases where no chance of damage or heat stress
    # if max > 0 then use damage probability from wave exposure
    A[:, 4] .= maximum(wave_stress) != 0.0 ? (wave_stress .- minimum(wave_stress)) ./ (maximum(wave_stress) - minimum(wave_stress)) : 0.0

    # risk from heat exposure
    A[:, 5] .= maximum(heat_stress) != 0.0 ? (heat_stress .- minimum(heat_stress)) ./ (maximum(heat_stress) - minimum(heat_stress)) : 0.0

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
    create_seed_matrix(A, min_area, inconn_seed, outconn_seed, waves, heat, predec, low_cover)

Create seeding specific decision matrix from criteria matrix. The weight criteria and filter.

# Arguments
- A : Criteria matrix
- min_area : Minimum available area for a site to be considered
- inconn_seed : Seed connectivity weight for seeding
- outconn_seed : Seed connectivity weight for seeding
- waves : Wave stress weight
- heat : Heat stress weight
- predec : Priority predecessor weight
- low_cover : Weighting for low coral cover (coral real estate), when seeding

# Returns
Tuple (SE, wse)
- SE : Matrix of shape [n sites considered, 7]
    1. Site index ID
    2. Incoming Centrality
    3. Outgoing Centrality
    4. Wave risk (higher values = less risk)
    5. Damage risk (higher values = less risk)
    6. Priority predecessors relating to coral real estate relative to max capacity
    7. Available space
- wse : 5-element vector of criteria weights
    1. incoming connectivity
    2. outgoing connectivity
    3. wave
    4. heat
    5. seed predecessors (weights importance of sites highly connected to priority sites for seeding)
    6. low cover (weights importance of sites with low cover/high available real estate to plant corals)
"""
function create_seed_matrix(A, min_area, inconn_seed, outconn_seed, waves, heat, predec, predec_zones_seed, low_cover)
    # Define seeding decision matrix, based on copy of A
    SE = copy(A)
    #remove consideration site depth
    SE = SE[:, 1:end-1]

    wse = [inconn_seed, outconn_seed, waves, heat, predec, predec_zones_seed, low_cover, heat]
    wse .= mcda_normalize(wse)

    SE[:, 4] = (1 .- SE[:, 4]) # compliment of wave risk
    SE[:, 5] = (1 .- SE[:, 5]) # compliment of heat risk

    # coral real estate as total area, sites with =<20% of area to be seeded available filtered out
    SE[vec(A[:, 8] .<= min_area), 8] .= NaN
    SE = SE[vec(.!any(isnan.(SE), dims=2)), :]

    return SE, wse
end


"""
    create_shade_matrix(A, wtconshade, wtwaves, wtheat, wtpredecshade, wthicover)

Create shading specific decision matrix and apply weightings.

# Arguments
- A : Criteria  matrix
- wtconshade : Shading connectivity weight
- wtwaves : Wave stress weight
- wtheat : Heat stress weight
- wtpredecshade : Priority predecessor weight for shading
- wthicover : Weighting for high coral cover when shading

# Returns
Tuple (SH, wsh)
- SH : Matrix of shape [n sites considered, 7]
    1. Site index ID
    2. Incoming Centrality
    3. Outgoing Centrality
    4. Wave risk (higher values = less risk)
    5. Damage risk (higher values = less risk)
    6. Priority predecessors relating to coral real estate relative to max capacity
    7. Available space
- wsh : 5-element vector of criteria weights
    1. shade connectivity
    2. wave
    3. heat
    4. shade predecessors (weights importance of sites highly connected to priority sites for shading)
    5. high cover (weights importance of sites with high cover of coral to shade)
"""
function create_shade_matrix(A, max_area, conn_shade, waves, heat, predec, predec_zones_shade, high_cover)
    # Set up decision matrix to be same size as A
    SH = copy(A)

    wsh = [conn_shade, conn_shade, waves, heat, predec, predec_zones_shade, high_cover, heat]
    wsh .= mcda_normalize(wsh)

    SH[:, 4] = (1.0 .- A[:, 4]) # complimentary of wave damage risk
    SH[:, 8] = (max_area .- A[:, 8]) # total area of coral cover

    SH[SH[:, 8].<0, 8] .= 0  # if any negative, scale back to zero
    return SH, wsh
end


"""
    dMCDA(d_vars::DMCDA_vars, alg_ind::Int64, log_seed::Bool, log_shade::Bool, prefseedsites::AbstractArray{Int}, prefshadesites::AbstractArray{Int}, rankingsin::Matrix{Int64})

# Arguments
- d_vars : DMCDA_vars type struct containing weightings and criteria values for site selection.
- alg_ind : integer indicating MCDA aggregation method to use (0: none, 1: order ranking, 2:topsis, 3: vikor)
- log_seed : boolean ideicating whether seeding sites are being re-assesed at current time
- log_shade : boolean ideicating whether shading/fogging sites are being re-assesed at current time
- prefshadesites : previous time step's selection of sites for shading
- prefseedsites : previous time step's selection of sites for seeding
- rankingsin : storage for site rankings

# Returns
Tuple :
    - prefseedsites : nsiteint highest ranked seeding sites
    - prefshadesites : nsiteint highest ranked shading/fogging sites
    - number of seed sites : nprefseedsites
    - nprefshadesites : number of shade sites
    - rankings : nsitesx3 matrix holding [site_id, seeding_rank, shading_rank],
        0 indicates sites that were not considered
"""
function dMCDA(d_vars::DMCDA_vars, alg_ind::Int64, log_seed::Bool, log_shade::Bool,
    prefseedsites::AbstractArray{Int}, prefshadesites::AbstractArray{Int},
    rankingsin::Matrix{Int64})::Tuple

    site_ids::Array{Int64} = copy(d_vars.site_ids)

    dist_thresh::Float64 = d_vars.dist_thresh

    # Force different sites to be selected
    site_ids = setdiff(site_ids, vcat(prefseedsites, prefshadesites))
    mod_n_ranks = min(size(rankingsin, 1), length(site_ids))
    if mod_n_ranks < length(d_vars.site_ids) && length(rankingsin) != 0
        rankingsin = rankingsin[in.(rankingsin[:, 1], [site_ids]), :]
        site_ids = rankingsin[:, 1]
    elseif length(rankingsin) != 0
        rankingsin = [site_ids zeros(Int64, length(site_ids)) zeros(Int64, length(site_ids))]
    end

    nsites::Int64 = length(site_ids)

    # if no sites are available, abort
    if nsites == 0
        return prefseedsites, prefshadesites, rankingsin
    end

    nsiteint::Int64 = d_vars.nsiteint
    prioritysites::Array{Int64} = d_vars.prioritysites[in.(d_vars.prioritysites, [site_ids])]
    priorityzones::Array{String} = d_vars.priorityzones

    strongpred = d_vars.strongpred[site_ids, :]
    in_conn = d_vars.in_conn[site_ids]
    out_conn = d_vars.out_conn[site_ids]
    zones = d_vars.zones[site_ids]
    wave_stress = d_vars.damprob[site_ids]
    heat_stress = d_vars.heatstressprob[site_ids]
    site_depth = d_vars.sitedepth[site_ids]
    sum_cover = d_vars.sumcover[site_ids]
    max_cover = d_vars.maxcover[site_ids]
    area = d_vars.area[site_ids]

    risk_tol = d_vars.risktol
    w_inconn = d_vars.wtinconnseed
    w_outconn = d_vars.wtoutconnseed
    w_shade_conn = d_vars.wtconshade
    w_waves = d_vars.wtwaves
    w_heat = d_vars.wtheat
    w_high_cover = d_vars.wthicover
    w_low_cover = d_vars.wtlocover
    w_predec_seed = d_vars.wtpredecseed
    w_predec_shade = d_vars.wtpredecshade
    w_predec_zones_seed = d_vars.wtzonesseed
    w_predec_zones_shade = d_vars.wtzonesshade

    # site_id, seeding rank, shading rank
    rankings = Int64[site_ids zeros(Int64, nsites) zeros(Int64, nsites)]

    # work out which priority predecessors are connected to priority sites
    predec::Array{Float64} = zeros(nsites, 3)
    predec[:, 1:2] .= strongpred
    predprior = predec[in.(predec[:, 1], [prioritysites']), 2]
    predprior = [x for x in predprior if !isnan(x)]

    predec[predprior, 3] .= 1.0

    # for zones, find sites which are zones and strongest predecessors of sites in zones
    zone_ids = intersect(priorityzones, unique(zones))
    zone_weights = mcda_normalize(collect(length(zone_ids):-1:1))
    zone_preds = zeros(nsites, 1)
    zone_sites = zeros(nsites, 1)

    for k in axes(zone_ids, 1)
        # find sites which are strongest predecessors of sites in the zone
        zone_preds_temp = strongpred[zones.==zone_ids[k]]
        for s in unique(zone_preds_temp)
            # for each predecessor site, add zone_weights* (no. of zone sites the site is a strongest predecessor for)
            zone_preds[site_ids.==s] .= zone_preds[site_ids.==s] .+ (zone_weights[k]) .* sum(zone_preds_temp .== s)
        end
        # add zone_weights for sites in the zone (whether a strongest predecessor of a zone or not)        
        zone_sites[zones.==zone_ids[k]] .= zone_weights[k]
    end

    # add weights for strongest predecessors and zones to get zone criteria
    zones_criteria = zone_preds .+ zone_sites

    A, filtered_sites = create_decision_matrix(site_ids, in_conn, out_conn, sum_cover, max_cover, area, wave_stress, heat_stress, site_depth, predec, zones_criteria, risk_tol)
    if isempty(A)
        # if all rows have nans and A is empty, abort mission
        return prefseedsites, prefshadesites, rankingsin
    end

    # cap to number of sites left after risk filtration
    nsiteint = min(nsiteint, length(A[:, 1]))

    # if seeding, create seeding specific decision matrix
    if log_seed
        SE, wse = create_seed_matrix(A, d_vars.min_area, w_inconn, w_outconn, w_waves, w_heat, w_predec_seed, w_predec_zones_seed, w_low_cover)
    end

    # if shading, create shading specific decision matrix
    if log_shade
        max_area = (area.*max_cover)[filtered_sites]
        SH, wsh = create_shade_matrix(A, max_area, w_shade_conn, w_waves, w_heat, w_predec_shade, w_predec_zones_shade, w_high_cover)
    end

    if alg_ind == 1
        mcda_func = order_ranking
    elseif alg_ind == 2
        mcda_func = topsis
    elseif alg_ind == 3
        mcda_func = vikor
    else
        error("Unknown MCDA algorithm selected. Valid options are 1 (Order Ranking), 2 (TOPSIS) and 3 (VIKOR).")
    end

    if log_seed && isempty(SE)
        prefseedsites = repeat([0], nsiteint)
    elseif log_seed
        prefseedsites, s_order_seed = rank_seed_sites!(SE, wse, rankings, nsiteint, mcda_func)
        if dist_thresh != 1.0
            prefseedsites .= distance_sorting(prefseedsites, s_order_seed[:, 1], d_vars.dist, dist_thresh, Int64(d_vars.top_n))
        end
    end

    if log_shade && isempty(SH)
        prefshadesites = repeat([0], nsiteint)
    elseif log_shade
        prefshadesites, s_order_shade = rank_shade_sites!(SH, wsh, rankings, nsiteint, mcda_func)
        if dist_thresh != 1.0
            prefshadesites .= distance_sorting(prefshadesites, s_order_shade[:, 1], d_vars.dist, dist_thresh, Int64(d_vars.top_n))
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
- pref_sites : original n highest ranked sites selected for seeding or shading.
- site_order : current order of ranked sites in terms of numerical site ID.
- dist : Matrix of unique distances between sites.
- dist_thresh : threshold for minimum deviance below the median distance between sites for selected sites.
- top_n : number of top ranked sites to re-select from.

# Returns
- prefsites : new set of selected sites for seeding or shading.
"""
function distance_sorting(pref_sites::AbstractArray{Int}, site_order::AbstractVector, dist::Array{Float64}, dist_thresh::Float64, top_n::Int64)::AbstractArray{Int}
    # set-up
    nsites = length(pref_sites)
    # sites to select alternatives from
    alt_sites = setdiff(site_order, pref_sites)[1:min(top_n, length(site_order) - nsites)]
    # minimum distance for filtering selected sites (within dist_thesh of median)
    min_dist = median(dist[.!isnan.(dist)]) - dist_thresh * median(dist[.!isnan.(dist)])

    # find all selected sites closer than the min distance
    pref_dists = findall(dist[pref_sites, pref_sites] .< min_dist)
    # indices to replace
    inds_rep = sort(unique(reinterpret(Int64, pref_dists)))
    # number of sites to replace
    select_n = length(inds_rep)
    # indices to keep
    inds_keep = [k for k in 1:length(pref_sites)]
    inds_keep = setdiff(inds_keep, inds_rep)

    # storage for new set of sites
    rep_sites = pref_sites

    while (length(alt_sites) .>= select_n)
        rep_sites = [rep_sites[inds_keep[:]]; alt_sites[1:select_n]]

        # Find all sites within these highly ranked but unselected sites which are further apart
        alt_dists = dist[rep_sites, rep_sites] .> min_dist

        # Select from these sites those far enough away from all sites
        inds_keep = sum(alt_dists, dims=2) .== nsites - 1

        # Keep sites that were far enough away last iteration
        inds_keep[1:end-select_n] .= true
        if length(inds_keep) == nsites
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

    return rep_sites
end

"""
    order_ranking(S::Array{Float64, 2})

Uses simple summation as aggregation method for decision criteria.
Then orders sites from highest aggregate score to lowest.

# Arguments
- S : Decision matrix (seeding or shading)

# Returns
- s_order : nsites × 3 matrix with columns
    1. site ids
    2. calculated site rank score (higher values = higher ranked)
    3. site order id
"""
function order_ranking(S::Array{Float64,2})::Array{Union{Float64,Int64},2}
    n::Int64 = size(S, 1)
    s_order::Array = Union{Float64,Int64}[zeros(Int, n) zeros(Float64, n) zeros(Int, n)]

    # Simple ranking - add criteria weighted values for each sites
    # Third column is derived from the number of sites for situations where
    # a subset of sites are being investigated (and so using their site IDs
    # will be inappropriate)
    @views s_order[:, 1] .= Int.(S[:, 1])
    @views s_order[:, 2] .= sum(S[:, 2:end], dims=2)

    # Reorder ranks (highest to lowest)
    s_order .= sortslices(s_order, dims=1, by=x -> x[2], rev=true)

    @views s_order[:, 3] .= Int.(1:size(S, 1))

    return s_order
end


"""
    topsis(S::Array{Float64, 2})

Calculates ranks using the aggregation method of the TOPSIS MCDA algorithm.
Rank for a particular site is calculated as a ratio

    C = S_n/(S_p + S_n)

S_n = √{∑(criteria .- NIS)²}
    is the squareroot of the summed differences between the criteria for a site and the
    Negative Ideal Solution (NIS), or worst performing site in each criteria.
S_p  = √{∑(criteria .- NIS)²}
    is the squareroot of the summed differences between the criteria for a site and the
    Positive Ideal Solution (PIS), or best performing site in each criteria.

    Details of this aggregation method in, for example [1].

# References
1. Opricovic, Serafim & Tzeng, Gwo-Hshiung. (2004) European Journal of Operational Research.
    Vol. 156. pp. 445.
    https://doi.org/10.1016/S0377-2217(03)00020-1.

# Arguments
- S : Decision matrix (seeding or shading)

# Returns
- s_order :
    1. site ids
    2. calculated site rank score (higher values = higher ranked)
    3. site order id
"""
function topsis(S::Array{Float64,2})::Array{Union{Float64,Int64},2}

    # compute the set of positive ideal solutions for each criteria
    PIS = maximum(S[:, 2:end], dims=1)

    # compute the set of negative ideal solutions for each criteria
    NIS = minimum(S[:, 2:end], dims=1)

    # calculate separation distance from the ideal and non-ideal solutions
    S_p = sqrt.(sum((S[:, 2:end] .- PIS) .^ 2, dims=2))
    S_n = sqrt.(sum((S[:, 2:end] .- NIS) .^ 2, dims=2))

    # final ranking measure of relative closeness C
    C = S_n ./ (S_p + S_n)

    # Create matrix where rank ids are integers (for use as indexers later)
    # Third column is derived from the number of sites for situations where
    # a subset of sites are being investigated (and so using their site IDs
    # will be inappropriate)
    s_order = Union{Float64,Int64}[Int.(S[:, 1]) C 1:size(S, 1)]

    # Reorder ranks (highest to lowest)
    s_order .= sortslices(s_order, dims=1, by=x -> x[2], rev=true)
    @views s_order[:, 3] .= Int.(1:size(S, 1))

    return s_order
end


"""
    vikor(S; v=0.5)

Calculates ranks using the aggregation method of the VIKOR MCDA algorithm.
Rank for a particular site is calculated as a linear combination of ratios,
weighted by v:
    Q = v(Sr - S_h) / (S_s - S_h) + (1 - v)(R - R_h) / (R_s - R_h)

where
- Sr = ∑(PIS-criteria) for each site, summed over criteria.
- R = max(PIS-criteria) for each site, with the max over criteria.
- S_h = min(∑(PIS-criteria)) over sites, the minimum summed distance from
    the positive ideal solution.
- S_s = max(∑(PIS-criteria)) over sites, maximum summed distance from
    the positive ideal solution.
- R_h = min(max(PIS-criteria)) over sites, the minimum max distance from
    the positive ideal solution.
- R_s = max(max(PIS-criteria)) over sites, the maximum max distance from
    the positive ideal solution.
- v = weighting, representing different decision-making strategies,
    or level of compromise between utility (overall performance)
    and regret (risk of performing very badly in one criteria despite
    exceptional performance in others)
    - v = 0.5 is consensus
    - v < 0.5 is minimal regret
    - v > 0.5 is max group utility (majority rules)

Details of this aggregation method in, for example [1]

# References
1. Alidrisi H. (2021) Journal of Risk and Financial Management.
    Vol. 14. No. 6. pp. 271.
    https://doi.org/10.3390/jrfm14060271

# Arguments
- S : Matrix
- v : Real

# Returns
- s_order :
    1. site ids
    2. calculated site rank score (higher values = higher ranked)
    3. site order id
"""
function vikor(S::Array{Float64,2}; v::Float64=0.5)::Array{Union{Float64,Int64},2}

    F_s = maximum(S[:, 2:end])

    # Compute utility of the majority Sr (Manhatten Distance)
    # Compute individual regret R (Chebyshev distance)
    sr_arg = (F_s .- S[:, 2:end])
    Sr = [S[:, 1] sum(sr_arg, dims=2)]
    R = [S[:, 1] maximum(sr_arg, dims=2)]

    # Compute the VIKOR compromise Q
    S_s, S_h = maximum(Sr[:, 2]), minimum(Sr[:, 2])
    R_s, R_h = maximum(R[:, 2]), minimum(R[:, 2])
    Q = @. v * (Sr[:, 2] - S_h) / (S_s - S_h) + (1 - v) * (R[:, 2] - R_h) / (R_s - R_h)
    Q .= 1.0 .- Q  # Invert rankings so higher values = higher rank

    # Create matrix where rank ids are integers (for use as indexers later)
    # Third column is necessary as a subset of sites will not match their Index IDs.
    s_order = Union{Float64,Int64}[Int.(S[:, 1]) Q zeros(size(Q, 1))]

    # sort Q by rank score in descending order
    s_order .= sortslices(s_order, dims=1, by=x -> x[2], rev=true)
    @views s_order[:, 3] .= Int.(1:size(S, 1))

    return s_order
end


"""
    unguided_site_selection(prefseedsites, prefshadesites, seed_years, shade_years, nsiteint, max_cover)

Randomly select seed/shade site locations for the given year, constraining to sites with max. carrying capacity > 0.
Here, `max_cover` represents the max. carrying capacity for each site (the `k` value).

# Arguments
- prefseedsites : Previously selected sites
- seed_years : bool, indicating whether to seed this year or not
- shade_years : bool, indicating whether to shade this year or not
- nsiteint : int, number of sites to intervene on
- available_space : vector/matrix : space available at each site (`k` value)
- depth : vector of site ids found to be within desired depth range
"""
function unguided_site_selection(prefseedsites, prefshadesites, seed_years, shade_years, nsiteint, available_space, depth)
    # Unguided deployment, seed/shade corals anywhere so long as available_space > 0.1
    # Only sites that have available space are considered, otherwise a zero-division error may occur later on.

    # Select sites (without replacement to avoid duplicate sites)
    candidate_sites = depth[(available_space.>0.0)[depth]]  # Filter down to site ids to be considered
    num_sites = length(candidate_sites)
    s_nsiteint = num_sites < nsiteint ? num_sites : nsiteint

    if seed_years
        prefseedsites = zeros(Int64, nsiteint)
        prefseedsites[1:s_nsiteint] .= StatsBase.sample(candidate_sites, s_nsiteint; replace=false)
    end

    if shade_years
        prefshadesites = zeros(Int64, nsiteint)
        prefshadesites[1:s_nsiteint] .= StatsBase.sample(candidate_sites, s_nsiteint; replace=false)
    end

    return prefseedsites[prefseedsites.>0], prefshadesites[prefshadesites.>0]
end
