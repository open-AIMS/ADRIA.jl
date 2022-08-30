"""Objects and methods for Dynamic Multi-Criteria Decision Analysis/Making"""


using StatsBase


struct DMCDA_vars  # {V, I, F, M} where V <: Vector
    site_ids  # ::V
    nsiteint  # ::I
    prioritysites  # ::V
    strongpred  # ::V
    centr  # ::V
    damprob  # ::A
    heatstressprob  # ::A
    sumcover  # ::F
    maxcover  # ::V
    area  # ::M
    risktol  # ::F
    wtconseed  # ::F
    wtconshade  # ::F
    wtwaves  # ::F
    wtheat  # ::F
    wthicover  # ::F
    wtlocover  # ::F
    wtpredecseed  # ::F
    wtpredecshade  # ::F
end


"""
    mcda_normalize(x::Union{Matrix, Vector})::Union{Matrix, Vector}

Normalize a Matrix (SE/SH) or Vector (wse/wsh) for MCDA.
"""
function mcda_normalize(x::Union{Matrix, Vector})::Union{Matrix, Vector}
    return x ./ sqrt(sum(x .^ 2))
end


"""
    align_rankings!(rankings::Array, s_order::Matrix, col::Int64)::Nothing

Align a vector of site rankings to match the indicated order in `s_order`.
"""
function align_rankings!(rankings::Array, s_order::Matrix, col::Int64)::Nothing
    # match site_ids by given order
    match_idx::Vector = findall(in.(rankings[:, 1], (s_order[:, 1], )))

    # Fill target ranking column
    rankings[match_idx, col] = s_order[:, 3]

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
"""
function rank_sites!(S, weights, rankings, nsiteint, mcda_func, rank_col)::Vector
    # Filter out all non-preferred sites
    selector = vec(.!all(S[:, 2:end] .== 0, dims=1))
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

    return prefsites
end
function rank_seed_sites!(S, weights, rankings, nsiteint, mcda_func)::Vector
    rank_sites!(S, weights, rankings, nsiteint, mcda_func, 2)
end
function rank_shade_sites!(S, weights, rankings, nsiteint, mcda_func)::Vector
    rank_sites!(S, weights, rankings, nsiteint, mcda_func, 3)
end

"""
    create_decision_matrix(site_ids, centr, sumcover, maxcover, area, damprob, heatstressprob, predec)

Creates criteria matrix `A`, where each column is a selection criterium and each row is a site.
Sites are then filtered based on heat and wave stress risk.

Where no sites are filtered, size of ``A := n_sites × 6 criteria``.

Columns indicate:
1. Site ID
2. Node Connectivity Centrality
3. Wave Damage Probability
4. Heat Stress Probability
5. Priority Predecessors
6. Available Area (relative to max cover)

# Arguments
- site_ids : vector of site ids
- centr : site centrality (relative strength of connectivity) (0 <= centr <= 1.0)
- sumcover : vector, sum of coral cover (across species) for each site (i.e., [x₁, x₂, ..., xₙ] where x_{1:n} <= 1.0)
- maxcover : maximum possible proportional coral cover (k) for each site, relative to total site area (k <= 1.0)
- area : total absolute area (in m²) for each site
- damprob : Probability of wave damage
- heatstressprob : Probability of site being affected by heat stress
- predec : list of priority predecessors (sites strongly connected to priority sites)
- risktol : tolerance for wave and heat risk (∈ [0,1]). Sites with heat or wave risk> risktol are filtered out.
"""
function create_decision_matrix(site_ids, centr, sumcover, maxcover, area, damprob, heatstressprob, predec, risktol)
    A = zeros(length(site_ids), 6)

    A[:, 1] .= site_ids  # Column of site ids

    # Account for cases where no coral cover
    c_cov_area = centr .* sumcover .* area

    # node connectivity centrality, need to instead work out strongest predecessors to priority sites
    A[:, 2] .= maximum(c_cov_area) != 0.0 ? c_cov_area / maximum(c_cov_area) : c_cov_area

    # Wave damage, account for cases where no chance of damage or heat stress
    # if max > 0 then use damage probability from wave exposure
    A[:, 3] .= maximum(damprob) != 0 ? (damprob .- minimum(damprob)) ./ (maximum(damprob) - minimum(damprob)) : damprob

    # risk from heat exposure
    A[:, 4] .= maximum(heatstressprob) != 0 ? (heatstressprob .- minimum(heatstressprob)) ./ (maximum(heatstressprob) - minimum(heatstressprob)) : heatstressprob

    # priority predecessors
    A[:, 5] .= predec[:, 3]

    # Proportion of empty space (no coral) compared to max possible cover
    A[:, 6] = (maxcover - sumcover) ./ maxcover

    # set any infs to zero
    A[maxcover .== 0, 6] .= 0.0

    # Filter out sites that have high risk of wave damage, specifically
    # exceeding the risk tolerance
    A[A[:, 3] .> risktol, 3] .= NaN
    rule = (A[:, 3] .<= risktol) .& (A[:, 4] .> risktol)
    A[rule, 4] .= NaN

    # remove rows with NaNs
    A = A[vec(.!any(isnan.(A), dims=2)), :]

    return A
end


"""
    create_seed_matrix(SE, A, wtconseed, wtwaves, wtheat, wtpredecseed, wtlocover)

Create seeding specific decision matrix from criteria matrix. The weight criteria and filter.

# Arguments
- SE : Seeding decision matrix, containing criteria specific to seeding
- A : Criteria  matrix
- wtconseed : Seed connectivity weight for seeding
- wtwaves : Wave stress weight
- wtheat : Heat stress weight
- wtpredecseed : Priority predecessor weight
- wtlocover : Weighting for low coral cover (coral real estate), when seeding 

# Returns
Tuple (SE, wse)
    - SE : Matrix of shape [n sites considered, 6]
           1. Site index ID
           2. Centrality
           3. Damage risk (higher values = less risk)
           4. Wave risk (higher values = less risk)
           5. Priority predecessors relating to coral real estate relative to max capacity
           6. Available space
    - wse : 5-element vector of criteria weights
        1. seed connectivity
        2. wave
        3. heat
        4. seed predecessors (?)
        5. low cover (?)
"""
function create_seed_matrix(SE, A, wtconseed, wtwaves, wtheat, wtpredecseed, wtlocover)
    wse = [1, wtconseed, wtwaves, wtheat, wtpredecseed, wtlocover]
    wse[2:end] .= mcda_normalize(wse[2:end])

    # Define seeding decision matrix
    SE[:, 1:2] .= A[:, 1:2]  # sites column (remaining), centrality

    SE[:, 3] .= 1.0 .- A[:, 3]  # compliment of damage risk
    SE[:, 4] .= 1.0 .- A[:, 4]  # compliment of wave risk
    SE[:, 5:6] .= A[:, 5:6]  # priority predecessors, coral real estate relative to max capacity

    # remove sites at maximum carrying capacity, take inverse log to emphasize importance of space for seeding
    SE = SE[vec(A[:, 6] .> 0), :]
    SE[:, 6] .= (10 .^ SE[:, 6]) ./ maximum(10 .^ SE[:, 6])

    return SE, wse
end


"""
    create_shade_matrix(SH, A, wtconshade, wtwaves, wtheat, wtpredecshade, wthicover)

Create shading specific decision matrix and apply weightings.

# Arguments
- SH : Shading decision matrix, containing criteria specific to seeding
- A : Criteria  matrix
- wtconshade : Shading connectivity weight
- wtwaves : Wave stress weight
- wtheat : Heat stress weight
- wtpredecshade : Priority predecessor weight for shading
- wthicover : Weighting for high coral cover when shading 

# Returns
Tuple (SH, wsh)
    - SH : Matrix of shape [n sites considered, 6]
           1. Site index ID
           2. Centrality
           3. Damage risk (higher values = less risk)
           4. Wave risk (higher values = less risk)
           5. Priority predecessors relating to coral real estate relative to max capacity
           6. Available space
    - wsh : 5-element vector of criteria weights
        1. shade connectivity
        2. wave
        3. heat
        4. shade predecessors (?)
        5. high cover (?)
"""
function create_shade_matrix(SH, A, wtconshade, wtwaves, wtheat, wtpredecshade, wthicover)
    wsh = [1, wtconshade, wtwaves, wtheat, wtpredecshade, wthicover]
    wsh[2:end] .= mcda_normalize(wsh[2:end])
    
    SH[:, 1:2] = A[:, 1:2] # sites column (remaining), absolute centrality
    SH[:, 3] = (1.0 .- A[:, 3]) # complimentary of wave damage risk
    SH[:, 4:5] = A[:, 4:5] # complimentary of heat damage risk, priority predecessors
    SH[:, 6] = (1.0 .- A[:, 6]) # coral cover relative to max capacity
    SH[SH[:,6] .> 1.0, 6] .= 1  # scale any sites above capacity back to 1
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
    site_ids::Array{Int64} = d_vars.site_ids
    nsites::Int64 = length(site_ids)
    nsiteint::Int64 = d_vars.nsiteint
    prioritysites::Array{Int64} = d_vars.prioritysites[in.(d_vars.prioritysites, [site_ids])]

    strongpred = d_vars.strongpred[site_ids, :]
    centr = d_vars.centr[site_ids]
    damprob = d_vars.damprob[site_ids]
    heatstressprob = d_vars.heatstressprob[site_ids]
    sumcover = d_vars.sumcover[site_ids]
    maxcover = d_vars.maxcover[site_ids]
    area = d_vars.area[site_ids]
    risktol = d_vars.risktol
    wtconseed = d_vars.wtconseed
    wtconshade = d_vars.wtconshade
    wtwaves = d_vars.wtwaves
    wtheat = d_vars.wtheat
    wthicover = d_vars.wthicover
    wtlocover = d_vars.wtlocover
    wtpredecseed = d_vars.wtpredecseed
    wtpredecshade = d_vars.wtpredecshade

    # site_id, seeding rank, shading rank
    rankings = Int64[site_ids zeros(Int64, nsites) zeros(Int64, nsites)]

    # work out which priority predecssors are connected to priority sites
    predec::Array{Float64} = zeros(nsites, 3)
    predec[:, 1:2] .= strongpred
    predprior = predec[in.(predec[:, 1], [prioritysites']), 2]
    predprior = [x for x in predprior if !isnan(x)]

    predec[predprior, 3] .= 1.0

    A = create_decision_matrix(site_ids, centr, sumcover, maxcover, area, damprob, heatstressprob, predec, risktol)
    if isempty(A)
        # if all rows have nans and A is empty, abort mission
        return prefseedsites, prefshadesites, rankings
    end

    # Set up SE and SH to be same size as A
    SE = zeros(size(A, 1), 6)
    SH = zeros(size(A, 1), 6)

    # cap to number of sites left after risk filtration
    nsiteint = min(nsiteint, length(A[:, 1]))

    # if seeding, create seeding specific decision matrix
    if log_seed
        SE, wse = create_seed_matrix(SE, A, wtconseed, wtwaves, wtheat, wtpredecseed, wtlocover)
    end

    # if shading, create shading specific decision matrix
    if log_shade
        SH, wsh = create_shade_matrix(SH, A, wtconshade, wtwaves, wtheat, wtpredecshade, wthicover)
    end

    if alg_ind == 1
        mcda_func = order_ranking
    elseif alg_ind == 2
        mcda_func = topsis
    elseif alg_ind == 3
        mcda_func = vikor
        # elseif alg_ind == 4
        #     mcda_func = multi_GA
    else
        error("Unknown MCDA algorithm selected. Valid options are 1 (Order Ranking), 2 (TOPSIS) and 3 (VIKOR).")
    end

    if isempty(SE)
        prefseedsites = repeat([0], nsiteint)
    elseif log_seed
        prefseedsites = rank_sites!(SE, wse, rankings, nsiteint, mcda_func)
    end

    if isempty(SH)
        prefshadesites = repeat([0], nsiteint)
    elseif log_shade
        prefshadesites = rank_sites!(SH, wsh, rankings, nsiteint, mcda_func)
    end

    # Replace with input rankings if seeding or shading rankings have not been filled
    if (sum(rankings[:, 2]) == 0.0) && (length(prefseedsites) != 0)
        rankings[:, 2] .= rankingsin[:, 2]
    end

    if (sum(rankings[:, 3]) == 0.0) && (length(prefshadesites) != 0)
        rankings[:, 3] .= rankingsin[:, 3]
    end

    return prefseedsites, prefshadesites, rankings
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
    2. calculated site rank
    3. site ranking order (lower values = higher ranked)
"""
function order_ranking(S::Array{Float64, 2})::Array{Union{Float64, Int64}, 2}
    n::Int64 = size(S,1)
    s_order::Array = Union{Float64, Int64}[zeros(Int, n) zeros(Float64, n) zeros(Int, n)]

    s_order[:, 3] .= Int.(1:size(S, 1))

    # simple ranking - add criteria weighted values for each sites
    @views s_order[:, 1] .= Int.(S[:, 1])
    @views s_order[:, 2] .= sum(S[:, 2:end], dims=2)
    s_order .= sortslices(s_order, dims=1, by=x->x[2], rev=true)

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
    2. calculated site rank
    3. site ranking order (lower values = higher ranked)
"""
function topsis(S::Array{Float64, 2})::Array{Union{Float64, Int64}, 2}

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
    s_order = Union{Float64, Int64}[Int.(S[:, 1]) C 1:size(S, 1)]

    # Reorder ranks
    s_order .= sortslices(s_order, dims=1, by=x->x[2], rev=true)

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
    2. calculated site rank
    3. site ranking order (lower values = higher ranked)
"""
function vikor(S::Array{Float64, 2}; v::Float64=0.5)::Array{Union{Float64, Int64}, 2}

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

    # Create matrix where rank ids are integers (for use as indexers later)
    s_order = Union{Float64, Int64}[Int.(S[:, 1]) Q Int.(1:size(Q, 1))]

    # sort Q in ascending order rows
    s_order .= sortslices(s_order, dims=1, by=x->x[2], rev=false)

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
"""
function unguided_site_selection(prefseedsites, prefshadesites, seed_years, shade_years, nsiteint, available_space)
    # Unguided deployment, seed/shade corals anywhere so long as available_space > 0.1
    # Only sites that have available space are considered, otherwise a zero-division error may occur later on.

    # Select sites (without replacement to avoid duplicate sites)
    candidate_sites = findall(available_space .> 0.0)
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

    return prefseedsites[prefseedsites .> 0], prefshadesites[prefshadesites .> 0]
end
