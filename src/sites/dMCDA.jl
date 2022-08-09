"""Objects and methods for Dynamic Multi-Criteria Decision Analysis/Making"""

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
    rank_sites!(S, weights, rankings, nsiteint)

# Arguments
- S : Matrix, Site preference values
- weights : weights to apply
- rankings : vector of site ranks to update
- nsiteint : number of sites to select for interventions
"""
function rank_sites!(S, weights, rankings, nsiteint)::Vector
    # Filter out all non-preferred sites
    selector = vec(.!all(S .== 0, dims=1))
    weights = weights[selector]
    S = S[:, selector]

    S[:, 2:end] = mcda_normalize(S[:, 2:end])
    S .= S .* repeat(weights', size(S, 1), 1)
    s_order = mcda_func(S)

    last_idx = min(nsiteint, size(s_order, 1))
    prefshadesites = Int.(s_order[1:last_idx, 1])

    # Match by site_id and assign rankings to log
    align_rankings!(rankings, s_order, 3)

    return prefshadesites
end


"""
    create_decision_matrix(site_ids, centr, sumcover, maxcover, area, damprob, heatstressprob, predec)

Creates decision matrix `A`, with sites filtered based on risk.

If no sites are filtered the size of ``A := n_sites x 6``.

Columns indicate:
1. Site ID
2. Node Connectivity Centrality
3. Wave Damage Probability
4. Heat Stress Probability
5. Priority Predecessors
6. Available Area (relative to max cover)

# Arguments
- site_ids : vector of site ids
- centr : site centrality (relative strength of connectivity)
- sumcover : vector, sum of coral cover (across species) for each site (i.e., [x₁, x₂, ..., xₙ] where x_{1:n} <= 1.0)
- maxcover : maximum possible proportional coral cover (k) for each site, relative to total site area (k <= 1.0)
- area : absolute area (in m²) for each site
- damprob : Probability of wave damage
- heatstressprob : Probability of site being affected by heat stress
- predec : list of priority predecessors (sites strongly connected to priority sites)
- risktol : tolerance for wave and heat risk (0-1). Sites with heat or wave risk> risktol are filtered out.
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
    A[:, 3] .= maximum(damprob) != 0 ? damprob / maximum(damprob) : damprob

    # risk from heat exposure
    A[:, 4] .= maximum(heatstressprob) != 0 ? heatstressprob / maximum(heatstressprob) : heatstressprob

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
    filter_seed_sites(SE, A, wtconseed, wtwaves, wtheat, wtpredecseed, wtlocover)

Filter seed sites

# Arguments
- SE : Pre-existing seed decision matrix to populate
- A : Criteria decision matrix
- wtconseed : Seed connectivity weight
- wtwaves : Wave stress weight
- wtheat : heat stress weight
- wtpredecseed : priority predecessor weight
- wtlocover : weighting for low coral cover (coral real estate), when seeding
- wthicover : weighting for hight coral cover, when shading
"""
function filter_seed_sites(SE, A, wtconseed, wtwaves, wtheat, wtpredecseed, wtlocover)
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

    return SE
end


"""
    dMCDA(d_vars::DMCDA_vars, alg_ind::Int64, log_seed::Bool, log_shade::Bool, prefseedsites::AbstractArray{Int}, prefshadesites::AbstractArray{Int}, rankingsin::Matrix{Int64})

# Returns
Tuple : preferred seed sites, preferred shade/fog sites, number of seed sites, number of shade sites, rankings

        `rankings` is an Nx3 matrix holding: site_id, seeding_rank, shading_rank
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

    ## Seeding - Filtered set
    # define seeding weights
    if log_seed
        SE = filter_seed_sites(SE, A, wtconseed, wtwaves, wtheat, wtpredecseed, wtlocover)
    end

    if log_shade
        ## Shading filtered set
        # define shading weights
        wsh = [1, wtconshade, wtwaves, wtheat, wtpredecshade, wthicover]
        wsh[2:end] .= mcda_normalize(wsh[2:end])

        SH[:, 1:2] = A[:, 1:2] # sites column (remaining), absolute centrality
        SH[:, 3] = 1.0 .- A[:, 3] # complimentary of wave damage risk
        SH[:, 4:5] = A[:, 4:5] # complimentary of heat damage risk, priority predecessors
        SH[:, 6] = 1.0 .- A[:, 6] # coral cover relative to max capacity
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
        prefseedsites = rank_sites!(SE, wse, rankings, nsiteint)
    end

    if isempty(SH)
        prefshadesites = repeat([0], nsiteint)
    elseif log_shade
        prefshadesites = rank_sites!(SH, wsh, rankings, nsiteint)
    end

    # Replace with input rankings if seeding or shading rankings have not been filled
    if (sum(rankings[:, 2]) == 0.0) && (length(prefseedsites) != 0)
        rankings[:, 2] .= @view rankingsin[:, 2]
    end

    if (sum(rankings[:, 3]) == 0.0) && (length(prefshadesites) != 0)
        rankings[:, 3] .= @view rankingsin[:, 3]
    end

    return prefseedsites, prefshadesites, rankings
end


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


function topsis(S::Array{Float64, 2})::Array{Union{Float64, Int64}, 2}

    # compute the set of positive ideal solutions for each criteria (max for
    # good criteria, min for bad criteria). Max used as all criteria
    # represent preferred attributes not costs or negative attributes
    PIS = maximum(S[:, 2:end], dims=1)

    # compute the set of negative ideal solutions for each criteria
    # (min for good criteria, max for bad criteria).
    # Min used as all criteria represent preferred attributes not
    # costs or negative attributes
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

# Arguments
- S : Matrix
- v : Real, level of compromise (utility vs. regret).
        - v = 0.5 is consensus
        - v < 0.5 is minimal regret
        - v > 0.5 is max group utility (majority rules)

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
    unguided_site_selection!(prefseedsites, prefshadesites, seed_years, shade_years, nsiteint, max_cover)

Randomly select seed/shade site locations for the given year, constraining to sites with max. carrying capacity > 0.
Here, `max_cover` represents the max. carrying capacity for each site (the `k` value).

# Arguments
- prefseedsites : Previously selected sites
- seed_years : bool, indicating whether to seed this year or not
- shade_years : bool, indicating whether to shade this year or not
- nsiteint : int, number of sites to intervene on
- max_cover : vector/matrix : maximum carrying capacity of each site (`k` value)
"""
function unguided_site_selection!(prefseedsites, prefshadesites, seed_years, shade_years, nsiteint, max_cover)::Nothing
    # Unguided deployment, seed/shade corals anywhere so long as max_cover > 0
    # otherwise, may select an empty site, causing zero-division error later on.

    # `unique()` used to catch cases where number of eligible sites < `nsiteint`
    # `resize!()` to recreate vector of `nsiteint` entries if duplicate sites are filtered out
    if seed_years
        prefseedsites .= resize!(unique(rand(findall(max_cover .> 0.0), nsiteint)), nsiteint)
    end

    if shade_years
        prefshadesites .= resize!(unique(rand(findall(max_cover .> 0.0), nsiteint)), nsiteint)
    end

    return nothing
end
