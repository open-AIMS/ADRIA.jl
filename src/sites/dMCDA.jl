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
    mcda_normalize(x)

Normalize a Matrix (SE/SH) or Vector (wse/wsh) for MCDA.
"""
function mcda_normalize(x::Union{Matrix, Vector})::Union{Matrix, Vector}
    return x ./ sqrt(sum(x .^ 2))
end


"""
"""
function align_rankings!(rankings::Array, s_order::Matrix, col::Int64)::Nothing
    # Add ranking column
    # s_order[:, 3] = 1:length(s_order[:, 1])

    # match site_ids by given order
    match_idx = findall(in.(rankings[:, 1], (s_order[:, 1], )))

    # Fill target ranking column
    rankings[match_idx, col] = s_order[:, 3]

    return
end


"""
prefseedsites, prefshadesites, nprefseedsites, nprefshadesites, rankings

# Returns
Tuple : preferred seed sites, preferred shade/fog sites, number of seed sites, number of shade sites, rankings

        `rankings` is an Nx3 matrix holding: site_id, seeding_rank, shading_rank
        0 indicates sites that were not considered
"""
function dMCDA(d_vars, alg_ind, log_seed, log_shade, prefseedsites, prefshadesites, rankingsin)

    site_ids = d_vars.site_ids
    nsites = length(site_ids)
    nsiteint = d_vars.nsiteint
    prioritysites = d_vars.prioritysites[in.(d_vars.prioritysites, [site_ids])]

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
    rankings = [site_ids zeros(Int, nsites) zeros(Int, nsites)]

    predec = zeros(nsites, 3)
    predec[:, 1:2] .= strongpred
    predprior = predec[in.(predec[:, 1], [prioritysites']), 2]
    predprior = [x for x in predprior if !isnan(x)]

    predec[predprior, 3] .= 1

    # Combine data into matrix
    A = zeros(length(site_ids), 6)
    SE = zeros(length(site_ids), 6)
    SH = zeros(length(site_ids), 6)

    A[:, 1] = site_ids  # column of site IDs

    # Account for cases where no coral cover
    c_cov_area = centr .* sumcover .* area

    # node connectivity centrality, need to instead work out strongest predecessors to priority sites
    A[:, 2] = maximum(c_cov_area) != 0.0 ? c_cov_area / maximum(c_cov_area) : c_cov_area

    # Account for cases where no chance of damage or heat stress
    # if max > 0 then use damage probability from wave exposure
    A[:, 3] .= maximum(damprob) != 0 ? damprob / maximum(damprob) : damprob

    # risk from heat exposure
    A[:, 4] = maximum(heatstressprob) != 0 ? heatstressprob / maximum(heatstressprob) : heatstressprob

    # priority predecessors
    A[:, 5] = predec[:, 3]

    A[:, 6] = (maxcover - sumcover) ./ maxcover # proportion of cover compared to max possible cover

    # set any infs to zero
    A[maxcover .== 0, 6] .= 0.0

    # Filter out sites that have high risk of wave damage, specifically
    # exceeding the risk tolerance
    A[A[:, 3] .> risktol, 3] .= NaN
    rule = (A[:, 3] .<= risktol) .& (A[:, 4] .> risktol)
    A[rule, 4] .= NaN

    # remove rows with NaNs
    A .= A[vec(.!any(isnan.(A), dims=2)), :]

    if isempty(A)
        # if all rows have nans and A is empty, abort mission
        nprefseedsites = 0
        nprefshadesites = 0
        return prefseedsites, nprefseedsites, prefshadesites, nprefshadesites, rankings
    end

    # number of sites left after risk filtration
    if nsiteint > length(A[:, 1])
        nsiteint = length(A[:, 1])
    end

    ## Seeding - Filtered set
    # define seeding weights
    if log_seed
        wse = [1, wtconseed, wtwaves, wtheat, wtpredecseed, wtlocover]
        wse[2:end] .= mcda_normalize(wse[2:end])

        # define seeding decision matrix
        SE[:, 1] = A[:, 1]  # sites column (remaining)
        SE[:, 2] = A[:, 2]  # centrality
        SE[:, 3] = (1.0 - A[:, 3])  # complementary of damage risk
        SE[:, 4] = (1.0 - A[:, 4])  # complimetary of wave risk
        SE[:, 5] = A[:, 5]  # priority predecessors
        SE[:, 6] = A[:, 6]  # coral real estate relative to max capacity

        # remove sites at maximum carrying capacity
        SE .= SE[vec(A[:, 6] .<= 0), :]
    end

    if log_shade
        ## Shading filtered set
        # define shading weights
        wsh = [1, wtconshade, wtwaves, wtheat, wtpredecshade, wthicover]
        wsh[2:end] .= mcda_normalize(wsh[2:end])

        SH[:, 1] = A[:, 1] # sites column (remaining)
        SH[:, 2] = A[:, 2] # absolute centrality
        SH[:, 3] = (1.0 .- A[:, 3]) # complimentary of wave damage risk
        SH[:, 4] = A[:, 4] # complimentary of heat damage risk
        SH[:, 5] = A[:, 5] # priority predecessors
        SH[:, 6] = (1.0 .- A[:, 6]) # coral cover relative to max capacity
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
        error("Unknown MCDA algorithm selected. Valid options are 1 to 3.")
    end

    if isempty(SE)
        prefseedsites = 0
    elseif log_seed
        # Remove cols that are all 0
        selector = vec(.!all(SE .== 0, dims=1))
        wse = wse[selector]
        SE = SE[:, selector]

        # normalisation
        SE[:, 2:end] = mcda_normalize(SE[:, 2:end])
        SE .= SE .* repeat(wse', size(SE, 1), 1)
        s_order = mcda_func(SE)

        last_idx = min(nsiteint, size(s_order, 1))
        prefseedsites = s_order[1:last_idx, 1]

        # Match by site_id and assign rankings to log
        align_rankings!(rankings, s_order, 2)
    end

    if isempty(SH)
        prefshadesites = 0
    elseif log_shade
        # Remove cols that are all 0
        selector = vec(.!all(SH .== 0, dims=1))
        wsh = wsh[selector]
        SH = SH[:, selector]

        # normalisation
        SH[:, 2:end] = mcda_normalize(SH[:, 2:end])
        SH .= SH .* repeat(wsh', size(SH, 1), 1)
        s_order = mcda_func(SH)

        last_idx = min(nsiteint, size(s_order, 1))
        prefshadesites = s_order[1:last_idx, 1]

        # Match by site_id and assign rankings to log
        align_rankings!(rankings, s_order, 3)
    end

    nprefseedsites = length(prefseedsites)
    nprefshadesites = length(prefshadesites)

    # Replace with input rankings if seeding or shading rankings have not been filled
    if (sum(rankings[:, 2]) == 0.0) && (nprefseedsites != 0)
        rankings[:, 2] .= rankingsin[:, 2]
    end

    if (sum(rankings[:, 3]) == 0.0) && (nprefshadesites != 0)
        rankings[:, 3] .= rankingsin[:, 3]
    end

    return prefseedsites, prefshadesites, nprefseedsites, nprefshadesites, rankings
end


function order_ranking(S)
    n = size(S,1)
    s_order = Union{Float64, Int64}[zeros(Int, n) zeros(Float64, n) zeros(Int, n)]

    s_order[:, 3] .= Int.(1:size(S, 1))

    # simple ranking - add criteria weighted values for each sites
    @views s_order[:, 1] .= Int.(S[:, 1])
    @views s_order[:, 2] .= sum(S[:, 2:end], dims=2)
    s_order .= sortslices(s_order, dims=1, by=x->x[2], rev=true)

    return s_order
end


function topsis(S)

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

Parameters
----------
S : Matrix
v : Real, level of compromise (utility vs. regret).
        - v = 0.5 is consensus
        - v < 0.5 is minimal regret
        - v > 0.5 is max group utility (majority rules)

"""
function vikor(S; v=0.5)

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
