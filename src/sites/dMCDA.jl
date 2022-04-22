struct DMCDA_vars{A, I, F}
    site_ids::A
    nsiteint::I
    prioritysites::A
    strongpred::A
    centr::A
    damprob::A
    heatstressprob::A
    sumcover::F
    maxcover::F
    area::A
    risktol::F
    wtconseed::F
    wtconshade::F
    wtwaves::F
    wtheat::F
    wthicover::F
    wtlocover::F
    wtpredecseed::F
    wtpredecshade::F

    # dMCDA_vars = struct('site_ids', depth_priority, 'nsiteint', nsiteint, 'prioritysites', sim_params.prioritysites, ...
    #     'strongpred', strongpred, 'centr', site_ranks.C1, 'damprob', 0, 'heatstressprob', 0, ...
    #     'sumcover', 0, 'maxcover', max_cover, 'area', site_data.area, 'risktol', risktol, 'wtconseed', wtconseed, 'wtconshade', wtconshade, ...
    #     'wtwaves', wtwaves, 'wtheat', wtheat, 'wthicover', wthicover, 'wtlocover', wtlocover, 'wtpredecseed', wtpredecseed, 'wtpredecshade', wtpredecshade);
end

mcda_normalize = (x) -> x[:, 2:end] ./ sqrt(sum(x[:, 2:end] .^ 2))


function align_rankings!(rankings::Array, s_order::Array)::Nothing
    # Add ranking column
    s_order[:, 3] = 1:length(s_order[:, 1])

    # [~, ii] = ismember(s_order[:, 1], rankings[:, 1], "rows")
    # align = ii[ii.!=0]
    align = rankings[:, 1] .== eachrow(s_order[:, 1])
    rankings[align, 2] .= s_order[:, 3]

    return
end


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
    rankings = [site_ids, zeros(nsites, 1), zeros(nsites, 1)]

    predec = zeros(nsites, 3)
    predec[:, 1:2] .= strongpred
    predprior = predec[in.(predec[:, 1], [prioritysites']), 2]
    deleteat!(predprior, findall(isnan.(predprior)))
    predec[predprior, 3] = 1

    # Combine data into matrix
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
    A[maxcover==0, 6] = 0.0

    # Filter out sites that have high risk of wave damage, specifically
    # exceeding the risk tolerance
    A[A[:, 3].>risktol, 3] = NaN
    rule = (A[:, 3] .<= risktol) & (A[:, 4] .> risktol)
    A[rule, 4] = NaN

    A .= A[vec(.!any(isnan.(A), dims=2)), :]  # if a row has a nan, delete it

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
        # wse[2:end] .= wse[2:end] ./ sqrt(sum(wse[2:end].^2));
        wse[2:end] .= mcda_normalize(wse)

        # define seeding decision matrix
        SE[:, 1] = A[:, 1]  # sites column (remaining)
        SE[:, 2] = A[:, 2]  # centrality
        SE[:, 3] = (1.0 - A[:, 3])  # complementary of damage risk
        SE[:, 4] = (1.0 - A[:, 4])  # complimetary of wave risk
        SE[:, 5] = A[:, 5]  # priority predecessors
        SE[:, 6] = A[:, 6]  # coral real estate relative to max capacity
        # SE[A[:, 6].<=0, :] = []  # remove sites at maximum carrying capacity
        deleteat!(SE, A[:, 6] .<= 0)
    end

    if log_shade
        ## Shading filtered set
        # define shading weights
        wsh = [1, wtconshade, wtwaves, wtheat, wtpredecshade, wthicover]
        # wsh[2:end] = wsh[2:end] ./ sqrt(sum(wsh[2:end].^2));
        wsh[2:end] .= mcda_normalize(wsh)

        SH[:, 1] = A[:, 1] # sites column (remaining)
        SH[:, 2] = A[:, 2] # absolute centrality
        SH[:, 3] = (1.0 - A[:, 3]) # complimentary of wave damage risk
        SH[:, 4] = A[:, 4] # complimentary of heat damage risk
        SH[:, 5] = A[:, 5] # priority predecessors
        SH[:, 6] = (1.0 - A[:, 6]) # coral cover relative to max capacity
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
        selector = all(SE .== 0, dims=1)
        deleteat!(wse, selector)
        SE = SE[:, .!vec(selector)]

        # normalisation
        SE[:, 2:end] = mcda_normalize(SE)
        SE .= SE .* repeat(wse, size(SE, 1), 1)
        s_order = mcda_func(SE)

        last_idx = min(nsiteint, height(s_order))
        prefseedsites = s_order[1:last_idx, 1]

        # Match by site_id and assign rankings to log
        # [~, ii] = ismember(s_order[:, 1], rankings[:, 1], "rows")
        # align = ii[ii.!=0]
        align_rankings!(rankings, s_order)
    end

    if isempty(SH)
        prefshadesites = 0
    elseif log_shade
        # Remove cols that are all 0
        selector = all(SH .== 0, dims=1)
        deleteat!(wsh, selector)
        SH = SH[:, .!vec(selector)]

        # normalisation
        SH[:, 2:end] = mcda_normalize(SH)
        SH .= SH .* repeat(wsh, size(SH, 1), 1)
        s_order = mcda_func(SH)

        last_idx = min(nsiteint, height(s_order))
        prefshadesites = s_order[1:last_idx, 1]

        # [~, ii] = ismember(s_order[:, 1], rankings[:, 1], "rows")
        # align = ii[ii.!=0]
        align_rankings!(rankings, s_order)
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

    s_order = zeros(nrow(S), 3)
    s_order[:, 3] .= 1:nrow(S)

    # simple ranking - add criteria weighted values for each sites
    @views @. s_order[:, 1] = S[:, 1]
    @views @. s_order[:, 2] = sum(S[:, 2:end], 2)
    sort!(s_order, dims=2, rev=true)

    return s_order
end


function topsis(S)

    # compute the set of positive ideal solutions for each criteria (max for
    # good crieteria, min for bad criteria). Max used as all crieteria
    # represent preferred attributes not costs or negative attributes
    PIS = maximum(S[:, 2:end])

    # compute the set of negative ideal solutions for each criteria
    # (min for good criteria, max for bad criteria).
    # Min used as all criteria represent preferred attributes not
    # costs or negative attributes
    NIS = minimum(S[:, 2:end])

    # calculate separation distance from the ideal and non-ideal solns
    S_p = sqrt(sum((S[:, 2:end] - PIS) .^ 2, 2))
    S_n = sqrt(sum((S[:, 2:end] - NIS) .^ 2, 2))

    # final ranking measure of relative closeness C
    C = S_n ./ (S_p + S_n)
    S_wt = [S[:, 1], C]
    s_order = sort(S_wt, dims=2, rev=true)

    return s_order
end


function vikor(S)
    F_s = max(S[:, 2:end])

    # Compute utility of the majority Sr (Manhatten Distance)
    # Compute individual regret R (Chebyshev distance)
    sr_arg = (F_s - S[:, 2:end])
    Sr = sum(sr_arg, 2)
    Sr = [S[:, 1], Sr]

    R = max(sr_arg, [], 2)  # TODO: Fix this matlab max() function use
    R = [S[:, 1], R]

    # Compute the VIKOR compromise Q
    S_s = max(Sr[:, 2])
    S_h = min(Sr[:, 2])
    R_s = max(R[:, 2])
    R_h = min(R[:, 2])
    Q = v * (Sr[:, 2] - S_h) / (S_s - S_h) + (1 - v) * (R[:, 2] - R_h) / (R_s - R_h)
    Q = [S[:, 1], Q]

    # sort Q in ascending order rows
    s_order = sort(Q, dims=2)

    return s_order
end
