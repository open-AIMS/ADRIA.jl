"""
    order_ranking(S::Array{Float64, 2})

Uses simple summation as aggregation method for decision criteria.
Then orders sites from highest aggregate score to lowest.

# Arguments
- `S` : Decision matrix (seeding or shading)

# Returns
- `s_order` : nsites × 3 matrix with columns
    1. site ids
    2. calculated site rank score (higher values = higher ranked)
    3. site order id
"""
function order_ranking(S::Array{Float64,2})::Array{Union{Float64,Int64},2}
    return sum(S, dims=2)
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
- `S` : Decision matrix (seeding or shading)

# Returns
- `s_order` :
    1. site ids
    2. calculated site rank score (higher values = higher ranked)
    3. site order id
"""
function mcda_topsis(S::Array{Float64,2})::Array{Union{Float64,Int64},2}

    # compute the set of positive ideal solutions for each criteria
    PIS = maximum(S, dims=1)

    # compute the set of negative ideal solutions for each criteria
    NIS = minimum(S, dims=1)

    # calculate separation distance from the ideal and non-ideal solutions
    S_p = sqrt.(sum((S .- PIS) .^ 2, dims=2))
    S_n = sqrt.(sum((S .- NIS) .^ 2, dims=2))

    # final ranking measure of relative closeness C
    C = S_n ./ (S_p + S_n)

    return C
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
- `S` : Matrix
- `v` : Real

# Returns
- `s_order` :
    1. site ids
    2. calculated site rank score (higher values = higher ranked)
    3. site order id
"""
function mcda_vikor(S::Array{Float64,2}; v::Float64=0.5)::Array{Union{Float64,Int64},2}

    F_s = maximum(S)

    # Compute utility of the majority Sr (Manhatten Distance)
    # Compute individual regret R (Chebyshev distance)
    sr_arg = (F_s .- S)
    Sr = sum(sr_arg, dims=2)
    R = maximum(sr_arg, dims=2)

    # Compute the VIKOR compromise Q
    S_s, S_h = maximum(Sr), minimum(Sr)
    R_s, R_h = maximum(R), minimum(R)
    Q = @. v * (Sr - S_h) / (S_s - S_h) + (1 - v) * (R - R_h) / (R_s - R_h)
    Q .= 1.0 .- Q  # Invert rankings so higher values = higher rank

    return Q
end