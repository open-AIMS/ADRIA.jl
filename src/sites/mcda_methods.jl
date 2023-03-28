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
    n::Int64 = size(S, 1)
    s_order::Array = Union{Float64,Int64}[zeros(Int, n) zeros(Float64, n) zeros(Int, n)]

    # Simple ranking - add criteria weighted values for each sites
    # Third column is derived from the number of sites for situations where
    # a subset of sites are being investigated (and so using their site IDs
    # will be inappropriate)
    @views s_order[:, 1] .= Int.(S[:, 1])
    @views s_order[:, 2] .= sum(S[:, 2:end], dims=2)

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
- `S` : Decision matrix (seeding or shading)

# Returns
- `s_order` :
    1. site ids
    2. calculated site rank score (higher values = higher ranked)
    3. site order id
"""
function mcda_topsis(S::Array{Float64,2})::Array{Union{Float64,Int64},2}

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
- `S` : Matrix
- `v` : Real

# Returns
- `s_order` :
    1. site ids
    2. calculated site rank score (higher values = higher ranked)
    3. site order id
"""
function mcda_vikor(S::Array{Float64,2}; v::Float64=0.5)::Array{Union{Float64,Int64},2}

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
    return s_order
end