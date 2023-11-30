module decision
"""Objects and methods for Dynamic Multi-Criteria Decision Analysis/Making"""

using StatsBase
using Distances, DataFrames
using Combinatorics
using JMcDM
using InteractiveUtils: subtypes
using ADRIA: Domain, EcoModel, n_locations, site_area, site_k_area, site_k

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
    leftover_space  # ::F
    k_area  # ::V
    min_area # ::F
    risk_tol  # ::F
    use_spatial_group # ::M
    n_spatial_grp # ::F
    reefs #::V
    wt_in_conn_seed  # ::F
    wt_out_conn_seed  # ::F
    wt_conn_fog  # ::F
    wt_waves_seed # ::F
    wt_waves_fog # ::F
    wt_heat_seed  # ::F
    wt_heat_fog  # ::F
    wt_depth_seed # ::F
    wt_hi_cover  # ::F
    wt_lo_cover  # ::F
    wt_predec_seed  # ::F
    wt_predec_fog  # ::F
    wt_zones_seed # ::F
    wt_zones_fog # ::F
end

include("mcda_methods.jl")
include("location_selection.jl")

function supported_jmcdm_methods()
    return [
        JMcDM.ARAS.ArasMethod,
        JMcDM.COCOSO.CocosoMethod,
        JMcDM.CODAS.CodasMethod,
        JMcDM.EDAS.EdasMethod,
        JMcDM.GREY.GreyMethod,
        JMcDM.MABAC.MabacMethod,
        JMcDM.MAIRCA.MaircaMethod,
        JMcDM.MARCOS.MarcosMethod,
        JMcDM.MOORA.MooraMethod,
        JMcDM.PIV.PIVMethod,
        JMcDM.PSI.PSIMethod,
        JMcDM.SAW.SawMethod,
        JMcDM.WASPAS.WaspasMethod,
        JMcDM.WPM.WPMMethod
    ]
end

function mcda_methods()
    return [
        order_ranking,
        adria_vikor,
        adria_topsis,
        supported_jmcdm_methods()...
    ]
end

"""
    DMCDA_vars(domain::Domain, criteria::NamedDimsArray,
               site_ids::AbstractArray, leftover_space::AbstractArray, area_to_seed::Float64,
               waves::AbstractArray, dhws::AbstractArray)::DMCDA_vars
    DMCDA_vars(domain::Domain, criteria::NamedDimsArray, site_ids::AbstractArray,
                leftover_space::AbstractArray, area_to_seed::Float64)::DMCDA_vars
    DMCDA_vars(domain::Domain, criteria::DataFrameRow, site_ids::AbstractArray,
                leftover_space::AbstractArray, area_to_seed::Float64)::DMCDA_vars
    DMCDA_vars(domain::Domain, criteria::DataFrameRow, site_ids::AbstractArray,
                leftover_space::AbstractArray, area_to_seed::Float64,
               waves::AbstractArray, dhw::AbstractArray)::DMCDA_vars

Constuctors for DMCDA variables.
"""
function DMCDA_vars(
    domain::Domain,
    criteria::NamedDimsArray,
    site_ids::AbstractArray,
    leftover_space::AbstractArray,
    area_to_seed::Float64,
    waves::AbstractArray,
    dhws::AbstractArray
)::DMCDA_vars

    # Site Data
    site_d = domain.site_data

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
        leftover_space,
        site_k_area(domain),
        criteria("coral_cover_tol") .* area_to_seed,
        criteria("deployed_coral_risk_tol"),
        criteria("use_spatial_group"),
        domain.sim_constants.n_spatial_grp,
        domain.site_data.Reef,
        criteria("seed_in_connectivity"),
        criteria("seed_out_connectivity"),
        criteria("fog_connectivity"),
        criteria("seed_wave_stress"),
        criteria("fog_wave_stress"),
        criteria("seed_heat_stress"),
        criteria("fog_heat_stress"),
        criteria("seed_depth"),
        criteria("coral_cover_high"),
        criteria("coral_cover_low"),
        criteria("seed_priority"),
        criteria("fog_priority"),
        criteria("seed_zone"),
        criteria("fog_zone"),
    )

    return mcda_vars
end
function DMCDA_vars(
    domain::Domain,
    criteria::NamedDimsArray,
    site_ids::AbstractArray,
    leftover_space::AbstractArray,
    area_to_seed::Float64
)::DMCDA_vars
    num_sites = n_locations(domain)
    return DMCDA_vars(
        domain,
        criteria,
        site_ids,
        leftover_space,
        area_to_seed,
        zeros(num_sites, 1),
        zeros(num_sites, 1),
    )
end
function DMCDA_vars(
    domain::Domain,
    criteria::DataFrameRow,
    site_ids::AbstractArray,
    leftover_space::AbstractArray,
    area_to_seed::Float64,
    waves::AbstractArray,
    dhw::AbstractArray
)::DMCDA_vars
    criteria_vec::NamedDimsArray = NamedDimsArray(collect(criteria), rows=names(criteria))
    return DMCDA_vars(
        domain, criteria_vec, site_ids, leftover_space, area_to_seed, waves, dhw
    )
end
function DMCDA_vars(
    domain::Domain,
    criteria::DataFrameRow,
    site_ids::AbstractArray,
    leftover_space::AbstractArray,
    area_to_seed::Float64
)::DMCDA_vars

    criteria_vec::NamedDimsArray = NamedDimsArray(collect(criteria), rows=names(criteria))
    return DMCDA_vars(domain, criteria_vec, site_ids, leftover_space, area_to_seed)
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
    mcda_normalize(x::DataFrame)::DataFrame

Normalize weights for a set of scenarios (wse/wsh) for MCDA.
"""
function mcda_normalize(x::DataFrame)::DataFrame
    return x ./ sum(Matrix(x); dims=2)
end

"""
    align_rankings!(rankings::Array, s_order::Matrix, col::Int64)::Nothing

Align a vector of site rankings to match the indicated order in `s_order`.
"""
function align_rankings!(rankings::Array, s_order::Matrix, col::Int64)::Nothing
    # Fill target ranking column
    for (i, site_id) in enumerate(s_order[:, 1])
        rankings[rankings[:, 1].==site_id, col] .= i
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
- `rank_col` : column to fill with rankings (2 for seed, 3 for fog)

# Returns
Sites in order of their rankings
"""
function rank_sites!(
    S::Matrix{Float64},
    weights::Vector{Float64},
    rankings::Matrix{Int64},
    n_site_int::Int64,
    mcda_func::Union{Function,Type{<:MCDMMethod}},
    rank_col)::Tuple{Vector{Int64},Matrix{Union{Float64,Int64}}}
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
    retrieve_ranks(S::Matrix, site_ids::Vector, weights::Vector{Float64}, mcda_func::Function)::Matrix{Union{Float64,Int64}}
    retrieve_ranks(S::Matrix, site_ids::Vector, weights::Vector{Float64}, mcda_func::Type{<:MCDMMethod})::Matrix{Union{Float64,Int64}}
    retrieve_ranks(site_ids::Vector, scores::Vector, maximize::Bool)::Matrix{Union{Float64,Int64}}

Get location ranks using mcda technique specified in mcda_func, weights and a decision matrix S.

# Arguments
- `S` : decision matrix containing criteria values for each location (n locations)*(m criteria)
- `site_ids` : array of site ids still remaining after filtering.
- `weights` : importance weights for each criteria.
- `mcda_func` : function/JMcDM DataType to use for mcda, specified as an element from methods_mcda.
- `scores` : set of scores derived from applying an mcda ranking method.
- `maximize` : Boolean indicating whether a mcda method is maximizing score (true), or minimizing (false).

# Returns
Matrix, containing site_ids, criteria values, ranks
"""
function retrieve_ranks(
    S::Matrix{Float64},
    site_ids::Vector{Float64},
    weights::Vector{Float64},
    mcda_func::Function)::Matrix{Union{Float64,Int64}}
    S = mcda_normalize(S) .* weights'
    scores = mcda_func(S)

    return retrieve_ranks(site_ids, vec(scores), true)
end
function retrieve_ranks(
    S::Matrix{Float64},
    site_ids::Vector{Float64},
    weights::Vector{Float64},
    mcda_func::Type{<:MCDMMethod},
)::Matrix{Union{Float64,Int64}}
    fns = fill(maximum, length(weights))
    results = mcdm(MCDMSetting(S, weights, fns), mcda_func())
    maximize = results.bestIndex == argmax(results.scores)

    return retrieve_ranks(site_ids, results.scores, maximize)
end
function retrieve_ranks(
    site_ids::Vector{Float64},
    scores::Vector,
    maximize::Bool
)::Matrix{Union{Float64,Int64}}
    s_order::Vector{Int64} = sortperm(scores, rev=maximize)
    return Union{Float64,Int64}[Int64.(site_ids[s_order]) scores[s_order]]
end


"""
    create_decision_matrix(site_ids, in_conn, out_conn, leftover_space, wave_stress, heat_stress, predec, risk_tol)

Creates criteria matrix `A`, where each column is a selection criterium and each row is a site.
Sites are then filtered based on heat and wave stress risk.

Where no sites are filtered, size of ``A := n_sites × 9 criteria``.

Columns indicate:
1. Site ID
2. Incoming node connectivity centrality
3. Outgoing node connectivity centrality
4. Wave stress
5. Heat stress
6. Priority predecessors
7. GBRMPA zoning criteria
8. Available area (relative to max cover)
9. Location depth

# Arguments
- `site_ids` : Unique site ids.
- `in_conn` : Site incoming centrality (relative strength of connectivity) (0 <= c <= 1.0).
- `out_conn` : Site outgoing centrality (relative strength of connectivity) (0 <= c <= 1.0).
- `leftover_space` : Available space for coral to grow (m²).
- `wave_stress` : Probability of wave damage.
- `heat_stress` : Probability of site being affected by heat stress
- `predec` : List of priority predecessors (sites strongly connected to priority sites).
- `risk_tol` : Tolerance for wave and heat risk (∈ [0,1]). Sites with heat or wave risk> risk_tol are filtered out.
"""
function create_decision_matrix(
    site_ids::Vector{Int64},
    in_conn::T,
    out_conn::T,
    leftover_space::Union{NamedDimsArray,T},
    wave_stress::T,
    heat_stress::T,
    site_depth::T,
    predec::Matrix{Float64},
    zones_criteria::T,
    risk_tol::Float64
)::Tuple{Matrix{Float64}, BitVector} where {T<:Vector{Float64}}
    A = zeros(length(site_ids), 9)
    A[:, 1] .= site_ids  # Column of site ids

    # node connectivity centrality, need to instead work out strongest predecessors to priority sites
    A[:, 2] .= maximum(in_conn) != 0.0 ? in_conn / maximum(in_conn) : 0.0
    A[:, 3] .= maximum(out_conn) != 0.0 ? out_conn / maximum(out_conn) : 0.0

    # Wave damage, account for cases where no chance of damage or heat stress
    # if max > 0 then use damage probability from wave exposure
    max_ws = maximum(wave_stress)
    min_ws = minimum(wave_stress)
    A[:, 4] .= max_ws == 0.0 ? 0.0 : (wave_stress .- min_ws) ./ (max_ws - min_ws)

    # risk from heat exposure
    max_hs = maximum(heat_stress)
    min_hs = minimum(heat_stress)
    A[:, 5] .= max_hs == 0.0 ? 0.0 : (heat_stress .- min_hs) ./ (max_hs - min_hs)

    # priority predecessors
    A[:, 6] .= predec[:, 3]

    # priority zone predecessors and sites
    A[:, 7] .= zones_criteria

    # Area of space for coral cover in m²
    A[:, 8] = leftover_space

    A[:, 9] = site_depth

    # Filter out sites that have high risk of wave damage, specifically
    # exceeding the risk tolerance
    A[A[:, 4].>risk_tol, 4] .= NaN
    rule = (A[:, 4] .<= risk_tol) .& (A[:, 5] .> risk_tol)
    A[rule, 5] .= NaN

    filtered = vec(.!any(isnan.(A), dims=2))

    # Remove rows with NaNs
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
- `wt_depth_seed` : Median depth weight
- `wt_predec_seed` : Priority predecessor weight
- `wt_predec_zones_seed` : Priority zones weight for seeding
- `wt_low_cover` : Weighting for low coral cover (coral real estate), when seeding

# Returns
Tuple (SE, wse)
- `SE` : Matrix of shape [n sites considered, 7]
    1. Site index ID
    2. Incoming centrality
    3. Outgoing centrality
    4. Wave risk (higher values = less risk)
    5. Damage risk (higher values = less risk)
    6. Priority predecessors
    7. GBRMPA zoning criteria
    8. Available space
    9. Location depth

- `wse` : 8-element vector of criteria weights
    1. incoming connectivity
    2. outgoing connectivity
    3. wave
    4. heat
    5. seed predecessors (weights importance of sites highly connected to priority sites for seeding)
    6. seed zones (weights importance of sites highly connected to or within priority zones for seeding)
    7. low cover (weights importance of sites with low cover/high available real estate to plant corals)
    8. depth
"""
function create_seed_matrix(
    A::Matrix{Float64},
    min_area::T,
    wt_in_conn_seed::T,
    wt_out_conn_seed::T,
    wt_waves_seed::T,
    wt_heat_seed::T,
    wt_predec_seed::T,
    wt_predec_zones_seed::T,
    wt_low_cover::T,
    wt_depth_seed::T,
)::Tuple{Matrix{Float64}, Vector{Float64}} where {T<:Float64}
    # Define seeding decision matrix, based on copy of A
    SE = copy(A)

    wse = [
        wt_in_conn_seed,
        wt_out_conn_seed,
        wt_waves_seed,
        wt_heat_seed,
        wt_predec_seed,
        wt_predec_zones_seed,
        wt_low_cover,
        wt_depth_seed,
    ]

    SE[:, 4] = (1 .- SE[:, 4])  # compliment of wave risk
    SE[:, 5] = (1 .- SE[:, 5])  # compliment of heat risk

    # Coral real estate as total area, sites with ≤ min_area to be seeded available filtered out
    # This will also filter out sites with 0 space
    SE[SE[:, 8] .<= min_area, 8] .= NaN

    # Mark "hot" locations for filter
    #SE[vec(A[:, 5] .>= 0.75), 5] .= NaN

    # Filter out identified locations
    SE = SE[vec(.!any(isnan.(SE), dims=2)), :]

    return SE, wse
end


"""
    create_fog_matrix(A, wt_conn_fog , wt_waves_fog, wt_heat_fog, wt_predec_fog, wt_hi_cover)

Create shading specific decision matrix and apply weightings.

# Arguments
- `A` : Criteria  matrix
- `k_area`: Carrying capacity (m²) for coral.
- `wt_conn_fog` : Shading connectivity weight
- `wt_waves_fog` : Wave stress weight
- `wt_heat_fog` : Heat stress weight
- `wt_predec_zones_fog` : Priority zones weight for fogging
- `wt_predec_fog` : Priority predecessor weight for fogging
- `wt_hi_cover` : Weighting for high coral cover when fogging

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
    1. fog connectivity
    2. wave
    3. heat
    4. fog predecessors (weights importance of sites highly connected to priority sites for fogging)
    4. fog zones (weights importance of sites highly connected to or within priority zones)
    5. high cover (weights importance of sites with high cover of coral to fog)
"""
function create_fog_matrix(
    A::Matrix{Float64},
    k_area::Vector{T},
    wt_conn_fog::T,
    wt_waves_fog::T,
    wt_heat_fog::T,
    wt_predec_fog::T,
    wt_predec_zones_fog::T,
    wt_hi_cover,
)::Tuple{Matrix{Float64},Vector{Float64}} where {T<:Float64}

    # Define weights vector
    wsh = [
        wt_conn_fog,
        wt_conn_fog,
        wt_waves_fog,
        wt_heat_fog,
        wt_predec_fog,
        wt_predec_zones_fog,
        wt_hi_cover,
    ]

    # Set up decision matrix to be same size as A
    SH = copy(A)

    # Remove consideration of site depth as shading not accounted for in bleaching model yet
    SH = SH[:, 1:(end - 1)]

    # Get area of coral cover (m²) by subtracting leftover space from total area for coral.
    SH[:, 8] = k_area .- SH[:, 8]

    SH[:, 4] = (1.0 .- SH[:, 4])  # complimentary of wave damage risk

    SH[SH[:, 8].<0, 8] .= 0  # if any negative, scale back to zero
    return SH, wsh
end


"""
    guided_site_selection(d_vars::DMCDA_vars, alg_ind::Int64, log_seed::Bool, log_fog::Bool, pref_seed_locs::AbstractArray{Int64}, pref_fog_locs::AbstractArray{Int64}, rankings_in::Matrix{Int64})

# Arguments
- `d_vars` : DMCDA_vars type struct containing weightings and criteria values for site selection.
- `alg_ind` : integer indicating MCDA aggregation method to use (0: none, 1: order ranking, 2:topsis, 3: vikor)
- `log_seed` : boolean indicating whether seeding sites are being re-assesed at current time
- `log_fog` : boolean indicating whether shading/fogging sites are being re-assesed at current time
- `pref_fog_locs` : previous time step's selection of sites for shading
- `pref_seed_locs` : previous time step's selection of sites for seeding
- `rankings_in` : pre-allocated store for site rankings
- `in_conn` : in-degree centrality
- `out_conn` : out-degree centrality
- `strong_pred` : strongest predecessors

# Returns
Tuple :
    - `pref_seed_locs` : Vector, Indices of preferred seeding locations
    - `pref_fog_locs` : Vector, Indices of preferred shading locations
    - `rankings` : Matrix[n_sites ⋅ 3] where columns are site_id, seeding_rank, shading_rank
        Values of 0 indicate sites that were not considered
"""
function guided_site_selection(
    d_vars::DMCDA_vars,
    alg_ind::T,
    log_seed::B,
    log_fog::B,
    pref_seed_locs::IA,
    pref_fog_locs::IB,
    rankings_in::Matrix{T},
    in_conn::Vector{Float64},
    out_conn::Vector{Float64},
    strong_pred::Vector{Int64};
    methods_mcda=mcda_methods()
)::Tuple{Vector{T}, Vector{T}, Matrix{T}} where {
    T<:Int64,IA<:AbstractArray{<:Int64},IB<:AbstractArray{<:Int64},B<:Bool
}
    use_spatial_group::Int64 = d_vars.use_spatial_group
    site_ids = copy(d_vars.site_ids)
    n_sites::Int64 = length(site_ids)

    # if no sites are available, abort
    if n_sites == 0
        return (
            zeros(Int64, length(pref_seed_sites)),
            zeros(Int64, length(pref_fog_sites)),
            rankings_in,
        )
    end

    n_iv_locs::Int64 = d_vars.n_site_int
    n_spatial_grp::Int64 = d_vars.n_spatial_grp
    priority_sites::Array{Int64} = d_vars.priority_sites[in.(d_vars.priority_sites, [site_ids])]
    priority_zones::Array{String} = d_vars.priority_zones

    zones = d_vars.zones[site_ids]

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

    A, filtered_sites = create_decision_matrix(
        site_ids,
        in_conn,
        out_conn,
        d_vars.leftover_space[site_ids],
        d_vars.dam_prob[site_ids],
        d_vars.heat_stress_prob[site_ids],
        d_vars.site_depth[site_ids],
        predec,
        zones_criteria,
        d_vars.risk_tol,
    )
    if isempty(A)
        # if all rows have nans and A is empty, abort mission
        return (
            zeros(Int64, length(pref_seed_locs)),
            zeros(Int64, length(pref_fog_locs)),
            rankings_in
        )
    end

    # cap to number of sites left after risk filtration
    n_iv_locs = min(n_iv_locs, length(A[:, 1]))

    # if seeding, create seeding specific decision matrix
    if log_seed
        SE, wse = create_seed_matrix(
            A,
            d_vars.min_area,
            d_vars.wt_in_conn_seed,
            d_vars.wt_out_conn_seed,
            d_vars.wt_waves_seed,
            d_vars.wt_heat_seed,
            d_vars.wt_predec_seed,
            d_vars.wt_zones_seed,
            d_vars.wt_lo_cover,
            d_vars.wt_depth_seed,
        )
    end

    # if shading, create fogging specific decision matrix
    if log_fog
        SH, wsh = create_fog_matrix(
            A,
            d_vars.k_area[site_ids][filtered_sites],
            d_vars.wt_conn_fog,
            d_vars.wt_waves_fog,
            d_vars.wt_heat_fog,
            d_vars.wt_predec_fog,
            d_vars.wt_zones_fog,
            d_vars.wt_hi_cover,
        )
    end

    if log_seed && isempty(SE)
        pref_seed_locs = zeros(Int64, n_iv_locs)
    elseif log_seed
            pref_seed_locs, rankings = constrain_reef_group(
                d_vars.reefs,
        pref_seed_locs, s_order_seed = rank_sites!(
            SE,
            wse,
            rankings,
            mcda_func,
            n_iv_locs .* 2,
            2,
        )
                s_order_seed,
                pref_seed_locs,
                rankings,
                n_spatial_grp,
                2,
            )
        end
    end

    if log_fog && isempty(SH)
        pref_fog_locs = zeros(Int64, n_iv_locs)
    elseif log_fog
        pref_fog_locs, s_order_fog = rank_sites!(SH, wsh, rankings, n_iv_locs, mcda_func, 3)
        if use_spatial_group == 1
            pref_fog_locs, rankings = constrain_spatial_group(
                d_vars.reefs,
                s_order_fog,
                pref_fog_locs,
                rankings,
                n_spatial_grp,
                3,
            )
        end
    end

    # Replace with input rankings if seeding or shading rankings have not been filled
    if sum(pref_seed_locs) == 0
        rankings[:, 2] .= rankings_in[:, 2]
    end

    if sum(pref_fog_locs) == 0
        rankings[:, 3] .= rankings_in[:, 3]
    end

    return pref_seed_locs, pref_fog_locs, rankings
end

"""
    constrain_spatial_group(loc_reefs::Vector{String}, s_order::Matrix{Union{Float64,Int64}}, pref_locs::Vector{Float}, 
        rankings::Matrix{Int64},  n_spatial_grp::Int64, rank_col::Int64)
# Arguments
- `loc_reefs` : List of the the reefs/cluster each location sites within
- `s_order` : Ordered set of locations and their aggregate criteria score
- `pref_locs` : Set of location ids to intervene at
- `rankings` : Current ranks of the set of locations
- `n_spatial_grp` : Number of selected locations to allow in the same reef/cluster.
- `rank_col` : Column under which the ranks for the intervention of interest are stored in rankings.

# Returns
Tuple :
    - `pref_locs` : Vector, Indices of preferred intervention locations
    - `rankings` : Matrix[n_sites ⋅ 3] where columns are site_id, seeding_rank, shading_rank
        Values of 0 indicate sites that were not considered
"""
function constrain_spatial_group(
    loc_reefs::Vector{String},
    s_order::Matrix{Union{Float64,Int64}},
    pref_locs::Vector{Int64},
    rankings::Matrix{Int64},
    n_spatial_grp::Int64,
    rank_col::Int64,
)
    # Get full ordering of locations
    loc_order = s_order[:, 1]
    n_loc_iv = length(pref_locs)
    # Get unique reefs/clusters in location set
    unique_reefs = unique(loc_reefs)

    for kk in 1:ceil(Int64, length(loc_order) / 2)
        pref_reefs = loc_reefs[pref_locs] # Reefs/clusters that selected locations sit within

        # Number of times a location appers within each reef/cluster
        sum_reef_pref_locs = [sum(pref_reefs .== rr) for rr in unique_reefs]
        # If more than n_spatial_grp in a reef/cluster, swap out the worst locations
        reef_swap = unique_reefs[findall((sum_reef_pref_locs .> n_spatial_grp))]

        # Locations to replace (lowest ranked sites where too many sites in the same reef/cluster)
        replace_locs = [
            pref_locs[pref_reefs .== reef][(n_spatial_grp + 1):end] for reef in reef_swap
        ]

        if !isempty(replace_locs)
            replace_locs = replace_locs[1:end...]

            pref_locs = setdiff(pref_locs, replace_locs) # Remove locations from preferred sites
            loc_order = setdiff(loc_order, replace_locs) # Remove locations from site order
            add_locs = setdiff(loc_order, pref_locs) # Get next highly ranked to replace removed locations with

            # New preferred location set
            pref_locs = vec([pref_locs... add_locs[1:(n_loc_iv - length(pref_locs))]...])
        else
            break
        end
    end

    # Addremoved sites at end of location order
    removed_sites = setdiff(s_order[:, 1], loc_order)
    s_order[:, 1] .= [
        loc_order[1:n_loc_iv]..., removed_sites..., loc_order[(n_loc_iv + 1):end]...
    ]

    # Match by site_id and assign rankings to log
    align_rankings!(rankings, s_order, rank_col)
    return pref_locs, rankings
end

"""
    unguided_site_selection(pref_seed_sites, pref_fog_sites, seed_years, fog_years, n_site_int, available_space, depth)

Randomly select seed/fog site locations for the given year, constraining to sites with max. carrying capacity > 0.

# Arguments
- `pref_seed_locs` : Previously selected seeding locations
- `pref_fog_locs` : Previously selected fogging locations
- `seed_years` : bool, indicating whether to seed this year or not
- `fog_years` : bool, indicating whether to fog this year or not
- `n_site_int` : int, number of sites to intervene on
- `available_space` : vector/matrix : space available at each site (`k` value)
- `depth` : vector of site ids found to be within desired depth range

# Returns
Tuple, of vectors indicating preferred seeding and shading locations by location index
"""
function unguided_site_selection(
    pref_seed_locs,
    pref_fog_locs,
    seed_years,
    fog_years,
    n_site_int,
    available_space,
    depth
)::Tuple{Vector, Vector}
    # Unguided deployment, seed/fog corals anywhere so long as available_space > 0.0
    # Only sites that have available space are considered, otherwise a zero-division error may occur later on.

    # Select sites (without replacement to avoid duplicate sites)
    candidate_sites = depth[(available_space.>0.0)[depth]]  # Filter down to site ids to be considered
    num_sites = length(candidate_sites)
    s_n_site_int = num_sites < n_site_int ? num_sites : n_site_int

    if seed_years
        pref_seed_locs = zeros(Int64, n_site_int)
        pref_seed_locs[1:s_n_site_int] .= StatsBase.sample(candidate_sites, s_n_site_int; replace=false)
    end

    if fog_years
        pref_fog_locs = zeros(Int64, n_site_int)
        pref_fog_locs[1:s_n_site_int] .= StatsBase.sample(
            candidate_sites, s_n_site_int; replace=false
        )
    end

    return pref_seed_locs[pref_seed_locs .> 0], pref_fog_locs[pref_fog_locs .> 0]
end


"""
    summary_stat_env(env_layer::NamedDimsArray dims::Union{Symbol,Tuple{Symbol,Symbol}}; w=0.5)::Vector{Float64}

Calculates mean over specified dimensions plus half the standard deviation.

# Arguments
- `env_layer` : Environmental data layer to calculate the mean of.
- `dims` : Dimensions to aggregate over.
- `w` : Weighting for std offset to mean.

# Returns
Weighted combination of mean and standard deviation of the projected environmental
conditions (e.g., DHWs, wave stress, etc):
    (μ * w) + (σ * (1 - w))
"""
function summary_stat_env(
    env_layer::AbstractArray,
    dims::Union{Int64,Symbol,Tuple{Symbol,Symbol}};
    w=0.5,
)::Vector{Float64}
    return vec((mean(env_layer, dims=dims).* w) .+ (std(env_layer, dims=dims) .* (1.0 - w)))
end

"""
    within_depth_bounds(loc_depth::Vector{T}, depth_max::T, depth_min::T)::BitVector{T} where {T<:Float64}

Determines whether a location is within the min/max depth bounds.
Used to filter locations based on their depth for location selection.

# Arguments
- `loc_depth` : Depths of considered locations (typically the median depth)
- `depth_max` : Maximum depth for each considered location
- `depth_min` : Minimum depth for each considered location

# Returns
BitVector, of logical indices indicating locations which satisfy the depth criteria.
"""
function within_depth_bounds(
    loc_depth::Vector{T}, depth_max::T, depth_min::T
)::BitVector where {T<:Float64}
    return (loc_depth .<= depth_max) .& (loc_depth .>= depth_min)
end

end
