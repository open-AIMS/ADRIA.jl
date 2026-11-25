"""Coral growth functions"""

using Distributions
using SpecialFunctions

using SparseArrays

"""
    proportional_adjustment!(coral_cover::Matrix{T}, loc_cover_cache::Vector{T})::Nothing where {T<:Float64}

Helper method to proportionally adjust coral cover, such that:
- `coral_cover` ∈ [0, 1].
- `sum(coral_cover, dims=1)` ≯  1.0

Modifies arrays in-place.

# Arguments
- `coral_cover` : Coral cover ∈ [0, 1]
- `loc_cover_cache` : temporary cache

# Returns
nothing
"""
function proportional_adjustment!(
    coral_cover::Union{SubArray{T},Matrix{T}}, loc_cover_cache::Vector{T}
)::Nothing where {T<:Float64}
    loc_cover_cache .= vec(sum(coral_cover; dims=1))
    loc_cover_cache[loc_cover_cache .≈ 1.0] .= 1.0
    if any(loc_cover_cache .> 1.0)
        exceeded::BitVector = vec(loc_cover_cache .> 1.0)
        msg = """
            Cover exceeded bounds, constraining to be within available space
            This indicates an issue with the model.
            Cover - Max Cover: $(sum(loc_cover_cache[exceeded] .- 1.0))
        """
        @debug msg

        # @views @. coral_cover[:, exceeded] = (
        #     coral_cover[:, exceeded] / loc_cover_cache[exceeded]'
        # )
        for loc in findall(exceeded)
            scale_factor = 1.0 / loc_cover_cache[loc]
            for sp_sc in axes(coral_cover, 1)
                coral_cover[sp_sc, loc] *= scale_factor
            end
        end

        coral_cover .= max.(coral_cover, 0.0)
    end

    return nothing
end
function proportional_adjustment!(
    coral_cover::Union{SubArray{T,3},Array{T,3}}, loc_cover_cache::Vector{T}
)::Nothing where {T<:Float64}
    loc_cover_cache .= vec(sum(coral_cover; dims=(1, 2)))
    loc_cover_cache[loc_cover_cache .≈ 1.0] .= 1.0
    if any(loc_cover_cache .> 1.0)
        exceeded::Vector{Int64} = findall(vec(loc_cover_cache .> 1.0))
        msg = """
            Cover exceeded bounds, constraining to be within available space
            This indicates an issue with the model.
            Cover - Max Cover: $(sum(loc_cover_cache[exceeded] .- 1.0))
        """
        @debug msg

        for idx in exceeded
            @views @. coral_cover[:, :, idx] = (
                coral_cover[:, :, idx] / loc_cover_cache[idx]'
            )
        end

        coral_cover .= max.(coral_cover, 0.0)
    end

    return nothing
end
"""
    proportional_adjustment!(coral_cover::Union{SubArray{T},Matrix{T}})::Nothing where {T<:Float64}

Adjust relative coral cover based on the proportion each size class contributes to area
covered. Assumes 1.0 represents 100% of available location area.
"""
function proportional_adjustment!(
    coral_cover::Union{SubArray{T},Matrix{T}}
)::Nothing where {T<:Float64}
    loc_cover_cache = zeros(size(coral_cover, 2))

    return proportional_adjustment!(coral_cover, loc_cover_cache)
end

"""
    fecundity_scope!(
        fec_groups::AbstractMatrix{T},
        fecundity_per_m²::AbstractMatrix{T},
        C_cover_t::AbstractArray{T,3},
        loc_area::AbstractMatrix{T};
        apply_allee_effect::Bool=true
    )::Nothing where {T<:Float64}

The scope that different coral groups and size classes have for
producing larvae without consideration of environment.

Coral fecundity per coral area of the different size classes.
When multiplied by the relative cover of each size class within taxa,
this produces an estimate of the relative fecundity of each coral group and size.
Total relative fecundity of a group is then calculated as the sum of
fecundities across size classes.

# Arguments
- `fec_groups` : Matrix[n_groups, n_locs], memory cache to place results into
- `fecundity_per_m²` : Matrix[n_groups, n_sizes], coral fecundity parameters (in per m²) for each species/size class
- `C_cover_t` : Matrix[n_groups, n_sizes, n_locs], of coral cover values for the previous time step
- `loc_area` : Vector[n_locs], total location area in m²
- `apply_allee_effect` : Whether to suppress fecundity at low density via `allee_effect()`
"""
function fecundity_scope!(
    fec_groups::AbstractMatrix{T},
    fecundity_per_m²::AbstractMatrix{T},
    C_cover_t::AbstractArray{T,3},
    loc_area::AbstractMatrix{T};
    apply_allee_effect::Bool=true
)::Nothing where {T<:Float64}
    ae_coef = 0.0
    for loc in axes(C_cover_t, 3)
        ae_coef = apply_allee_effect ? allee_effect(sum(C_cover_t[:, :, loc])) : 1.0
        for grp in axes(C_cover_t, 1)
            for sz in axes(C_cover_t, 2)
                @views fec_groups[grp, loc] +=
                    fecundity_per_m²[grp, sz] *
                    C_cover_t[grp, sz, loc] *
                    loc_area[1, loc]
            end
            @views fec_groups[grp, loc] *= ae_coef
        end
    end

    return nothing
end

"""
# Arguments
- `C_cover` : Location relative coral cover
"""
function allee_effect(relative_cover::Float64)
    # Below 0.0001 this will quickly start to yield negative values
    min_rel_cover = 0.0001
    x = relative_cover < min_rel_cover ? log10(min_rel_cover) : log10(relative_cover)
    return 0.999 + 0.0169 * x + 0.00197 * x^2 + 0.0125 * x^3
end

"""
    settler_density(α, β, L)

Density potential for settlers estimated with a Beverton-Holt (B-H) function, following [1]
and as detailed in [2].

Note for β: "For corals, the actual number of 6-month old recruits for each coral group
    is generated [...] following a Poisson distribution with recruitment event rate λ.

# Examples
```julia
settler_density(2.5, 5000.0, L)
```

# Arguments
- `α` : Maximum achievable density (settlers/m²) for a 100% free space (set to 2.5 in [1] for Corymbose)
- `β` : Stock of larvae required to produce one-half the maximum settlement (larvae/m²),
        i.e., α/2(m²), set to 5000 in [1].
- `L` : Available larval pool

# Returns
Settler density (settlers / m²)

# References
1. Bozec, Y.-M., Hock, K., Mason, R. A. B., Baird, M. E.,
     Castro-Sanguino, C., Condie, S. A., Puotinen, M.,
     Thompson, A., & Mumby, P. J. (2022).
   Cumulative impacts across Australia's Great Barrier Reef:
     A mechanistic evaluation.
   Ecological Monographs, 92(1), e01494.
   https://doi.org/10.1002/ecm.1494

2. Haddon, M. (2011). Modelling and quantitative methods in fisheries. CRC Press/Chapman and
    Hall, Boca Raton, Florida, USA.
"""
function settler_density(α::T, β::T, L::T)::Float64 where {T<:Float64}
    return (α .* L) ./ (β .+ L)
end

"""
    recruitment_rate(larval_pool::AbstractArray{T,2}, A::AbstractArray{T}; α=2.5, β=5000.0)

Calculates coral recruitment for each species/group and location.

# Arguments
- `larval_pool` : Available larval pool
- `A` : Available space, i.e., the substratum that is suitable for coral recruitment, in m²
- `α` : Maximum achievable density (settlers/m²) for a 100% free space
- `β` : Stock of larvae required to produce 50% of the maximum settlement

# Returns
λ, total coral recruitment for each coral taxa and location based on a Poisson distribution.
"""
function recruitment_rate(
    larval_pool::AbstractArray{T,2}, A::AbstractArray{T};
    α::Union{T,Vector{T}}=2.5, β::Union{T,Vector{T}}=5000.0, rng=Random.GLOBAL_RNG
)::Matrix{T} where {T<:Float64}
    sd = settler_density.(α, β, larval_pool) .* A'

    @views sd[sd .> 0.0] .= rand.(rng, Poisson.(sd[sd .> 0.0]))

    return sd
end

"""
    settler_cover(fec_scope::T, conn::SparseMatrixCSC{Float64,Int64}, leftover_space::V, α::V, β::V, basal_area_per_settler::V, potential_settlers::T)::T where {T<:Matrix{Float64},V<:Vector{Float64}}

Determine area settled by recruited larvae.

Note: Units for all areas are assumed to be in m².

# Arguments
- `fec_scope` : Fecundity scope
- `conn` : Connectivity between locations (rows: source locations; cols: sink locations)
- `leftover_space` : Difference between locations' maximum carrying capacity and current
    coral cover (in m²)
- `α` : max number of settlers / m²
- `β` : larvae / m² required to produce 50% of maximum settlement
- `basal_area_per_settler` : area taken up by a single settler
- `potential_settlers` : matrix to (re)use as cache to avoid memory allocations

# Returns
Area covered by recruited larvae (in m²)
"""
function settler_cover(
    fec_scope::T,
    conn::SparseMatrixCSC{Float64,Int64},
    leftover_space::V,
    α::V,
    β::V,
    basal_area_per_settler::V,
    potential_settlers::T;
    rng::AbstractRNG=Random.GLOBAL_RNG
)::T where {T<:AbstractMatrix{Float64},V<:Vector{Float64}}

    # Send larvae out into the world (reuse potential_settlers to reduce allocations)
    # Note, conn rows need not sum to 1.0 as this missing probability accounts for larvae
    # which do not settle. Pers comm with C. Ani (2023-01-29 13:24 AEST).
    # [Larval pool for each location in larvae/m²] * [survival rate]
    # this is known as in-water mortality.
    # Set to 0.0 as it is now taken care of by connectivity data.
    # Mwater::Float64 = 0.0

    # Efficient Dense * Sparse multiplication
    mul!(potential_settlers, fec_scope, conn)

    # Larvae have landed, work out how many are recruited
    # Determine area covered by recruited larvae (settler cover) per m^2
    # recruits per m^2 per site multiplied by area per settler
    return recruitment_rate(potential_settlers, leftover_space; α=α, β=β, rng=rng) .*
           basal_area_per_settler
end

"""
    constrain_recruitment(recruitment, agg_cover_above_threshold_mask, C_cover_t, habitable_loc_areas; round_threshold=0.95)

To prevent overgrowth when adding recruitment to C_cover_t, set a threshold above which
recruits are either set to zero or rescaled.
"""
function constrain_recruitment!(
    recruitment::Matrix{Float64},
    agg_cover_above_threshold_mask::BitVector,
    C_cover_t::Array{Float64,3},
    habitable_loc_areas::Vector{Float64};
    round_threshold::Float64=0.95
)
    # Relative coral cover and recruits cover for each location
    C_rel_cover = loc_coral_cover(C_cover_t) ./ habitable_loc_areas
    loc_recruits_rel_cover = loc_recruits_cover(recruitment) ./ habitable_loc_areas

    # Set recruitment to 0 where C_cover_t is already above round_threshold
    cover_above_threshold_mask = C_rel_cover .> round_threshold
    recruitment[:, cover_above_threshold_mask] .= 0.0

    # Mask for locations with (C_cover + recruits) > threshold and C_cover <= threshold
    cover_with_recruits = C_rel_cover .+ loc_recruits_rel_cover
    rec_above_threshold = cover_with_recruits .> round_threshold
    agg_cover_above_threshold_mask .= rec_above_threshold .&& .!cover_above_threshold_mask

    if any(agg_cover_above_threshold_mask)
        is_test_env() ? nothing : @debug "Constraining recruits within error bounds."

        exceeded_agg_cover_locs = habitable_loc_areas[agg_cover_above_threshold_mask]
        threshold_cap = round_threshold .* exceeded_agg_cover_locs
        exceeded_locs = loc_coral_cover(C_cover_t[:, :, agg_cover_above_threshold_mask])
        recruit_cover = loc_recruits_cover(recruitment[:, agg_cover_above_threshold_mask])

        recruits_scale_factor = (threshold_cap .- exceeded_locs) ./ recruit_cover

        @assert !any(recruits_scale_factor .<= 0)

        # Rescale recruits in target locations to fit within error bounds
        n_groups = size(C_cover_t, 1)
        recruitment[:, agg_cover_above_threshold_mask] .*= repeat(
            recruits_scale_factor', n_groups
        )
    end
end
