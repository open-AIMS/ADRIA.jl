"""Coral growth functions"""

using FLoops
using Distributions
using SpecialFunctions

"""
    growth_rate(linear_extension::Matrix{Float64}, diam_bin_widths::Matrix{Float64})::Matrix{Float64}

Determine the rate of growth representing the proportion of each size class that moves
up a size class each (yearly) time step. Values > 1 indicate transitions to higher size
classes occurs more than once per time step.

# Arguments
- `linear_extension` : Linear extension in cm/year
- `diam_bin_widths` : diameter of each size class (bin) in cm

# Returns
Matrix, of size \$[n_{species} ⋅ n_{classes}]\$ indicating proportional growth rates
for each.
"""
function growth_rate(
    linear_extension::Matrix{Float64}, diam_bin_widths::Matrix{Float64}
)::Matrix{Float64}
    return ((2.0 .* linear_extension) ./ diam_bin_widths)'
end

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

        @views @. coral_cover[:, exceeded] = (
            coral_cover[:, exceeded] / loc_cover_cache[exceeded]'
        )

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

        @floop for idx in exceeded
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
    depth_coefficient(d::Union{Int64,Float64})::Float64

Model by Baird et al., [1] providing an indication of a relationship between bleaching
and depth.

# Arguments
- `d` : median depth of location in meters

# Returns
Proportion of population affected by bleaching at depth `x`.
Values are constrained such that 0.0 <= x <= 1.0

# References
1. Baird, A., Madin, J., Álvarez-Noriega, M., Fontoura, L., Kerry, J., Kuo, C.,
     Precoda, K., Torres-Pulliza, D., Woods, R., Zawada, K., & Hughes, T. (2018).
   A decline in bleaching suggests that depth can provide a refuge from global
     warming in most coral taxa.
   Marine Ecology Progress Series, 603, 257-264.
   https://doi.org/10.3354/meps12732
"""
function depth_coefficient(d::Union{Int64,Float64})::Float64
    return max(min(exp(-0.07551 * (d - 2.0)), 1.0), 0.0)
end

"""
    bleaching_mortality!(Y::AbstractArray{Float64,2}, capped_dhw::AbstractArray{Float64,2},
        depth_coeff::AbstractArray{Float64}, tstep::Int64, depth::Vector{Float64},
        s::Vector{Float64}, dhw::AbstractArray{Float64}, a_adapt::Vector{Float64},
        n_adapt::Real)::Nothing

Calculates and applies bleaching mortality, taking into account depth and bleaching
sensitivity of corals. Model is adapted from Bozec et al., [2], itself based on data
from Hughes et al., [3] (bleaching sensitivity) and Baird et al., [1] (relationship
between bleaching and depth).

# Arguments
- `Y` : Matrix to save results into, of shape \$n_{species} ⋅ n_{locations}\$
- `tstep` : Current time step
- `depth` : Mean site depth (m) for each site
- `s` : Bleaching sensitivity of corals (relative values) for each taxa/size class
- `dhw` : Degree Heating Week experienced at site
- `a_adapt` : Level of assisted adaptation (DHW reduction)
- `n_adapt` : Level of natural adaptation (DHW reduction linearly scaled over time)

# Returns
Nothing

# References
1. Baird, A., Madin, J., Álvarez-Noriega, M., Fontoura, L., Kerry, J., Kuo, C.,
     Precoda, K., Torres-Pulliza, D., Woods, R., Zawada, K., & Hughes, T. (2018).
   A decline in bleaching suggests that depth can provide a refuge from global
     warming in most coral taxa.
   Marine Ecology Progress Series, 603, 257-264.
   https://doi.org/10.3354/meps12732

2. Bozec, Y.-M., Hock, K., Mason, R. A. B., Baird, M. E., Castro-Sanguino, C.,
     Condie, S. A., Puotinen, M., Thompson, A., & Mumby, P. J. (2022).
   Cumulative impacts across Australia's Great Barrier Reef: A mechanistic evaluation.
   Ecological Monographs, 92(1), e01494.
   https://doi.org/10.1002/ecm.1494

3. Hughes, T. P., Kerry, J. T., Baird, A. H., Connolly, S. R., Dietzel, A., Eakin, C. M.,
     Heron, S. F., Hoey, A. S., Hoogenboom, M. O., Liu, G., McWilliam, M. J., Pears, R. J.,
     Pratchett, M. S., Skirving, W. J., Stella, J. S., & Torda, G. (2018).
   Global warming transforms coral reef assemblages.
   Nature, 556(7702), 492-496.
   https://doi.org/10.1038/s41586-018-0041-2
"""
function bleaching_mortality!(Y::AbstractArray{Float64,2},
    capped_dhw::AbstractArray{Float64,2}, depth_coeff::AbstractArray{Float64},
    tstep::Int64, depth::Vector{Float64}, s::Vector{Float64},
    dhw::AbstractArray{Float64}, a_adapt::Vector{Float64}, n_adapt::Real)::Nothing

    # Initial mortality
    # Bozec et al., (2022) derive their model from Hughes et al., (2018, Fig 2C).
    # They digitized data from the figure and fit a linear model to obtain
    # intercept of 0.168 and slope of 0.347 (rounded below).
    # In exponential form, this is: `InitMort(%) = exp(0.168+0.347*DHW) - 1.`
    # Bozec et al., state this is valid for DHW values 0 - 10, but also
    # assume it remains valid for values > 10.0
    # A depth coefficient is also introduced to account for the reduced
    # experienced heat at greater depths.

    # The model is modified to incorporate adaptation effect but maximum
    # reduction is to capped to 0.
    @. capped_dhw =
        min.(ℯ^(0.17 + 0.35 * max(0.0, dhw' - (a_adapt + (tstep * n_adapt)))) - 1.0, 100.0)
    @. depth_coeff = ℯ^(-0.07551 * (depth - 2.0))

    # Estimate long-term bleaching mortality with an estimated depth coefficient and
    # initial bleaching mortality (models from Bozec et al., 2022)
    #
    # BleachMort(%) = 100*(1-(1-capped_dhw/100)^6)
    #
    # As we want values between 0 - 1 rather than %, we drop the `100 *`
    # We also want remaining population, so we also drop the initial `1 - `
    # End result is how much coral survives a bleaching event.
    @. Y = (1.0 - ((depth_coeff' * s) * (capped_dhw / 100.0)))^6.0

    return nothing
end

"""
    bleaching_mortality!(cover::AbstractArray{Float64, 3}, dhw::Vector{Float64},
        depth_coeff::Vector{Float64}, stdev::AbstractMatrix{Float64},
        dist_t_1::AbstractArray{Float64, 3}, dist_t::AbstractArray{Float64, 3},
        prop_mort::SubArray{Float64})::Nothing

Applies bleaching mortality by assuming critical DHW thresholds are normally distributed for
all non-juvenile (> 5cm diameter) size classes.

Distributions are informed by learnings from Bairos-Novak et al., [1] and (unpublished)
data referred to in Hughes et al., [2]. Juvenile mortality is assumed to be primarily
represented by other factors (i.e., background mortality; see Álvarez-Noriega et al., [3]).
The proportion of the population which bleached is estimated with the Cumulative Density
Function. Bleaching mortality is then estimated with a depth-adjusted coefficient
(from Baird et al., [4]).

# Arguments
- `cover` : Coral cover for current timestep
- `dhw` : DHW for all represented locations
- `depth_coeff` : Pre-calculated depth coefficient for all locations
- `stdev` : Standard deviation of DHW tolerance
- `dist_t_1` : Critical DHW threshold distribution for current timestep, for all species and
    locations
- `dist_t` : Critical DHW threshold distribution for next timestep, for all species and
    locations
- `prop_mort` : Cache to store records of bleaching mortality (stores mortalities for
    t-1 and t)
- `n_sizes` : Number of size classes within each functional group

# References
1. Bairos-Novak, K.R., Hoogenboom, M.O., van Oppen, M.J.H., Connolly, S.R., 2021.
     Coral adaptation to climate change: Meta-analysis reveals high heritability across
     multiple traits.
   Global Change Biology 27, 5694-5710.
   https://doi.org/10.1111/gcb.15829

2. Hughes, T.P., Kerry, J.T., Baird, A.H., Connolly, S.R., Dietzel, A., Eakin, C.M.,
     Heron, S.F., Hoey, A.S., Hoogenboom, M.O., Liu, G., McWilliam, M.J., Pears, R.J.,
     Pratchett, M.S., Skirving, W.J., Stella, J.S., Torda, G., 2018.
   Global warming transforms coral reef assemblages.
   Nature 556, 492-496.
   https://doi.org/10.1038/s41586-018-0041-2

3. Álvarez-Noriega, M., Baird, A.H., Bridge, T.C.L., Dornelas, M., Fontoura, L.,
     Pizarro, O., Precoda, K., Torres-Pulliza, D., Woods, R.M., Zawada, K.,
     Madin, J.S., 2018.
   Contrasting patterns of changes in abundance following a bleaching event between
     juvenile and adult scleractinian corals.
   Coral Reefs 37, 527-532.
   https://doi.org/10.1007/s00338-018-1677-y

4. Baird, A., Madin, J., Álvarez-Noriega, M., Fontoura, L., Kerry, J., Kuo, C.,
     Precoda, K., Torres-Pulliza, D., Woods, R., Zawada, K., & Hughes, T. (2018).
   A decline in bleaching suggests that depth can provide a refuge from global
     warming in most coral taxa.
   Marine Ecology Progress Series, 603, 257-264.
   https://doi.org/10.3354/meps12732
"""
function bleaching_mortality!(
    cover::AbstractArray{Float64,3},
    dhw::Vector{Float64},
    depth_coeff::Vector{Float64},
    stdev::AbstractMatrix{Float64},
    dist_t_1::AbstractArray{Float64,3},
    dist_t::AbstractArray{Float64,3},
    prop_mort::SubArray{Float64}
)::Nothing
    n_groups, n_sizes, n_locs = size(cover)

    non_juveniles = 2:n_sizes

    # Adjust distributions for each functional group over all locations, ignoring juveniles
    # we assume the high background mortality of juveniles includes DHW mortality
    for loc in 1:n_locs
        # Skip locations that have no heat stress
        if dhw[loc] == 0.0
            continue
        end

        # Determine bleaching mortality for each non-juvenile species/size class
        for grp in 1:n_groups
            for sc in non_juveniles
                # Skip location if there is no population
                if cover[grp, sc, loc] == 0.0
                    continue
                end

                μ::Float64 = dist_t_1[grp, sc, loc]
                affected_pop::Float64 = truncated_normal_cdf(
                    # Use previous mortality threshold as minimum
                    dhw[loc], μ, stdev[grp, sc], prop_mort[1, grp, sc, loc], μ + HEAT_UB
                )

                mort_pop::Float64 = 0.0
                if affected_pop > 0.0
                    # Calculate depth-adjusted bleaching mortality
                    mort_pop = (affected_pop * depth_coeff[loc])

                    # Set values close to 0.0 (e.g., 1e-214) to 0.0
                    # https://github.com/JuliaLang/julia/issues/23376#issuecomment-324649815
                    if (mort_pop + one(mort_pop)) ≈ one(mort_pop)
                        mort_pop = 0.0
                    end
                end

                prop_mort[2, grp, sc, loc] = mort_pop
                if mort_pop > 0.0
                    # Re-create distribution
                    # Use same stdev as target size class to maintain genetic variance
                    # pers comm K.B-N (2023-08-09 16:24 AEST)
                    dist_t[grp, sc, loc] = truncated_normal_mean(
                        μ, stdev[grp, sc], mort_pop, μ + HEAT_UB
                    )

                    # Update population
                    cover[grp, sc, loc] = cover[grp, sc, loc] * (1.0 - mort_pop)
                end
            end
        end
    end

    return nothing
end

"""
    breeders(μ_o::T, μ_s::T, h²::T)::T where {T<:Float64}

Apply Breeder's equation.

```
S = μ_s - μ_o
μ_t1 = μ_o + (S * h²)
```

# Arguments
- `μ_o` : Mean of original population
- `μ_s` : Mean of next generation
- `h²` : Narrow-sense heritability

# Returns
Updated distribution mean
"""
function breeders(μ_o::T, μ_s::T, h²::T)::T where {T<:Float64}
    return μ_o + ((μ_s - μ_o) * h²)
end

"""
    _shift_distributions!(cover::SubArray{F}, growth_rate::SubArray{F}, dist_t::SubArray{F})::Nothing where {F<:Float64}

Combines distributions between size classes > 1 to represent the shifts that occur as each
size class grows. Weights for the distributions are based on proportional cover and the
assumed growth rates for each size class.

i.e., \$(w_{i+1,1}, w_{i+1,2}) := (c_{i-1,i} / sum(c_{i-1,i})) * (g_{i-1,i} / sum(g_{i-1,i}))\$

where \$w\$ are the weights and \$g\$ is the growth rates.

# Arguments
- `cover` : Coral cover for \$t-1\$
- `growth_rate` : Growth rates for the given size classes/species
- `dist_t` : Critical DHW threshold distribution for timestep \$t\$
- `stdev` : Standard deviations of coral DHW tolerance
"""
function _shift_distributions!(
    cover::SubArray{F},
    growth_rate::SubArray{F},
    dist_t::SubArray{F}
)::Nothing where {F<:Float64}
    # Weight distributions based on growth rate and cover
    # Do from largest size class to size class 2
    # (values for size class 1 gets replaced by recruitment process)
    prop_growth = MVector{2,F}(0.0, 0.0)
    for i in length(growth_rate):-1:2
        # Skip size class if nothing is moving up
        sum(view(cover, (i - 1):i)) == 0.0 ? continue : false

        prop_growth .= @views (cover[(i - 1):i] ./ sum(cover[(i - 1):i])) .*
            (growth_rate[(i - 1):i] ./ sum(growth_rate[(i - 1):i]))
        if sum(prop_growth) == 0.0
            continue
        end

        # Weighted sum
        dist_t[i] =
            (dist_t[i - 1] * prop_growth[1] + dist_t[i] * prop_growth[2]) / sum(prop_growth)
    end

    return nothing
end

"""
    adjust_DHW_distribution!(
        cover::SubArray{F},
        n_sizes::Int64,
        dist_t::Matrix{F},
        growth_rate::Matrix{F}
    )::Nothing where {F<:Float64}

Adjust critical DHW thresholds for a given group/size class distribution as mortalities
affect the distribution over time, and corals mature (moving up size classes).

# Arguments
- `cover` : Coral cover (for timestep \$t-1\$ and \$t\$)
- `n_sizes` : Number of coral size classes represented
- `dist_t` : Distributions for timestep \$t\$
- `growth_rate` : Growth rates for each group/size class
"""
function adjust_DHW_distribution!(
    cover::SubArray{F},
    n_sizes::Int64,
    dist_t::Matrix{F},
    growth_rate::Matrix{F}
)::Nothing where {F<:Float64}
    _, n_sp_sc, n_locs = size(cover)

    step::Int64 = n_sizes - 1

    # Adjust population distribution
    for (grp, sc1) in enumerate(1:n_sizes:n_sp_sc)
        for loc in 1:n_locs
            sc_end::Int64 = sc1 + step

            # Skip if no cover
            if sum(@view(cover[1, sc1:sc_end, loc])) == 0.0
                continue
            end

            # Combine distributions using a MixtureModel for all non-juvenile size
            # classes (we pass in all relevant size classes for the functional group here).
            @views _shift_distributions!(
                cover[1, sc1:sc_end, loc], growth_rate[grp, :], dist_t[sc1:sc_end, loc]
            )
        end
    end

    return nothing
end
function adjust_DHW_distribution!(
    cover_t_1::SubArray{T,3},
    dist_t::AbstractArray{T,3},
    growth_rate::Matrix{T}
)::Nothing where {T<:Float64}
    groups, _, locs = axes(cover_t_1)

    for grp in groups
        for loc in locs
            if sum(@view(cover_t_1[grp, :, loc])) == 0.0
                continue
            end

            @views _shift_distributions!(
                cover_t_1[grp, :, loc],
                growth_rate[grp, :],
                dist_t[grp, :, loc]
            )
        end
    end

    return nothing
end

function settler_DHW_tolerance!(
    c_mean_t_1::Array{F,3},
    c_mean_t::Array{F,3},
    k_area::Vector{F},
    tp::YAXArray{F,2,Matrix{F}},
    settlers::Matrix{F},
    fec_params_per_m²::Matrix{F},
    h²::F
)::Nothing where {F<:Float64}
    groups, sizes, _ = axes(c_mean_t_1)

    sink_loc_ids::Vector{Int64} = findall(k_area .> 0.0)
    source_locs::BitVector = BitVector(undef, length(k_area))  # cache for source locations

    # Number of reproductive size classes for each group
    n_reproductive::Matrix{Int64} = sum(fec_params_per_m² .> 0.0; dims=2)

    # Cache to hold indication of which size classes are considered reproductive
    reproductive_sc::BitVector = falses(sizes)

    @inbounds for sink_loc in sink_loc_ids
        sink_settlers = @view(settlers[:, sink_loc])'
        if sum(sink_settlers) .== 0.0
            # Only update locations where recruitment occurred
            continue
        end

        # Find sources for each sink
        # Note: source locations can include the sink location (self-seeding)
        @inbounds source_locs .= @view(tp.data[:, sink_loc]) .> 0.0

        # Calculate contribution to cover to determine weights for each functional group
        @inbounds w::Matrix{F} = sink_settlers .* @view(tp.data[source_locs, sink_loc])
        w_per_group::Matrix{F} = w ./ sum(w; dims=1)

        # If there is any influence from another location for a group, the tolerance
        # values should be updated.
        update_group::BitMatrix = sum(w_per_group; dims=1) .> 0.0

        @inbounds for grp in groups
            # Determine weights based on contribution to recruitment.
            # This weights the recruited corals by the size classes and source locations
            # which contributed to recruitment.
            if update_group[1, grp]
                # Get distribution mean of reproductive size classes at source locations
                # recalling that source locations may include the sink location due to
                # self-seeding.
                reproductive_sc .= @view(fec_params_per_m²[grp, :]) .> 0.0

                # Determine combined mean
                # https://en.wikipedia.org/wiki/Mixture_distribution#Properties
                settler_means::SubArray{F} = @view(
                    c_mean_t_1[grp, reproductive_sc, source_locs]
                )

                recruit_μ::F = 0.0
                @inbounds for i in axes(settler_means, 1), j in axes(settler_means, 2)
                    recruit_μ += settler_means[i, j] * w_per_group[j, grp]
                end
                recruit_μ /= n_reproductive[grp]

                # Mean for generation t is determined through Breeder's equation
                @views c_mean_t[grp, 1, sink_loc] = breeders(
                    c_mean_t_1[grp, 1, sink_loc], recruit_μ, h²
                )
            end
        end
    end

    return nothing
end

"""
    fecundity_scope!(fec_groups::Array{Float64, 2}, fec_all::Array{Float64, 2},
                     fec_params::Array{Float64}, C_cover_t::Array{Float64, 2},
                     k_area::Array{Float64})::Nothing

    fecundity_scope!(fec_groups::AbstractMatrix{T}, fec_all::AbstractArray{T, 3},
                     fec_params::AbstractMatrix{T}, C_cover_t::AbstractArray{T, 3},
                     loc_area::AbstractMatrix{T})::Nothing where {T<:Float64}

The scope that different coral groups and size classes have for
producing larvae without consideration of environment.

Coral fecundity per coral area of the different size classes.
When multiplied by the relative cover of each size class within taxa,
this produces an estimate of the relative fecundity of each coral group and size.
Total relative fecundity of a group is then calculated as the sum of
fecundities across size classes.

# Arguments
- `fec_groups` : Matrix[n_classes, n_locs], memory cache to place results into
- `fec_all` : Matrix[n_taxa, n_locs], temporary cache to place intermediate fecundity values into
- `fec_params` : Vector, coral fecundity parameters (in per m²) for each species/size class
- `C_cover_t` : Matrix[n_taxa, n_locs], of coral cover values for the previous time step
- `loc_area` : Vector[n_locs], total location area in m²
"""
function fecundity_scope!(
    fec_groups::AbstractMatrix{T},
    fec_all::AbstractMatrix{T},
    fec_params::AbstractVector{T},
    C_cover_t::AbstractMatrix{T},
    loc_area::AbstractMatrix{T}
)::Nothing where {T<:Float64}
    n_groups::Int64 = size(fec_groups, 1)   # number of coral groups: 5
    n_group_and_size::Int64 = size(fec_params, 1)  # number of coral size classes: 35
    n_classes::Int64 = Int64(n_group_and_size / n_groups)

    fec_all .= fec_params .* C_cover_t .* loc_area
    for (i, (s, e)) in enumerate(
        zip(1:n_classes:n_group_and_size, n_classes:n_classes:(n_group_and_size + 1))
    )
        @views fec_groups[i, :] .= vec(sum(fec_all[s:e, :]; dims=1))
    end

    return nothing
end
function fecundity_scope!(
    fec_groups::AbstractMatrix{T},
    fec_all::AbstractArray{T,3},
    fec_params::AbstractMatrix{T},
    C_cover_t::AbstractArray{T,3},
    loc_area::AbstractMatrix{T}
)::Nothing where {T<:Float64}

    # Dimensions of fec all are [groups ⋅ sizes ⋅ locations]
    fec_all .= fec_params .* C_cover_t .* reshape(loc_area, (1, size(loc_area)...))
    # Sum over size classes
    @views fec_groups[:, :] .= dropdims(sum(fec_all; dims=2); dims=2)

    return nothing
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
function recruitment_rate(larval_pool::AbstractArray{T,2}, A::AbstractArray{T};
    α::Union{T,Vector{T}}=2.5, β::Union{T,Vector{T}}=5000.0)::Matrix{T} where {T<:Float64}
    sd = settler_density.(α, β, larval_pool) .* A'
    @views sd[sd .> 0.0] .= rand.(Poisson.(sd[sd .> 0.0]))

    return sd
end

"""
    settler_cover(fec_scope::T, conn::AbstractMatrix{Float64}, leftover_space::T, α::V, β::V, basal_area_per_settler::V, potential_settlers::T)::T where {T<:Matrix{Float64},V<:Vector{Float64}}

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
    conn::AbstractMatrix{Float64},
    leftover_space::V,
    α::V,
    β::V,
    basal_area_per_settler::V,
    potential_settlers::T,
    valid_sources::BitVector,
    valid_sinks::BitVector
)::Matrix{Float64} where {T<:AbstractMatrix{Float64},V<:Vector{Float64}}

    # Determine active sources and sinks
    valid_sources .= dropdims(sum(conn.data; dims=2) .> 0.0; dims=2)
    valid_sinks .= dropdims(sum(conn.data; dims=1) .> 0.0; dims=1)

    # Send larvae out into the world (reuse potential_settlers to reduce allocations)
    # Note, conn rows need not sum to 1.0 as this missing probability accounts for larvae
    # which do not settle. Pers comm with C. Ani (2023-01-29 13:24 AEST).
    # [Larval pool for each location in larvae/m²] * [survival rate]
    # this is known as in-water mortality.
    # Set to 0.0 as it is now taken care of by connectivity data.
    # Mwater::Float64 = 0.0
    @views @inbounds potential_settlers[:, valid_sinks] .= (
        fec_scope[:, valid_sources] * conn.data[valid_sources, valid_sinks]
    )

    # Larvae have landed, work out how many are recruited
    # Determine area covered by recruited larvae (settler cover) per m^2
    # recruits per m^2 per site multiplied by area per settler
    return recruitment_rate(potential_settlers, leftover_space; α=α, β=β) .*
           basal_area_per_settler
end
