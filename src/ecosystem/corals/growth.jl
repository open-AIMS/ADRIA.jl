"""Coral growth functions"""


using Distributions


# include("growth_expanded.jl")


"""
    growth_rate(linear_extension::Array{Float64}, diam_bin_widths::Array{Float64})::Vector{Float64}

Determine the rate of growth representing the proportion of each size class that moves
up a size class each (yearly) time step. Values > 1 indicate transitions to higher size
classes occurs more than once per time step.

# Arguments
- `linear_extension` : Linear extension in cm/year
- `diam_bin_widths` : diameter of each size class (bin) in cm

# Returns
Vector, of size \$N = [n_{species} ⋅ n_{classes}]\$ indicating proportional growth rates
for each.
"""
function growth_rate(linear_extension::Matrix{Float64}, diam_bin_widths::Vector{Float64})::Vector{Float64}
    return vec(((2.0 * linear_extension) ./ diam_bin_widths')')
end


"""
    proportional_adjustment!(covers::Matrix{T}, cover_tmp::Vector{T}, max_cover::Array{T})::Nothing where {T<:Real}

Helper method to proportionally adjust coral cover.
Modifies arrays in-place.

# Arguments
- `covers` : Coral cover result set
- `cover_tmp` : Temporary cache matrix used to hold sum over species. Avoids memory allocations
- `max_cover` : Maximum possible coral cover for each site
"""
function proportional_adjustment!(covers::Union{SubArray{T},Matrix{T}}, cover_tmp::Vector{T}, max_cover::Array{T})::Nothing where {T<:Float64}
    # Proportionally adjust initial covers
    cover_tmp .= vec(sum(covers, dims=1))
    if any(cover_tmp .> max_cover)
        exceeded::Vector{Int64} = findall(cover_tmp .> max_cover)
        @views @. covers[:, exceeded] = (covers[:, exceeded] / cover_tmp[exceeded]') * max_cover[exceeded]'
    end

    covers .= max.(covers, 0.0)

    return
end


"""
    growthODE(du, X, p, _)

Base coral growth function.

Proportion of corals within a size class are modeled to transition to the
next size class up. Assumes colony sizes are evenly distributed within each
size bin. Transitions are a ratio of the change in colony size to the width
of the bin. See `coral_spec()` for further detail.

Note that recruitment pertains to coral groups (\$n = 6\$) and represents
the contribution to the cover of the smallest size class within each
group.  While growth and mortality metrics pertain to groups (6) as well
as size classes (6) across all sites (total of 36 by \$n_sites\$), recruitment is
a 6 by \$n_sites\$ array.
"""
function growthODE(du::Array{Float64,2}, X::Array{Float64,2}, p::NamedTuple, _::Real)::Nothing
    # Indices

    # small = [1, 7, 13, 19, 25, 31]
    # mid = [2:5; 8:11; 14:17; 20:23; 26:29; 32:35]
    # large = [6, 12, 18, 24, 30, 36]

    # Intermediate values are now calculated outside of ODE function
    # To avoid repeat calculations
    # sXr : available space (sigma) * current cover (X) * growth rate (r)
    # X_mb : current cover (X) * background mortality (mb)
    # rec : recruitment factors for each coral group (6 by n_sites)

    @views @. du[p.small, :] = p.rec - p.sXr[p.small, :] - p.X_mb[p.small, :]
    @views @. du[p.mid, :] = p.sXr[p.mid-1, :] - p.sXr[p.mid, :] - p.X_mb[p.mid, :]
    @views @. du[p.large, :] = p.sXr[p.large-1, :] + p.sXr[p.large, :] - p.X_mb[p.large, :]

    return
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
    @. capped_dhw = min.(ℯ^(0.17 + 0.35 * max(0.0, dhw' - (a_adapt + (tstep * n_adapt)))) - 1.0, 100.0)
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

    return
end

"""
    bleaching_mortality!(cover::Matrix{Float64}, dhw::Vector{Float64},
        depth_coeff::Vector{Float64}, stdev::Vector{Float64}, dist_t_1::Matrix,
        dist_t::Matrix, prop_mort::SubArray{Float64})::Nothing

Applies bleaching mortality by assuming critical DHW thresholds are normally distributed for
all non-Juvenile (> 5cm diameter) size classes. Distributions are informed by learnings from
Bairos-Novak et al., [1] and (unpublished) data referred to in Hughes et al., [2]. Juvenile
mortality is assumed to be primarily represented by other factors (i.e., background
mortality; see Álvarez-Noriega et al., [3]). The proportion of the population which bleached
is estimated with the Cumulative Density Function. Bleaching mortality is then estimated
with a depth-adjusted coefficient (from Baird et al., [4]).

# Arguments
- `cover` : Coral cover for current timestep
- `dhw` : DHW for all represented locations
- `depth_coeff` : Pre-calculated depth coefficient for all locations
- `stdev` : Standard deviation of DHW tolerance
- `dist_t_1` : Critical DHW threshold distribution for current timestep, for all species and
           locations
- `dist_t` : Critical DHW threshold distribution for next timestep, for all species and
              locations
- `prop_mort` : Cache to store records of bleaching mortality

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
function bleaching_mortality!(cover::Matrix{Float64}, dhw::Vector{Float64},
    depth_coeff::Vector{Float64}, stdev::Vector{Float64}, dist_t_1::Matrix{Distribution},
    dist_t::Matrix{Distribution}, prop_mort::SubArray{Float64})::Nothing
    n_sp_sc, n_locs = size(cover)

    # Adjust distributions for all locations, ignoring juveniles
    # we assume the high background mortality of juveniles includes DHW mortality
    @floop for (sp_sc, loc) in Iterators.product(3:n_sp_sc, 1:n_locs)
        # Skip if location experiences no heat stress or there is no population
        if dhw[loc] == 0.0 || cover[sp_sc, loc] == 0.0
            continue
        end

        affected_pop::Float64 = cdf(dist_t_1[sp_sc, loc], dhw[loc])
        mort_pop::Float64 = 0.0
        if affected_pop > 0.0
            # Calculate depth-adjusted bleaching mortality
            mort_pop = affected_pop * depth_coeff[loc]

            # Set values close to 0.0 (e.g., 1e-214) to 0.0
            # https://github.com/JuliaLang/julia/issues/23376#issuecomment-324649815
            if (mort_pop + one(mort_pop)) ≈ one(mort_pop)
                mort_pop = 0.0
            end
        end

        prop_mort[sp_sc, loc] = mort_pop
        if mort_pop > 0.0
            # Re-create distribution
            d::Distribution = dist_t_1[sp_sc, loc]
            μ::Float64 = mean(d)
            dist_t[sp_sc, loc] = truncated(Normal(μ, stdev[sp_sc]), mort_pop, μ + HEAT_UB)

            # Update population
            cover[sp_sc, loc] = cover[sp_sc, loc] * (1.0 - mort_pop)
        end
    end
end

"""
    _shift_distributions!(cover::SubArray, growth_rate::SubArray, dist_t::SubArray, stdev::SubArray)::Nothing

Combines distributions between size classes > 1 to represent the shifts that occur as each
size class grows. Priors for the distributions are based on proportional cover and the
assumed growth rates for each size class.

i.e., (w_{i+1,1}, w_{i+1,2}) := (c_{i-1,i} / sum(c_{i-1,i})) * (g_{i-1,i} / sum(g_{i-1,i}))

where \$w\$ are the weights/priors and \$g\$ is the growth rates.

# Arguments
- `cover` : Coral cover for \$t-1\$
- `growth_rate` : Growth rates for the given size classes/species
- `dist_t` : Critical DHW threshold distribution for timestep \$t\$
- `stdev` : Standard deviations of coral DHW tolerance
- `tmp` : Reused cache for MixtureModels to avoid allocations
"""
function _shift_distributions!(cover::SubArray, growth_rate::SubArray, dist_t::SubArray, stdev::SubArray)::Nothing
    # Weight distributions based on growth rate and cover
    @floop for i in 6:-1:3
        sum(cover[i-1:i]) == 0.0 ? continue : false
        prop_growth = (cover[i-1:i] ./ sum(cover[i-1:i])) .* (growth_rate[i-1:i] ./ sum(growth_rate[i-1:i]))
        x::MixtureModel = MixtureModel([dist_t[i-1], dist_t[i]], prop_growth ./ sum(prop_growth))

        # Use same stdev as target size class to maintain genetic variance
        # pers comm K.B-N (2023-08-09 16:24 AEST)
        dist_t[i] = truncated(Normal(mean(x), stdev[i]), minimum(x), mean(x) + HEAT_UB)
    end

    # # Size class 1 all moves up to size class 2 anyway
    μ = mean(dist_t[1])
    dist_t[2] = truncated(Normal(μ, stdev[2]), minimum(dist_t[1]), μ + HEAT_UB)

    return nothing
end

"""
    adjust_DHW_distribution!(cover::SubArray, n_groups::Int64, dist_t_1::SubArray,
        dist_t::SubArray, growth_rate::SubArray, stdev::SubArray, h²::Float64)::Nothing

Adjust critical DHW thresholds for a given species/size class distribution as mortalities
affect the distribution over time, and corals mature (moving up size classes).

# Arguments
- `cover` : Coral cover (for timestep \$t-1\$ and \$t\$)
- `n_groups` : Number of coral groups represented
- `dist_t_1` : Distributions for timestep \$t-1\$
- `dist_t` : Distributions for timestep \$t\$
- `growth_rate` : Growth rates for each species/size class
- `stdev` : standard deviations of DHW tolerances for each size class
- `h²` : heritability value
"""
function adjust_DHW_distribution!(cover::SubArray, n_groups::Int64, dist_t_1::Matrix{Distribution},
    dist_t::Matrix{Distribution}, growth_rate::Matrix{Float64}, stdev::Vector{Float64}, h²::Float64)::Nothing
    _, n_sp_sc, n_locs = size(cover)

    step::Int64 = n_groups - 1
    weights::Vector{Float64} = zeros(3)

    # Adjust population distribution
    @floop for (sc1, loc) in Iterators.product(1:n_groups:n_sp_sc, 1:n_locs)
        sc4::Int64 = sc1 + 3
        sc6::Int64 = sc1 + step

        # Skip if no cover
        if sum(cover[1, sc1:sc6, loc]) == 0.0
            continue
        end

        # Combine distributions using a MixtureModel for all size
        # classes >= 4 (the correct size classes are selected within the
        # function).
        @views _shift_distributions!(cover[1, sc1:sc6, loc], growth_rate[sc1:sc6], dist_t[sc1:sc6, loc], stdev[sc1:sc6])

        # Reproduction. Size class >= 3 spawn larvae.
        # A new distribution for size class 1 is then determined by taking the
        # distribution of size classes >= 3 at time t+1 and size class 1 at time t, and
        # applying the S⋅h² calculation, where:
        # - $S$ is the distance between the means of the gaussian distributions
        # - $h$ is heritability (assumed to range from 0.25 to 0.5, nominal value of 0.3)
        dt_1::Truncated{Normal{Float64},Continuous,Float64,Float64,Float64} = dist_t_1[sc1, loc]
        if any(cover[2, sc4:sc6, loc] .!= 0.0)
            @views weights .= cover[2, sc4:sc6, loc] / sum(cover[2, sc4:sc6, loc])

            d_t::MixtureModel = MixtureModel(Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}[dist_t[sc4:sc6, loc]...], weights)

            S::Float64 = mean(d_t) - mean(dt_1)
            if S != 0.0
                # The new distribution mean for size class 1 is then: prev mean + (S⋅h²)
                # The standard deviation is assumed to be the same as parents
                # Note: Nominally, h² := 0.3, but ranges from 0.25 - 0.5
                μ_t::Float64 = mean(dt_1) + (S * h²)
                dist_t[sc1, loc] = truncated(Normal(μ_t, stdev[sc1]), minimum(d_t), μ_t + HEAT_UB)
            end
        end
    end

    return nothing
end


"""
    fecundity_scope!(fec_groups::Array{Float64, 2}, fec_all::Array{Float64, 2},
                     fec_params::Array{Float64}, Y_pstep::Array{Float64, 2},
                     k_area::Array{Float64})::Nothing

The scope that different coral groups and size classes have for
producing larvae without consideration of environment.

Coral fecundity per coral area of the different size classes.
When multiplied by the relative cover of each size class within taxa,
this produces an estimate of the relative fecundity of each coral group and size.
Total relative fecundity of a group is then calculated as the sum of
fecundities across size classes.

# Arguments
- `fec_groups` : Matrix[n_classes, n_sites], memory cache to place results into
- `fec_all` : Matrix[n_taxa, n_sites], temporary cache to place intermediate fecundity values into
- `fec_params` : Vector, coral fecundity parameters (in per m²) for each species/size class
- `Y_pstep` : Matrix[n_taxa, n_sites], of coral cover values for the previous time step
- `site_area` : Vector[n_sites], total site area in m²
"""
function fecundity_scope!(fec_groups::AbstractArray{T,2}, fec_all::AbstractArray{T,2}, fec_params::AbstractArray{T},
    Y_pstep::AbstractArray{T,2}, site_area::AbstractArray{T})::Nothing where {T<:Float64}
    ngroups::Int64 = size(fec_groups, 1)   # number of coral groups: 6
    nclasses::Int64 = size(fec_params, 1)  # number of coral size classes: 36

    fec_all .= fec_params .* Y_pstep .* site_area
    for (i, (s, e)) in enumerate(zip(1:ngroups:nclasses, ngroups:ngroups:nclasses+1))
        @views fec_groups[i, :] .= vec(sum(fec_all[s:e, :], dims=1))
    end
    return nothing
end


"""
    stressed_fecundity(tstep, a_adapt, n_adapt, stresspast, LPdhwcoeff, DHWmaxtot, LPDprm2, n_groups)

Estimate how scope for larval production by each coral type changes as a
function of last year's heat stress. The function is theoretical and is not
yet verified by data.

Stressed Fecundity (sf) is based on the proportion of baseline fecundity
that is unaffected by heat stress in the previous year - e.g., a value
of 0.9 inside sf(i, j) indicates that species i at site j can only produce
90% of its usual larval output.

# Arguments
- `tstep` : Current time step
- `a_adapt` : DHW reduction of enhanced corals
- `n_adapt` : DHWs reduction (linearly scales with `tstep`)
- `stresspast` : DHW at previous time step for each site
- `LPdhwcoeff` :
- `DHWmaxtot` : Maximum DHW
- `LPDprm2` : Larval production parameter 2
- `n_groups` : Number of groups

# Returns
`sf` : Array of values ∈ [0,1] indicating reduced fecundity from a baseline.
"""
function stressed_fecundity(tstep::Int64, a_adapt::Vector{T}, n_adapt::T,
    stresspast::Vector{T}, LPdhwcoeff::T, DHWmaxtot::T, LPDprm2::T, n_groups::Int64)::Matrix{T} where {T<:Float64}
    ad::Vector{Float64} = @. a_adapt + tstep * n_adapt

    # using half of DHWmaxtot as a placeholder
    # for the maximum capacity for thermal adaptation
    tmp_ad::Vector{Float64} = @. (1.0 - (ad / (DHWmaxtot / 2.0)))

    # One way around dimensional issue - tmp_ad for each class as the averaged
    # of the enhanced and unenhanced corals in that class
    # KA note: this works as it averages over size classes and not across groups.
    tmp_ad2::Vector{Float64} = vec(mean(reshape(tmp_ad, Int64(length(tmp_ad) / n_groups), n_groups), dims=1))

    return 1.0 .- exp.(.-(exp.(-LPdhwcoeff .* (stresspast' .* tmp_ad2 .- LPDprm2))))
end


"""
    settler_density(α, β, L)

Note for β: "For corals, the actual number of 6-month old recruits for each coral group
    is generated [...] following a Poisson distribution with recruitment event rate λ.

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

# Examples
settler_density(2.5, 5000.0, L)
"""
function settler_density(α, β, L)
    return (α .* L) ./ (β .+ L)
end


"""
    recruitment_rate(larval_pool::AbstractArray{T,2}, A::AbstractArray{T}; α=2.5, β=5000.0)

Calculates coral recruitment for each species/group and location.

# Arguments
- `larval_pool` : Available larval pool
- `A` : Available space (0 - 1) relative to maximum area covered by
      cropped algal turf, i.e., the substratum that is suitable
      for coral recruitment
- `α` : Maximum achievable density (settlers/m²) for a 100% free space
- `β` : Stock of larvae required to produce 50% of the maximum settlement

# Returns
λ, total coral recruitment for each coral taxa and location based on a Poisson distribution.
"""
function recruitment_rate(larval_pool::AbstractArray{T,2}, A::AbstractArray{T};
    α::Union{T,Vector{T}}=2.5, β::Union{T,Vector{T}}=5000.0)::AbstractArray{T} where {T<:Float64}
    sd = replace(settler_density.(α, β, larval_pool), Inf => 0.0, NaN => 0.0) .* A
    sd[sd.>0.0] .= rand.(Poisson.(sd[sd.>0.0]))
    return sd
end


"""
    settler_cover(fec_scope, sf, TP_data, leftover_space, α, β, basal_area_per_settler)

Determine area settled by recruited larvae.

Note: Units for all areas are assumed to be in m².

# Arguments
- `fec_scope` : fecundity scope
- `TP_data` : Transition probability
- `leftover_space` : difference between sites' maximum carrying capacity and current coral cover (\$k - C_s\$)
- `α` : max number of settlers / m²
- `β` : larvae / m² required to produce 50% of maximum settlement
- `basal_area_per_settler` : area taken up by a single settler

# Returns
Area covered by recruited larvae (in m²)
"""
function settler_cover(fec_scope::T,
    TP_data::T, leftover_space::Matrix{Float64},
    α::V, β::V, basal_area_per_settler::V)::Matrix{Float64} where {T<:AbstractArray{<:Float64,2},V<:Vector{Float64}}

    # Send larvae out into the world (reuse fec_scope to reduce allocations)
    # fec_scope .= (fec_scope .* sf)
    # fec_scope .= (fec_scope * TP_data) .* (1.0 .- Mwater)  # larval pool for each site (in larvae/m²)

    # As above, but more performant, less readable.
    Mwater::Float64 = 0.95
    fec_scope *= TP_data
    fec_scope .*= (1.0 .- Mwater)

    # Larvae have landed, work out how many are recruited
    # recruits per m^2 per site multiplied by area per settler
    fec_scope .= recruitment_rate(fec_scope, leftover_space; α=α, β=β) .* basal_area_per_settler

    # Determine area covered by recruited larvae (settler cover) per m^2
    return min.(fec_scope, leftover_space)
end
