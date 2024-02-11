"""Coral growth functions"""

using FLoops
using Distributions
using SpecialFunctions


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
    return vec(((2.0 .* linear_extension) ./ diam_bin_widths')')
end


"""
    proportional_adjustment!(coral_cover::Matrix{T}, cover_tmp::Vector{T})::Nothing where {T<:Float64}

Helper method to proportionally adjust coral cover, such that:
- `coral_cover` ∈ [0, 1].
- `sum(coral_cover, dims=1)` ≯  1.0

Modifies arrays in-place.

# Arguments
- `coral_cover` : Coral cover ∈ [0, 1]
- `cover_tmp` : temporary cache

# Returns
nothing
"""
function proportional_adjustment!(
    coral_cover::Union{SubArray{T},Matrix{T}}, cover_tmp::Vector{T}
)::Nothing where {T<:Float64}
    cover_tmp .= vec(sum(coral_cover; dims=1))
    cover_tmp[cover_tmp .≈ 1.0] .= 1.0
    if any(cover_tmp .> 1.0)
        exceeded::BitVector = vec(cover_tmp .> 1.0)
        msg = """
            Cover exceeded bounds, constraining to be within available space
            This indicates an issue with the model.
            Cover - Max Cover: $(sum(cover_tmp[exceeded] .- 1.0))
        """
        @debug msg
        @views @. coral_cover[:, exceeded] = (
            coral_cover[:, exceeded] / cover_tmp[exceeded]'
        )
    end

    coral_cover .= max.(coral_cover, 0.0)

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
    cover_tmp = zeros(size(coral_cover, 2))
    return proportional_adjustment!(coral_cover, cover_tmp)
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
function growthODE(du::Matrix{Float64}, X::Matrix{Float64}, p::NamedTuple, t::Real)::Nothing
    # Indices
    # small = [1, 7, 13, 19, 25, 31]
    # mid = [2:5; 8:11; 14:17; 20:23; 26:29; 32:35]
    # large = [6, 12, 18, 24, 30, 36]

    # sXr : available space (sigma) * current cover (X) * growth rate (r)
    # X_mb : current cover (X) * background mortality (mb)
    p.sXr .= (max.(1.0 .- sum(X, dims=1), 0.0) .* X .* p.r)
    p.X_mb .= X .* p.mb

    # For each size class, we determine the corals coming into size class due to growth,
    # and subtract those leaving the size class due to growth and background mortality
    # e.g., C = [corals coming in] - [corals going out]
    # The smallest size class only has corals leaving the size class.
    @views @. du[p.small, :] = -p.sXr[p.small, :] - p.X_mb[p.small, :]
    @views @. du[p.mid, :] = (p.sXr[p.mid-1, :] - p.X_mb[p.mid-1, :]) - (p.sXr[p.mid, :] + p.X_mb[p.mid, :])
    @views @. du[p.large, :] = (p.sXr[p.large-1, :] + p.sXr[p.large, :]) - p.X_mb[p.large, :]

    return nothing
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

    return nothing
end

"""
    bleaching_mortality!(cover::Matrix{Float64}, dhw::Vector{Float64},
        depth_coeff::Vector{Float64}, stdev::Vector{Float64}, dist_t_1::Matrix,
        dist_t::Matrix, prop_mort::SubArray{Float64})::Nothing

Applies bleaching mortality by assuming critical DHW thresholds are normally distributed for
all non-Juvenile (> 5cm diameter) size classes.

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
    depth_coeff::Vector{Float64}, stdev::Vector{Float64}, dist_t_1::Matrix{Float64},
    dist_t::Matrix{Float64}, prop_mort::SubArray{Float64}
)::Nothing
    n_sp_sc, n_locs = size(cover)

    # Determine non-juvenile size classes.
    # First two of each functional group are the juvenile size classes
    size_classes = 1:n_sp_sc
    juveniles = sort!(vec([size_classes[1:6:end] size_classes[2:6:end]]))
    non_juveniles = setdiff(size_classes, juveniles)

    # Adjust distributions for each functional group over all locations, ignoring juveniles
    # we assume the high background mortality of juveniles includes DHW mortality
    for loc in 1:n_locs
        # Skip locations that have no heat stress
        if dhw[loc] == 0.0
            continue
        end

        # Determine bleaching mortality for each non-juvenile species/size class
        for sp_sc in non_juveniles
            # Skip location if there is no population
            if cover[sp_sc, loc] == 0.0
                continue
            end

            μ::Float64 = dist_t_1[sp_sc, loc]
            affected_pop::Float64 = truncated_normal_cdf(
                # Use previous mortality threshold as minimum
                dhw[loc], μ, stdev[sp_sc], prop_mort[1, sp_sc, loc], μ + HEAT_UB
            )

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

            prop_mort[2, sp_sc, loc] = mort_pop
            if mort_pop > 0.0
                # Re-create distribution
                # Use same stdev as target size class to maintain genetic variance
                # pers comm K.B-N (2023-08-09 16:24 AEST)
                dist_t[sp_sc, loc] = truncated_normal_mean(
                    μ, stdev[sp_sc], mort_pop, μ + HEAT_UB
                )

                # Update population
                cover[sp_sc, loc] = cover[sp_sc, loc] * (1.0 - mort_pop)
            end
        end
    end

    clamp!(cover, 0.0, 1.0)

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
size class grows. Priors for the distributions are based on proportional cover and the
assumed growth rates for each size class.

i.e., \$(w_{i+1,1}, w_{i+1,2}) := (c_{i-1,i} / sum(c_{i-1,i})) * (g_{i-1,i} / sum(g_{i-1,i}))\$

where \$w\$ are the weights/priors and \$g\$ is the growth rates.

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
    for i in 6:-1:2
        # Skip size class if nothing is moving up
        sum(@view(cover[i-1:i])) == 0.0 ? continue : false

        prop_growth = @views (cover[i-1:i] ./ sum(cover[i-1:i])) .* (growth_rate[i-1:i] ./ sum(growth_rate[i-1:i]))

        dist_t[i] = sum(@view(dist_t[i-1:i]), Weights(prop_growth ./ sum(prop_growth)))
    end

    return nothing
end

"""
    adjust_DHW_distribution!(cover::SubArray{F}, n_groups::Int64, dist_t::Matrix{F},
        growth_rate::Matrix{F})::Nothing where {F<:Float64}

Adjust critical DHW thresholds for a given species/size class distribution as mortalities
affect the distribution over time, and corals mature (moving up size classes).

# Arguments
- `cover` : Coral cover (for timestep \$t-1\$ and \$t\$)
- `n_groups` : Number of coral groups represented
- `dist_t` : Distributions for timestep \$t\$
- `growth_rate` : Growth rates for each species/size class
"""
function adjust_DHW_distribution!(
    cover::SubArray{F},
    n_groups::Int64,
    dist_t::Matrix{F},
    growth_rate::Vector{F}
)::Nothing where {F<:Float64}
    _, n_sp_sc, n_locs = size(cover)

    step::Int64 = n_groups - 1

    # Adjust population distribution
    for sc1 in 1:n_groups:n_sp_sc
        for loc in 1:n_locs
            sc6::Int64 = sc1 + step

            # Skip if no cover
            if sum(@view(cover[1, sc1:sc6, loc])) == 0.0
                continue
            end

            # Combine distributions using a MixtureModel for all non-juvenile size
            # classes (we pass in all relevant size classes for the species group here).
            @views _shift_distributions!(cover[1, sc1:sc6, loc], growth_rate[sc1:sc6], dist_t[sc1:sc6, loc])
        end
    end

    return nothing
end

"""
    settler_DHW_tolerance!(c_mean_t_1::Matrix{F}, c_mean_t::Matrix{F}, k_area::Vector{F}, tp::AbstractMatrix{F}, settlers::Matrix{F}, fec_params_per_m²::Vector{F}, h²::F)::Nothing where {F<:Float64}

Update DHW tolerance means of newly settled corals into the distributions for the sink
locations.

# Arguments
- `c_mean_t_1` : DHW tolerance distribution for previous timestep (t - 1)
- `c_mean_t` : DHW tolerance distribution for current timestep (t)
- `k_area` : Absolute coral habitable area
- `tp` : Transition Probability matrix indicating the proportion of larvae that reaches a
    sink location (columns) from a given source location (rows)
    Columns should sum to 1
- `settlers` : recruited corals for each sink location
- `fec_params_per_m²` : Fecundity parameters for each species/size class combination
- `h²` : narrow-sense heritability
"""
function settler_DHW_tolerance!(
    c_mean_t_1::Matrix{F},
    c_mean_t::Matrix{F},
    k_area::Vector{F},
    tp::AbstractMatrix{F},
    settlers::Matrix{F},
    fec_params_per_m²::Vector{F},
    h²::F
)::Nothing where {F<:Float64}
    # Potential sink locations (TODO: pass in later)
    sink_loc_ids::Vector{Int64} = findall(k_area .> 0.0)

    source_locs::BitVector = falses(length(k_area))  # cache for source locations
    reproductive_sc::BitVector = falses(6)  # cache for reproductive size classes

    settler_sc::StepRange = 1:6:36
    for sink_loc in sink_loc_ids
        if sum(@views(settlers[:, sink_loc])) .== 0.0
            # Only update locations where recruitment occurred
            continue
        end

        # Note: source locations can include the sink location (self-seeding)
        source_locs .= @view(tp[:, sink_loc]) .> 0.0

        # Calculate contribution to cover to determine weights for each species/group
        w::NamedDimsArray = @views settlers[:, sink_loc]' .* tp[source_locs, sink_loc]
        w_per_group = w ./ sum(w, dims=1)
        replace!(w_per_group, NaN=>0.0)

        # Determine new distribution mean for each species at all locations
        for (sp, sc1) in enumerate(settler_sc)
            sc1_6::UnitRange{Int64} = sc1:sc1+5

            # Get distribution mean of reproductive size classes at source locations
            # recalling that source locations may include the sink location due to
            # self-seeding.
            reproductive_sc .= @view(fec_params_per_m²[sc1_6]) .> 0.0
            settler_means::SubArray{Float64} = @view(c_mean_t_1[sc1_6[reproductive_sc], source_locs])

            # Determine weights based on contribution to recruitment.
            # This weights the recruited corals by the size classes and source locations
            # which contributed to recruitment.
            if sum(@view(w_per_group[:, sp])) > 0.0
                ew::NamedDimsArray = repeat(@view(w_per_group[:, sp]), inner=count(reproductive_sc))

                # Determine combined mean
                # https://en.wikipedia.org/wiki/Mixture_distribution#Properties
                recruit_μ::Float64 = sum(settler_means, Weights(ew ./ sum(ew)))

                # Mean for generation t is determined through Breeder's equation
                c_mean_t[sc1, sink_loc] = breeders(c_mean_t_1[sc1, sink_loc], recruit_μ, h²)
            end
        end
    end

    return nothing
end


"""
    fecundity_scope!(fec_groups::Array{Float64, 2}, fec_all::Array{Float64, 2},
                     fec_params::Array{Float64}, C_t::Array{Float64, 2},
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
- `C_t` : Matrix[n_taxa, n_sites], of coral cover values for the previous time step
- `site_area` : Vector[n_sites], total site area in m²
"""
function fecundity_scope!(
    fec_groups::AbstractArray{T,2},
    fec_all::AbstractArray{T,2},
    fec_params::AbstractArray{T},
    C_t::AbstractArray{T,2},
    site_area::AbstractArray{T}
)::Nothing where {T<:Float64}
    ngroups::Int64 = size(fec_groups, 1)   # number of coral groups: 6
    nclasses::Int64 = size(fec_params, 1)  # number of coral size classes: 36

    fec_all .= fec_params .* C_t .* site_area
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
    sd = replace(settler_density.(α, β, larval_pool), Inf => 0.0, NaN => 0.0) .* A
    @views sd[sd.>0.0] .= rand.(Poisson.(sd[sd.>0.0]))
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
    leftover_space::T,
    α::V,
    β::V,
    basal_area_per_settler::V,
    potential_settlers::T
)::T where {T<:Matrix{Float64},V<:Vector{Float64}}

    # Determine active sources and sinks
    valid_sources::BitVector = vec(sum(conn, dims=2) .> 0.0)
    valid_sinks::BitVector = vec(sum(conn, dims=1) .> 0.0)

    # Send larvae out into the world (reuse potential_settlers to reduce allocations)
    # Note, conn rows need not sum to 1.0 as this missing probability accounts for larvae
    # which do not settle. Pers comm with C. Ani (2023-01-29 13:24 AEST).
    # [Larval pool for each location in larvae/m²] * [survival rate]
    # this is known as in-water mortality.
    # Set to 0.0 as it is now taken care of by connectivity data.
    # Mwater::Float64 = 0.0
    @views potential_settlers[:, valid_sinks] .= (
        fec_scope[:, valid_sources] * conn[valid_sources, valid_sinks]
    )

    # Larvae have landed, work out how many are recruited
    # Determine area covered by recruited larvae (settler cover) per m^2
    # recruits per m^2 per site multiplied by area per settler
    return recruitment_rate(potential_settlers, leftover_space; α=α, β=β) .* basal_area_per_settler
end
