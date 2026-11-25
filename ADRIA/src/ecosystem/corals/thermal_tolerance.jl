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
    effective_dhw_at_depth(dhw_surface, depth, κ_base, mixing_scale)

Compute effective DHW at a given depth using Beer-Lambert attenuation with
mixing-dependent attenuation decay.

At low surface DHW, thermal stratification is intact and deeper sites experience
substantially less thermal stress. At high surface DHW, convective/wind mixing
erodes stratification, reducing the effective attenuation coefficient towards zero
(well-mixed conditions) so that mortality approaches the surface rate at extreme DHW.

# Arguments
- `dhw_surface` : Surface (or reference depth ~2m) DHW
- `depth` : Site depth in meters
- `κ_base` : Base attenuation coefficient (m⁻¹) in the `dhw_surface -> 0` limit;
    Free parameter, not empirically anchored — see note below
- `mixing_scale` : DHW at which mixing halves the effective attenuation coefficient.
    Free parameter; treat as a scenario axis rather than a fitted value. Zero is the
    permanently well-mixed limit, in which depth provides no refuge at any DHW

Both are exposed as model factors by [`DepthAttenuation`](@ref) and are read from the
scenario `param_set` in `run_model`.

# Notes
Both parameters are unconstrained by Baird et al. [1], despite that reference being
cited here for the qualitative depth-refuge effect. They report no attenuation
coefficient; their only quantitative depth anchor is ~50% fewer bleached colonies at
12m than at 2m, which `κ_base = 0.04` does not reproduce (matching it needs 0.056-0.137
depending on assumed surface DHW). Their assemblage model is also logit-linear in
log10(depth) — a power law — so no exponential coefficient exists to extract, and the
implied κ varies with anchor depth (0.145 at 6m, 0.050 at 27m). They found no
significant depth × DHW interaction over DHW 4-10, so `mixing_scale` is likewise
unconstrained by that dataset.

# References
1. Baird, A., Madin, J., Álvarez-Noriega, M., Fontoura, L., Kerry, J., Kuo, C.,
     Precoda, K., Torres-Pulliza, D., Woods, R., Zawada, K., & Hughes, T. (2018).
   A decline in bleaching suggests that depth can provide a refuge from global
     warming in most coral taxa.
   Marine Ecology Progress Series, 603, 257-264.
   https://doi.org/10.3354/meps12732
"""
function effective_dhw_at_depth(
    dhw_surface::Float64, depth::Float64, κ_base::Float64, mixing_scale::Float64
)::Float64
    # `mixing_scale = 0` is the permanently well-mixed limit: no attenuation at any DHW.
    # Taken as a special case because the general expression is 0/0 when `dhw_surface` is
    # also zero, which would otherwise put a NaN into `bleaching_mortality!`.
    κ_eff = mixing_scale > 0.0 ? κ_base / (1.0 + dhw_surface / mixing_scale) : 0.0
    Δz = max(0.0, depth - 2.0)
    return dhw_surface * exp(-κ_eff * Δz)
end

"""
    ps_effective_dhw_at_depth(dhw_surface, depth; R, ζ₁, ζ₂, mixing_scale)

Compute effective DHW at a given depth using a dual-band exponential decay model
(Paulson & Simpson 1977), which separately accounts for infrared and visible-light
heat penetration.

The standard Beer-Lambert single-exponential overestimates heat penetration at depth
because it cannot simultaneously represent two physically distinct regimes:
- **Infrared** (~58% of solar energy) is absorbed within the first ~0.35 m. Deeper
  sites never receive this fraction regardless of water clarity or mixing.
- **Visible light** (~42%) penetrates to 5–25 m depending on water clarity (Jerlov
  water type) and is the dominant heat source for coral habitat.

At high surface DHW, wind/convective mixing erodes thermal stratification, redistributing
heat downward. This is captured by scaling `ζ₂` with DHW, so that under extreme events
the visible-band effective scale depth increases and deeper sites experience more stress.

DHW at depth is computed relative to the SST reference depth (2 m):

    fraction(z) = [R·exp(-z/ζ₁) + (1−R)·exp(-z/ζ₂_eff)] /
                  [R·exp(-z_ref/ζ₁) + (1−R)·exp(-z_ref/ζ₂_eff)]

    effective_dhw = dhw_surface × clamp(fraction, 0, 1)

# Arguments
- `dhw_surface` : Surface DHW at the SST reference depth (~2 m)
- `depth` : Site depth in metres
- `R` : Infrared fraction of solar energy (default 0.58; range 0.50–0.62)
- `ζ₁` : Infrared e-folding depth in metres (default 0.35); fixed — property of water,
    not affected by mixing
- `ζ₂` : Visible-light e-folding depth in metres under calm conditions (default 15.0,
    Jerlov Type II — typical of mid-shelf GBR); range 5 m (turbid inshore) to 25 m
    (clear oceanic)
- `mixing_scale` : DHW at which mixing doubles `ζ₂` (default 12.0)

# Notes
This function models accumulated thermal stress (DHW), not instantaneous solar flux.
The dual-band formulation is physically motivated by how solar energy is deposited at
depth; temperature anomaly propagation via mixing dynamics is approximated through the
`mixing_scale` parameter.

This method is hypothetical and untested. It requires further empirical data for comparison
and assessment.

# References
1. Paulson, C. A., & Simpson, J. J. (1977).
   Irradiance measurements in the upper ocean.
   Journal of Physical Oceanography, 7(6), 952–956.
   https://doi.org/10.1175/1520-0485(1977)007<0952:IMITUO>2.0.CO;2
2. Baird, A., Madin, J., Álvarez-Noriega, M., Fontoura, L., Kerry, J., Kuo, C.,
   Precoda, K., Torres-Pulliza, D., Woods, R., Zawada, K., & Hughes, T. (2018).
   A decline in bleaching suggests that depth can provide a refuge from global
   warming in most coral taxa.
   Marine Ecology Progress Series, 603, 257–264.
   https://doi.org/10.3354/meps12732
"""
function ps_effective_dhw_at_depth(
    dhw_surface::Float64, depth::Float64;
    R::Float64=0.58,            # infrared fraction of solar energy
    ζ₁::Float64=0.35,           # infrared e-folding depth (m) — fixed, water property
    ζ₂::Float64=15.0,           # visible e-folding depth (m) — Jerlov Type II
    mixing_scale::Float64=12.0  # DHW at which mixing doubles ζ₂
)::Float64
    z_ref = 2.0   # Reference depth of SST measurement
    ζ₂_eff = ζ₂ * (1.0 + dhw_surface / mixing_scale)  # Mixing increases visible penetration

    I_ref = R * exp(-z_ref / ζ₁) + (1.0 - R) * exp(-z_ref / ζ₂_eff)
    I_z = R * exp(-depth / ζ₁) + (1.0 - R) * exp(-depth / ζ₂_eff)

    fraction = I_z / I_ref
    return dhw_surface * clamp(fraction, 0.0, 1.0)
end

"""
    bleaching_mortality!(
        cover::AbstractArray{Float64,3},
        eff_dhw::Vector{Float64},
        stdev::AbstractMatrix{Float64},
        dist_t::AbstractArray{Float64,3},
        bleach_dhw::SubArray{Float64},
        tol_ceil::AbstractArray{Float64,3};
        mean_tol_update_method::Symbol=:blend,
        low_dhw_reduction::Float64=0.0,
        low_dhw_threshold::Float64=0.0
    )::Nothing

Applies bleaching mortality by assuming critical DHW thresholds are normally distributed for
all size classes, as colony size is a weak predictor of bleaching susceptibility
(Álvarez-Noriega et al., [3]).

Distributions are informed by learnings from Bairos-Novak et al., [1] and (unpublished)
data referred to in Hughes et al., [2]. The proportion of the population which bleached is
estimated with the Cumulative Density Function. Bleaching mortality is then estimated by
incorporating light absorption principles (with the Beer-Lambert Law) as a proxy for heat
dissipation at depths. The approach loosely aligns with the depth-adjusted coefficient from
Baird et al., [4].

# Arguments
- `cover` : Coral cover for current timestep
- `eff_dhw` : Effective DHW for all represented locations at their depths
- `stdev` : Standard deviation of DHW tolerance
- `dist_t` : Critical DHW threshold distribution for current timestep, for all species and
    locations
- `bleach_dhw` : Cache to store records of bleaching mortality (stores mortalities for
    t-1 and t)
- `tol_ceil` : Upper limit on the mean DHW tolerance
- `mean_tol_update_method` : How the mean tolerance of a bleached size class is updated.
    - `:blend` (default): average of the pre-bleaching mean and the survivors' mean,
      weighted by the number of individuals before and after bleaching, i.e.
      `(μ + (1 - m) * μ_surv) / (2 - m)` for mortality fraction `m`. Damps the large
      first-event jump of `:legacy`.
    - `:legacy`: the survivors' mean replaces the mean, and is re-used as the centre of
      the (re-truncated) distribution at the next bleaching event.
- `low_dhw_reduction` : Fraction in [0, 1] by which bleaching mortality is reduced at
    locations whose effective DHW is below `low_dhw_threshold` (see `LowDHWBleaching`).
    Mortality is scaled by `1 - low_dhw_reduction`. Values outside [0, 1] throw an
    `ArgumentError`. The default 0.0 leaves mortality unchanged.
- `low_dhw_threshold` : Effective DHW (degree heating weeks) below which the reduction
    applies. The default 0.0 means it never applies.

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

3. Álvarez-Noriega, M., Aston, E.A., Becker, M., Fabricius, K.E., Figueira, W.F.,
     Gordon, S.E., Krensel, R., Lechene, M.A., Remmers, T., Toor, M., Ferrari, R., 2025.
   Challenging Paradigms Around the Role of Colony Size, Taxa, and Environment on
     Bleaching Susceptibility.
   Global Change Biology 31, e70090.
   https://doi.org/10.1111/gcb.70090

4. Baird, A., Madin, J., Álvarez-Noriega, M., Fontoura, L., Kerry, J., Kuo, C.,
     Precoda, K., Torres-Pulliza, D., Woods, R., Zawada, K., & Hughes, T. (2018).
   A decline in bleaching suggests that depth can provide a refuge from global
     warming in most coral taxa.
   Marine Ecology Progress Series, 603, 257-264.
   https://doi.org/10.3354/meps12732

# See also
- `effective_dhw_at_depth()`
"""
function bleaching_mortality!(
    cover::AbstractArray{Float64,3},
    eff_dhw::Vector{Float64},
    stdev::AbstractMatrix{Float64},
    dist_t::AbstractArray{Float64,3},
    bleach_dhw::SubArray{Float64},
    tol_ceil::AbstractArray{Float64,3};
    mean_tol_update_method::Symbol=:blend,
    low_dhw_reduction::Float64=0.0,
    low_dhw_threshold::Float64=0.0
)::Nothing
    mean_tol_update_method in (:blend, :legacy) ||
        throw(ArgumentError("Unknown mean_tol_update_method: $mean_tol_update_method"))
    # Optional reduction of bleaching mortality where effective DHW is below
    # `low_dhw_threshold` (see `LowDHWBleaching`). The defaults leave mortality unchanged.
    0.0 <= low_dhw_reduction <= 1.0 ||
        throw(ArgumentError("low_dhw_reduction must lie in [0, 1]"))
    n_groups, n_sizes, _ = size(cover)

    # New evidence that all size classes bleach
    # Álvarez-Noriega et al., 2025.
    # Challenging Paradigms Around the Role of Colony Size, Taxa, and Environment on
    # Bleaching Susceptibility.
    # Global Change Biology 31, e70090.
    # https://doi.org/10.1111/gcb.70090
    all_sizes = 1:n_sizes

    # Potential bleaching locations (skip locations with no heat stress)
    # Use HEAT_LB (heat tolerance distribution lower bound) as the bleaching threshold
    active_locs = findall(eff_dhw .> HEAT_LB)

    # Adjust distributions for each functional group over all locations and size classes
    for loc in active_locs
        # Determine bleaching mortality for each species/size class
        for grp = 1:n_groups
            # Skip location if there is no population
            if sum(@view(cover[grp, :, loc])) == 0.0
                continue
            end

            for sc in all_sizes
                # Skip location if there is no population
                if cover[grp, sc, loc] == 0.0
                    continue
                end

                μ::Float64 = dist_t[grp, sc, loc]
                μ_ceil::Float64 = tol_ceil[grp, sc, loc]  # initial mean + HEAT_UB (fixed)
                affected_pop::Float64 = truncated_normal_cdf(
                    # Use the previous bleaching DHW as the distribution lower bound,
                    # with HEAT_LB (1.0 DHW-weeks) as the minimum susceptibility threshold
                    eff_dhw[loc], μ, stdev[grp, sc],
                    max(HEAT_LB, bleach_dhw[1, grp, sc, loc]),
                    μ_ceil
                )

                mort_pop::Float64 = 0.0
                if affected_pop > 0.0
                    # Depth effect is now handled via effective_dhw_at_depth (see model_phases.jl).
                    mort_pop = affected_pop
                    if low_dhw_reduction > 0.0 && eff_dhw[loc] < low_dhw_threshold
                        mort_pop *= 1.0 - low_dhw_reduction
                    end

                    # Set values close to 0.0 (e.g., 1e-214) to 0.0
                    # https://github.com/JuliaLang/julia/issues/23376#issuecomment-324649815
                    if (mort_pop + one(mort_pop)) ≈ one(mort_pop)
                        mort_pop = 0.0
                    end
                end

                # Store current DHW (DHW-weeks) as the bleaching threshold for this location.
                # Previously stored mort_pop (a dimensionless fraction 0-1), which caused a
                # dimensional mismatch: the value was used as a lower bound on the DHW
                # tolerance distribution (units: DHW-weeks) in the next timestep, making
                # the truncation effectively identical to truncating at 0.
                bleach_dhw[2, grp, sc, loc] = eff_dhw[loc]
                if mort_pop > 0.0
                    # 1. Re-create distribution truncated at current bleaching DHW: survivors
                    # must have had tolerance above eff_dhw[loc], so this is the correct lower
                    # bound.
                    #
                    # 2. Use same stdev as target size class to maintain genetic variance
                    # pers comm K.B-N (2023-08-09 16:24 AEST)
                    μ_surv::Float64 = truncated_normal_mean(
                        μ, stdev[grp, sc], max(HEAT_LB, eff_dhw[loc]), μ_ceil
                    )
                    dist_t[grp, sc, loc] = new_mean_heat_stress_tolerance(
                        μ, μ_surv, mort_pop; mean_tol_update_method
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
    new_mean_heat_stress_tolerance(
        μ::Float64, μ_surv::Float64, mort_pop::Float64; mean_tol_update_method::Symbol
    )::Float64

Mean heat stress tolerance of a size class after a bleaching event.

# Arguments
- `μ` : Mean tolerance before bleaching
- `μ_surv` : Mean tolerance of the survivors (truncated-normal mean)
- `mort_pop` : Fraction of the population lost to bleaching
- `mean_tol_update_method` : `:legacy` returns `μ_surv`; `:blend` returns the mean of `μ` and
    `μ_surv` weighted by the number of individuals before and after bleaching.
"""
function new_mean_heat_stress_tolerance(
    μ::Float64, μ_surv::Float64, mort_pop::Float64; mean_tol_update_method::Symbol
)::Float64
    return if mean_tol_update_method === :legacy
        μ_surv
    elseif mean_tol_update_method === :blend
        (μ + (1.0 - mort_pop) * μ_surv) / (2.0 - mort_pop)
    else
        throw(ArgumentError("Unknown mean_tol_update_method: $mean_tol_update_method"))
    end
end

"""
    breeders(μ_base::T, μ_target::T, h²::T)::T where {T<:Float64}

Shift a tolerance-distribution mean a fraction `h²` of the way from `μ_base` towards
`μ_target`, per the (univariate) Breeder's equation.

```
S = μ_target - μ_base
R = S * h²
μ_new = μ_base + R
```

# Arguments
- `μ_base` : Baseline tolerance-distribution mean (the value being shifted)
- `μ_target` : Tolerance-distribution mean being selected towards
- `h²` : Narrow-sense heritability

# Returns
`μ_base` shifted a fraction `h²` of the way towards `μ_target`.

# Notes
Used in `settler_DHW_tolerance!` to blend, for each source location/functional
group/size class, that population's own mean tolerance from the previous timestep
(`μ_base`) with its mean this timestep (`μ_target`). Each source population's heritable
response to selection is applied locally, before its larvae disperse via connectivity —
not after mixing tolerance means across the connectivity network.
"""
function breeders(μ_base::T, μ_target::T, h²::T)::T where {T<:Float64}
    return μ_base + ((μ_target - μ_base) * h²)
end

"""
    _shift_distributions!(
        cover::SubArray{F},
        growth_rate::SubArray{F},
        dist_t::SubArray{F}
    )::Nothing where {F<:Float64}

Combines distributions between size classes > 1 to represent the shifts that occur as each
size class grows. Weights for the distributions are based on proportional cover and the
assumed growth rates for each size class.

i.e., \$(w_{i+1,1}, w_{i+1,2}) := (c_{i-1,i} / sum(c_{i-1,i})) * (g_{i-1,i} / sum(g_{i-1,i}))\$

where \$w\$ are the weights and \$g\$ is the growth rates.

# Arguments
- `cover` : Coral cover for \$t-1\$
- `growth_rate` : Growth rates for the given size classes/species
- `dist_t` : Critical DHW threshold distribution for timestep \$t\$
"""
function _shift_distributions!(
    cover::SubArray{T},
    growth_rate::SubArray{T},
    dist_t::SubArray{T}
)::Nothing where {T}
    # Weight distributions based on growth rate and cover
    # Do from largest size class to size class 2
    # (values for size class 1 gets replaced by recruitment process)
    for i = length(growth_rate):-1:2
        # Pre-compute sums once
        cover_sum = sum(view(cover, (i - 1):i))

        # Skip size class if nothing is moving up
        cover_sum == 0.0 ? continue : false

        growth_sum = sum(view(growth_rate, (i - 1):i))

        if cover_sum > 0.0 && growth_sum > 0.0
            # Calculate weights directly
            w1 = (cover[i - 1] / cover_sum) * (growth_rate[i - 1] / growth_sum)
            w2 = (cover[i] / cover_sum) * (growth_rate[i] / growth_sum)

            cache_sum = w1 + w2
            if cache_sum > 0.0
                dist_t[i] = (dist_t[i - 1] * w1 + dist_t[i] * w2) / cache_sum
            end
        end
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
        for loc = 1:n_locs
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
"""
    adjust_DHW_distribution!(cover_t_1, dist_t, growth_rate, tol_ceil)

Adjust critical DHW distributions across all groups and locations as corals mature
into higher size classes. Applies a hard ceiling so tolerance cannot exceed
`tol_ceil` (initial population mean + HEAT_UB) regardless of size-class redistribution.

# Arguments
- `cover_t_1`  : Coral cover at the previous timestep (groups × sizes × locations).
- `dist_t`     : Tolerance distribution means for the current timestep; updated in place.
- `growth_rate`: Growth rates per (group, size, location).
- `tol_ceil`   : Per-(group, size, location) tolerance ceiling (initial mean + HEAT_UB).
"""
function adjust_DHW_distribution!(
    cover_t_1::SubArray{T,3},
    dist_t::AbstractArray{T,3},
    growth_rate::AbstractArray{T,3},
    tol_ceil::AbstractArray{T,3}
)::Nothing where {T<:Float64}
    groups, _, locs = axes(cover_t_1)

    for grp in groups
        for loc in locs
            if sum(@view(cover_t_1[grp, :, loc])) == 0.0
                continue
            end

            @views _shift_distributions!(
                cover_t_1[grp, :, loc],
                growth_rate[grp, :, loc],
                dist_t[grp, :, loc]
            )
        end
    end

    dist_t .= min.(dist_t, tol_ceil)
    return nothing
end

"""
    settler_DHW_tolerance!(
        c_mean_t_1::AbstractArray{F,3},
        c_mean_t::AbstractArray{F,3},
        C_cover_t::AbstractArray{F,3},
        k_area::Vector{F},
        tp::SparseMatrixCSC{F},
        settlers::AbstractMatrix{F},
        fecundity_per_m²::AbstractMatrix{F},
        h²::F
    )::Nothing where {F<:Float64}

# Arguments
- `c_mean_t_1` : Mean of heat tolerance distribution at t-1 for each functional group,
size class and location.
- `c_mean_t` : Mean of heat tolerance distribution at t for each functional group,
size class and location.
- `C_cover_t` : Coral cover at timestep t.
- `k_area` : Coral habitable available at each location (`k` value) in either relative or
absolute units.
- `tp` : Connectivity matrix.
- `settlers` : Recruitment cover matrix for each functional group (rows) and locations (cols).
- `fecundity_per_m²` : Fecundity in number of corals per area (in m²) each functional group.
- `h²`      : Heritability parameter.
- `tol_ceil`: Per-(group, size, location) tolerance ceiling (initial mean + HEAT_UB).
  Applied as a hard cap to `c_mean_t` after settler tolerances are written, preventing
  connectivity-driven drift from compounding across timesteps beyond the initial ceiling.
"""
function settler_DHW_tolerance!(
    c_mean_t_1::AbstractArray{F,3},
    c_mean_t::AbstractArray{F,3},
    C_cover_t::AbstractArray{F,3},
    k_area::Vector{F},
    tp::SparseMatrixCSC{F},
    settlers::AbstractMatrix{F},
    fecundity_per_m²::AbstractMatrix{F},
    h²::F,
    tol_ceil::AbstractArray{F,3}
)::Nothing where {F<:Float64}
    groups, _, locs = axes(c_mean_t_1)
    n_groups = length(groups)
    n_locs = length(locs)

    # Use Breeder's equation to calculate mean heat tolerance of offspring for each loc/grp
    # Reproductive size classes differ between functional groups
    fecund_idx_per_group = [findall(fecundity_per_m²[grp, :] .> 0.0) for grp in groups]
    max_n_fecund = maximum(length, fecund_idx_per_group)
    c_mean_settlers = zeros(n_groups, n_locs)       # will be the next c_mean_t[:,1,:]

    # Caches for the offspring c_mean and the cover-proportion weights of the reproductive
    # size classes (reused each iteration to avoid reallocating)
    c_mean_per_size_buf::Vector{Float64} = zeros(Float64, max_n_fecund)
    w_buf::Vector{Float64} = zeros(Float64, max_n_fecund)

    pos_settlers_loc_idx = findall(vec(sum(settlers; dims=1)) .> 0)
    for loc in pos_settlers_loc_idx
        for grp in groups
            fecund_idx = fecund_idx_per_group[grp]
            isempty(fecund_idx) && continue

            @views cover_col = C_cover_t[grp, fecund_idx, loc]
            cover_sum = sum(cover_col)
            (cover_sum == 0) ? continue : 0

            w_view = @view w_buf[1:length(fecund_idx)]
            c_mean_per_size = @view c_mean_per_size_buf[1:length(fecund_idx)]
            w_view .= cover_col ./ cover_sum
            w = StatsBase.weights(w_view)
            @views c_mean_per_size .= breeders.(
                c_mean_t_1[grp, fecund_idx, loc],
                c_mean_t[grp, fecund_idx, loc],
                h²
            )
            c_mean_settlers[grp, loc] = StatsBase.mean(c_mean_per_size, w)
        end
    end

    sink_loc_ids::Vector{Int64} = findall(k_area .> 0.0)

    # Number of reproductive size classes for each group
    n_reproductive::Vector{Int64} = vec(sum(fecundity_per_m² .> 0.0; dims=2))

    rows = rowvals(tp)
    vals = nonzeros(tp)

    @inbounds for sink_loc in sink_loc_ids
        sink_settlers = @view(settlers[:, sink_loc])
        if sum(sink_settlers) == 0.0
            # Only update locations where recruitment occurred
            continue
        end

        col_range = nzrange(tp, sink_loc)
        if isempty(col_range)
            continue
        end

        # Calculate sum of weights (connectivity strength) to normalize
        # We perform the summation over the column manually to avoid allocating a view
        sum_tp::F = 0.0
        for k in col_range
            sum_tp += vals[k]
        end

        if sum_tp == 0.0
            continue
        end

        inv_sum_tp::F = 1.0 / sum_tp

        for grp in groups
            if sink_settlers[grp] == 0.0
                continue
            end

            recruit_μ::F = 0.0

            # Fused loop: Iterate over sources (non-zeros in TP column)
            for k in col_range
                src = rows[k]
                w = vals[k] * inv_sum_tp
                recruit_μ += c_mean_settlers[grp, src] * w
            end

            c_mean_t[grp, 1, sink_loc] = recruit_μ
        end
    end

    c_mean_t .= min.(c_mean_t, tol_ceil)
    return nothing
end
