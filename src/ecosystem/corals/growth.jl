"""Coral growth functions"""


using Distributions


# include("growth_expanded.jl")


function growth_rate(linear_extension::Array{Float64}, diam_bin_widths::Array{Float64})::Vector{Float64}
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
function proportional_adjustment!(covers::Matrix{T}, cover_tmp::Vector{T}, max_cover::Array{T})::Nothing where {T<:Float64}
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
function growthODE(du::Array{Float64,2}, X::Array{Float64,2}, p::NamedTuple, _::Float64)::Nothing
    # Indices
    # p.small_massives := [26, 27, 28]
    # p.small := [1, 7, 13, 19, 25, 31]
    # p.mid := [2:4; 8:10; 14:17; 20:23; 29; 32:35]
    # p.large := [18, 24, 30, 36]
    # p.acr_5_11 := [5, 11]
    # p.acr_6_12 := [6, 12]

    # Intermediate values are now calculated outside of ODE function
    # To avoid repeat calculations
    # sXr : available space (sigma) * current cover (X) * growth rate (r)
    # X_mb : current cover (X) * background mortality (mb)
    # sm_comp : competition factor * area taken up by small massives (represents gain via competition with small massives)
    # M_sm : Mortality of small massives (background mortality + competition with acroporas)
    # rec : recruitment factors for each coral group (6 by n_sites)

    @views @. du[p.acr_5_11, :] = p.sXr[p.acr_5_11-1, :] - p.sXr[p.acr_5_11, :] + (p.sm_comp * X[p.acr_5_11, :]) - p.X_mb[p.acr_5_11, :]
    @views @. du[p.acr_6_12, :] = p.sXr[p.acr_6_12-1, :] + p.sXr[p.acr_6_12, :] + (p.sm_comp * X[p.acr_6_12, :]) - p.X_mb[p.acr_6_12, :]

    @views @. du[p.small_massives, :] = p.sXr[p.small_massives-1, :] - p.sXr[p.small_massives, :] - p.M_sm

    @views @. du[p.small, :] = p.rec - p.sXr[p.small, :] - p.X_mb[p.small, :]
    @views @. du[p.mid, :] = p.sXr[p.mid-1, :] - p.sXr[p.mid, :] - p.X_mb[p.mid, :]
    @views @. du[p.large, :] = p.sXr[p.large-1, :] + p.sXr[p.large, :] - p.X_mb[p.large, :]

    return
end


"""
    slow_ODE(Y, X, p, t)

Slow version of the growth model, ported directly from MATLAB.

Proportions of corals within a size class transitioning to the next size
class up (r) is based on the assumption that colony sizes within each size
bin are evenly distributed within bins. Transitions are then a simple
ratio of the change in colony size to the width of the bin. See
coralParms for further explanation of these coral metrics.

Note that recruitment pertains to coral groups (n = 6) and represents
the contribution to the cover of the smallest size class within each
group.  While growth and mortality metrics pertain to groups (6) as well
as size classes (6) across all sites (total of 36 by nsites), recruitment is
a 6 by nsites array.

Reshape flattened input from ODE back to expected matrix shape
Dims: (coral species, sites)
"""
function slow_ODE(Y, X, p, t)

    r, P, mb, rec, comp = p

    ## Density dependent growth and recruitment
    # P - sum over coral covers within each site
    # This sets the carrying capacity k := 0.0 <= k <= P
    # resulting in a matrix of (species * sites)
    # ensuring that density dependence is applied per site
    k = max(P - sum(X, 1), 0.0)

    # Total cover of small massives and encrusting
    X_sm = sum(X[26:28, :])

    # Total cover of largest tabular Acropora
    X_tabular = (X[6, :] + X[12, :]) # this is for enhanced and unenhanced

    k_X_r = k .* X .* r
    k_rec = k .* rec
    X_mb = X .* mb

    # Tabular Acropora Enhanced
    kX5 = k .* X[5, :]
    Y[5, :] = k_X_r[4, :] - kX5 .* (r[5] + comp .* X_sm) - X_mb[5, :]
    Y[6, :] = kX5 .* (r[5] + comp .* X_sm) + k_X_r[6, :] - X_mb[6, :]

    # Tabular Acropora Unenhanced
    kX11 = k .* X[11, :]
    Y[11, :] = k_X_r[10, :] - kX11 .* (r[11] + comp .* X_sm) - X_mb[11, :]
    Y[12, :] = kX11 .* (r[11] + comp .* X_sm) + k_X_r[12, :] - X_mb[12, :]

    # Corymbose Acropora Enhanced
    Y[13, :] = k_rec[3, :] - k .* X[13, :] .* r[13] - X_mb[13, :]

    # Small massives and encrusting Unenhanced
    Y[25, :] = k_rec[5, :] - k .* X[25, :] .* r[25] - X_mb[25, :]
    Y[26:28, :] = k_X_r[25:27, :] - k_X_r[26:28, :] - X[26:28, :] .* (mb[26:28] + comp .* X_tabular)

    # Small size classes
    Y[[1, 7, 19, 31], :] = k_rec[[1, 2, 4, 6], :] - k_X_r[[1, 7, 19, 31], :] - X_mb[[1, 7, 19, 31], :]

    # Mid size classes
    Y[[2:4, 8:10, 14:17, 20:23, 29, 32:35], :] = k_X_r[[1:3, 7:9, 13:16, 19:22, 28, 31:34], :] - k_X_r[[2:4, 8:10, 14:17, 20:23, 29, 32:35], :] - X_mb[[2:4, 8:10, 14:17, 20:23, 29, 32:35], :]

    # Larger size classes
    Y[[18, 24, 30, 36], :] = k_X_r[[17, 23, 29, 35], :] + k_X_r[[18, 24, 30, 36], :] - X_mb[[18, 24, 30, 36], :]

    # Ensure no non-negative values
    Y = max(Y, 0)
end


"""
    bleaching_mortality!(Y::Matrix{Float64}, tstep::Int64, depth::Vector{Float64},
        s::Vector{Float64}, dhw::Float64, a_adapt::Float64, n_adapt::Float64)

Calculates bleaching mortality taking into account depth and bleaching sensitivity of corals.
Model is adapted from Bozec et al., [2], itself based on data from Hughes et al., [3]
(bleaching sensitivity) and Baird et al., [1] (relationship between bleaching and depth).

# Arguments
- `Y` : Matrix to save results into
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
function bleaching_mortality!(Y::AbstractArray{Float64,2}, capped_dhw::AbstractArray{Float64,2}, depth_coeff::AbstractArray{Float64}, tstep::Int64, depth::Vector{Float64},
    s::Vector{Float64}, dhw::AbstractArray{Float64}, a_adapt::Vector{Float64}, n_adapt::Real)::Nothing

    # Incorporate adaptation effect but maximum reduction is to 0
    @. capped_dhw = ℯ^(0.17 + 0.35 * max(0.0, dhw' - (a_adapt + (tstep * n_adapt))))
    @. depth_coeff = ℯ^(-0.07551 * (depth - 2.0))

    # Estimate long-term bleaching mortality with an estimated depth coefficient and
    # initial bleaching mortality (models from Bozec et al., 2022)
    # Bozec et al., formulated the model to produce initial mortality (`m_init`) values
    # as a percentage (i.e., 0 - 100) and so we divide by 100 again to arrive at values 0 - 1.
    # m_init::Array{Float64} = min.(((depth_coeff .* s')' .* ℯ .^ (0.17 .+ 0.35 .* capped_dhw)) / 100.0 / 100.0, 1.0)

    # How much coral survives bleaching event
    # Y .= (1.0 .- m_init) .^ 6
    @. Y = (1.0 - min(((depth_coeff' * s) * capped_dhw) / 100.0 / 100.0, 1.0))^6.0

    return
end


"""
    fecundity_scope!(fec_groups::Array{Float64, 2}, fec_all::Array{Float64, 2}, fec_params::Array{Float64},
                     Y_pstep::Array{Float64, 2}, k_area::Array{Float64})::Nothing

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

# Returns
Matrix[n_classes, n_sites] : fecundity per m² of coral
"""
function fecundity_scope!(fec_groups::AbstractArray{T,2}, fec_all::AbstractArray{T,2}, fec_params::AbstractArray{T},
    Y_pstep::AbstractArray{T,2}, site_area::AbstractArray{T})::Nothing where {T<:Float64}
    ngroups::Int64 = size(fec_groups, 1)   # number of coral groups: 6
    nclasses::Int64 = size(fec_params, 1)  # number of coral size classes: 36

    fec_all .= fec_params .* Y_pstep .* site_area
    for (i, (s, e)) in enumerate(zip(1:ngroups:nclasses, ngroups:ngroups:nclasses+1))
        @views fec_groups[i, :] .= vec(sum(fec_all[s:e, :], dims=1))
    end

    # Above is equivalent to the below, but generic to any group/class size
    # @views fec_groups[1, :] = sum(fec_all[1:6, :], dims=1);   # Tabular Acropora enhanced
    # @views fec_groups[2, :] = sum(fec_all[7:12, :], dims=1);  # Tabular Acropora unenhanced
    # @views fec_groups[3, :] = sum(fec_all[13:18, :], dims=1); # Corymbose Acropora enhanced
    # @views fec_groups[4, :] = sum(fec_all[19:24, :], dims=1); # Corymbose Acropora unenhanced
    # @views fec_groups[5, :] = sum(fec_all[25:30, :], dims=1); # Small massives and encrusting
    # @views fec_groups[6, :] = sum(fec_all[31:36, :], dims=1); # Large massives

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
sf : Array of values ∈ [0,1] indicating reduced fecundity from a baseline.
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
    recruitment_rate(larval_pool, α=2.5, β=5000.0)

# Arguments
- `larval_pool` : Available larval pool
- `A` : Proportion of available area
- `α` : Maximum achievable density (settlers/m²) for a 100% free space
- `β` : Stock of larvae required to produce 50% of the maximum settlement

# Returns
λ, coral recruitment for each coral taxa based on a Poisson distribution.
"""
function recruitment_rate(larval_pool::AbstractArray{T,2}, A::AbstractArray{T}; α=2.5, β=5000.0)::AbstractArray{T} where {T<:Float64}
    sd = replace(settler_density.(α, β, larval_pool), Inf => 0.0, NaN => 0.0) .* A
    sel = sd .> 0.0
    sd[sel] .= rand.(Poisson.(sd[sel]))
    return sd
end


"""
    recruitment(larval_pool, A::Matrix{<:Real}; α=2.5, β=5000.0)

# Arguments
- `larval_pool` : Available larval pool
- `A` : Available space (0 - 1) relative to maximum area covered by
      cropped algal turf, i.e., the substratum that is suitable 
      for coral recruitment
- `α` : Maximum achievable density (settlers/m²) for a 100% free space
- `β` : Stock of larvae required to produce 50% of the maximum settlement

# Returns
Total coral recruitment for each coral taxa and site based on a Poisson distribution.
"""
function recruitment(larval_pool::AbstractArray{T,2}, A::Matrix{T}; α::Union{T,Vector{T}}=2.5, β::Union{T,Vector{T}}=5000.0)::Matrix{T} where {T<:Float64}
    # Minimum of recruited settler density (`recruitment_rate`) and max possible settler density (α)
    # return min.(recruitment_rate(larval_pool, A; α, β), α)
    return recruitment_rate(larval_pool, A; α, β)
end

"""
    settler_cover(fec_scope, sf, TP_data, leftover_space, α, β, basal_area_per_settler)

# Arguments
- `fec_scope` : fecundity scope
- `sf` : stressed fecundity
- `TP_data` : Transition probability
- `leftover_space` : difference between sites' maximum carrying capacity and current coral cover (k - C_s)
- `α` : max number of settlers / m²
- `β` : larvae / m² required to produce 50% of maximum settlement (default: 5000.0)
- `basal_area_per_settler` : area taken up by a single settler

# Returns
Area covered by recruited larvae (in m²)
"""
function settler_cover(fec_scope::T, sf::T,
    TP_data::T, leftover_space::Matrix{Float64},
    α::V, β::V, basal_area_per_settler::V)::Matrix{Float64} where {T<:AbstractArray{<:Float64,2},V<:Vector{Float64}}

    # Send larvae out into the world (reuse fec_scope to reduce allocations)
    # fec_scope .= (fec_scope .* sf)
    # fec_scope .= (fec_scope * TP_data) .* (1.0 .- Mwater)  # larval pool for each site (in larvae/m²)

    # As above, but more performant, less readable.
    Mwater::Float64 = 0.95
    mul!(fec_scope, (fec_scope .* sf), TP_data)
    fec_scope .= fec_scope .* (1.0 .- Mwater)

    # Larvae have landed, work out how many are recruited
    # recruits per m^2 per site multiplied by area per settler
    fec_scope .= recruitment(fec_scope, leftover_space; α=α, β=β) .* basal_area_per_settler

    # Determine area covered by recruited larvae (settler cover) per m^2
    return min.(fec_scope, leftover_space)
end
