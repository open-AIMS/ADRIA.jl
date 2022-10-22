"""Coral growth functions"""


using Distributions

include("growth_expanded.jl")


function growth_rate(linear_extension, diam_bin_widths, mean_colony_diameter_m)
    bin_shift = ((2.0 * linear_extension'[:]) ./ diam_bin_widths) .^ 2
    bin_shift[bin_shift.>1] .= 1
    class_6_diameter_ratio = (mean_colony_diameter_m[:, end] .+ (diam_bin_widths[end] / (2 * 100))) ./ mean_colony_diameter_m[:, end]
    mean_diameter_ratio = hcat((mean_colony_diameter_m[:, 2:end] ./ mean_colony_diameter_m[:, 1:end-1]) .^ 2, class_6_diameter_ratio .^ 2)

    return (mean_diameter_ratio'[:] .* bin_shift)
end


"""
    proportional_adjustment!(Yout::AbstractArray{<:Real}, cover_tmp::AbstractArray{<:Real}, max_cover::AbstractArray{<:Real})

Helper method to proportionally adjust coral cover.
Modifies arrays in-place.

# Arguments
- Yout : Coral cover result set
- cover_tmp : Temporary cache matrix used to hold sum over species. Avoids memory allocations
- max_cover : maximum possible coral cover for each site
"""
function proportional_adjustment!(Yout::AbstractArray{<:Real}, cover_tmp::AbstractArray{<:Real}, max_cover::AbstractArray{<:Real})
    # Proportionally adjust initial covers
    cover_tmp .= vec(sum(Yout, dims=1))
    if any(cover_tmp .> max_cover)
        exceeded::Vector{Int32} = findall(cover_tmp .> max_cover)

        @views Yout[:, exceeded] .= (Yout[:, exceeded] ./ cover_tmp[exceeded]') .* max_cover[exceeded]'
    end

    return max.(Yout, 0.0)
end


"""
    growthODE(du, X, p, _)

Base coral growth function.
"""
function growthODE(du::Array{Float64,2}, X::Array{Float64,2}, p::NamedTuple, _::Real)::Nothing
    s = @view p.sigma[:, :]

    # X is cover relative to `k` (max. carrying capacity)
    # So we subtract from 1.0 to get leftover/available space, relative to `k`
    s .= max.(1.0 .- sum(X, dims=1), 0.0)
    s[p.k'.==0.0] .= 0.0

    # Indices
    # p.small_massives := [26, 27, 28]
    # p.small := [1, 7, 13, 19, 25, 31]
    # p.mid := [2:4; 8:10; 14:17; 20:23; 29; 32:35]
    # p.large := [18, 24, 30, 36]
    # p.acr_5_11 := [5, 11]
    # p.acr_6_12 := [6, 12]

    # Use temporary caches
    sXr = @view p.sXr[:, :]
    X_mb = @view p.X_mb[:, :]
    M_sm = @view p.M_sm[:, :]
    r_comp = @view p.r_comp[:, :]

    @. sXr = s * X * p.r  # leftover space * current cover * growth_rate
    @. X_mb = X * p.mb    # current cover * background mortality

    @views @. M_sm = X[p.small_massives, :] * (p.mb[p.small_massives] + p.comp * (X[6, :] + X[12, :])')

    r_comp .= p.comp .* sum(X[p.small_massives, :], dims=1)
    @views @. du[p.acr_5_11, :] = sXr[p.acr_5_11-1, :] - sXr[p.acr_5_11, :] + r_comp * X[p.acr_5_11] - X_mb[p.acr_5_11, :]
    @views @. du[p.acr_6_12, :] = sXr[p.acr_6_12-1, :] + sXr[p.acr_6_12, :] + r_comp * X[p.acr_6_12] - X_mb[p.acr_6_12, :]

    @views @. du[p.small_massives, :] = sXr[p.small_massives-1, :] - sXr[p.small_massives, :] - M_sm

    @views @. du[p.small, :] = s * p.rec - sXr[p.small, :] - X_mb[p.small, :]
    @views @. du[p.mid, :] = sXr[p.mid-1, :] - sXr[p.mid, :] - X_mb[p.mid, :]
    @views @. du[p.large, :] = sXr[p.large-1, :] + sXr[p.large, :] - X_mb[p.large, :]

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
        s::Vector{Float64}, dhw::Float64, a_adapt::Float64, n_adapt::Float64,
        bleach_resist::Vector{Float64})

Calculates bleaching mortality taking into account depth and bleaching sensitivity of corals.
Model is adapted from Bozec et al., [2], itself based on data from Hughes et al., [3]
(bleaching sensitivity) and Baird et al., [1] (relationship between bleaching and depth).

# Arguments
- Y : Matrix to save results into
- tstep : current time step
- depth : mean site depth (m) for each site
- s : bleaching sensitivity of corals (relative values) for each taxa/size class
- dhw : Degree Heating Week experienced at site
- a_adapt : Level of assisted adaptation (DHW reduction)
- n_adapt : Level of natural adaptation (DHW reduction linearly scaled over time)
- bleach_resist : Level of bleaching resistance (inherent resilience to DHW)

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
function bleaching_mortality!(Y::Matrix{Float64}, tstep::Int64, depth::Vector{Float64},
    s::Vector{Float64}, dhw::Vector{Float64}, a_adapt::Vector{Float64}, n_adapt::Float64,
    bleach_resist::Vector{Float64})::Nothing

    # Incorporate adaptation effect but maximum reduction is to 0
    ad::Array{Float64} = a_adapt .+ bleach_resist .+ (tstep .* n_adapt)
    capped_dhw::Array{Float64} = max.(0.0, dhw' .- ad)

    # Estimate long-term bleaching mortality with an estimated depth coefficient and
    # initial bleaching mortality (models from Bozec et al., 2022)
    # `m_init` as initially formulated produces values as a percentage (i.e., 0 - 100)
    # and so we divide by 100 again to arrive at values 0 - 1.
    depth_coeff = ℯ .^ (-0.07551 .* (depth .- 2.0))
    m_init = min.(((depth_coeff .* s')' .* ℯ .^ (0.17 .+ 0.35 .* capped_dhw)) / 100.0 / 100.0, 1.0)

    # How much coral survives bleaching event
    Y .= (1.0 .- m_init) .^ 6

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
- fec_groups : Matrix[n_classes, n_sites], memory cache to place results into
- fec_all : Matrix[n_taxa, n_sites], temporary cache to place intermediate fecundity values into
- fec_params : Vector, coral fecundity parameters (in per m²) for each species/size class
- Y_pstep : Matrix[n_taxa, n_sites], of coral cover values for the previous time step
- site_area : Vector[n_sites], total site area in m²

# Returns
Matrix[n_classes, n_sites] : fecundity per m² of coral
"""
function fecundity_scope!(fec_groups::Array{Float64,2}, fec_all::Array{Float64,2}, fec_params::Array{Float64},
    Y_pstep::Array{Float64,2}, site_area::Array{Float64})::Nothing
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
- tstep : int,
- a_adapt : array, DHW
- n_adapt : array, DHWs per year for all species
- stresspast : array, DHW at previous time step for each site
- LPdhwcoeff : float,
- DHWmaxtot : int, maximum DHW
- LPDprm2 : int, larval production parameter 2
- n_groups : int, number of groups

# Returns
sf : Array of values ∈ [0,1] indicating reduced fecundity from a baseline.
"""
function stressed_fecundity(tstep, a_adapt, n_adapt, stresspast, LPdhwcoeff, DHWmaxtot, LPDprm2, n_groups)::Matrix{Float64}
    ad::Vector{Float64} = @. a_adapt + tstep * n_adapt

    # using half of DHWmaxtot as a placeholder
    # for the maximum capacity for thermal adaptation
    tmp_ad::Vector{Float64} = @. (1.0 - (ad / (DHWmaxtot / 2)))

    # One way around dimensional issue - tmp_ad for each class as the averaged
    # of the enhanced and unenhanced corals in that class
    # KA note: this works as it averages over size classes and not across groups.
    tmp_ad2::Array{Float64} = vec(mean(reshape(tmp_ad, Int64(length(tmp_ad) / n_groups), n_groups), dims=1))

    return 1.0 .- exp.(-(exp.(-LPdhwcoeff .* (stresspast' .* tmp_ad2 .- LPDprm2))))
end


"""
    settler_density(α, β, L)

Note for β: "For corals, the actual number of 6-month old recruits for each coral group
    is generated [...] following a Poisson distribution with recruitment event rate λ.

# Arguments
- α : maximum achievable density (settlers/m²) for a 100% free space (set to 2.5 in [1] for Corymbose)
- β : stock of larvae required to produce one-half the maximum settlement (larvae/m²),
        i.e., α/2(m²), set to 5000 in [1].
- L : available larval pool

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
- larval_pool : Available larval pool
- α : maximum achievable density (settlers/m²) for a 100% free space
- β : stock of larvae required to produce 50% of the maximum settlement

# Returns
λ, coral recruitment for each coral taxa based on a Poisson distribution.
"""
function recruitment_rate(larval_pool; α=2.5, β=5000.0)
    sd = replace(settler_density.(α, β, larval_pool), Inf => 0.0, NaN => 0.0)
    sd[sd.>0.0] .= rand.(Poisson.(sd[sd.>0.0]))
    return sd
end


"""
    recruitment(larval_pool, A::Matrix{<:Real}; α=2.5, β=5000.0)

# Arguments
- larval_pool : Available larval pool
- A : proportional space (in m²) covered by cropped algal turf,
        i.e., the substratum that is suitable for coral recruitment
- α : maximum achievable density (settlers/m²) for a 100% free space
- β : stock of larvae required to produce 50% of the maximum settlement

# Returns
Total coral recruitment for each coral taxa and site based on a Poisson distribution.
"""
function recruitment(larval_pool, A::Matrix{<:Real}; α=2.5, β=5000.0)
    # Minimum of recruited settler density (`recruitment_rate`) and max possible settler density (α)
    return min.(recruitment_rate(larval_pool; α, β), α) .* A
end

"""
    settler_cover(fec_scope, sf, TP_data, leftover_space, max_density, basal_area_per_settler)

# Arguments
- fec_scope : fecundity scope
- sf : stressed fecundity
- TP_data : Transition probability
- leftover_space : difference between sites' maximum carrying capacity and current coral cover (k - C_s)
- max_density : number of settlers / m²
- basal_area_per_settler : area taken up by a single settler

# Returns
Area covered by recruited larvae (in m²)
"""
function settler_cover(fec_scope, sf, TP_data, leftover_space, max_density, basal_area_per_settler)
    # Send larvae out into the world
    actual_fecundity = (fec_scope .* sf)

    Mwater = 0.95
    larval_pool = (actual_fecundity * TP_data) .* (1 - Mwater)  # larval pool for each site (in larvae/m²)

    # β is stock of larvae required to produce 50% of the maximum settlement
    β = replace((max_density .* leftover_space) ./ 2.0, Inf => 0.0, NaN => 0.0)

    # Larvae have landed, work out how many are recruited
    λ = recruitment(larval_pool, leftover_space; α=max_density, β=β)  # recruits per m^2 per site

    # Determine area covered by recruited larvae (settler cover) and constrain to available space
    return min.(λ .* basal_area_per_settler, leftover_space)
end
