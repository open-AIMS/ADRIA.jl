"""Coral growth functions"""


"""
    growthODE(du, X, p, _)

Base coral growth function.
"""
function growthODE(du::Array{Float64, 2}, X::Array{Float64, 2}, p::NamedTuple, _::Real)::Nothing
    k = @view p.k[:, :]
    k .= max.(p.P .- sum(X, dims=1), 0.0)

    # Use temporary caches
    k_X_r = @view p.kXr[:, :]
    k_rec = @view p.k_rec[:, :]
    X_mb = @view p.X_mb[:, :]
    kX_sel_en = @view p.kX_sel_en[:, :]
    X_tab = @view p.X_tab[:, :]
    @. k_X_r = k * X * p.r
    @. k_rec = k * p.rec
    @. X_mb = X * p.mb

    @views @. kX_sel_en = k * X[p.sel_en, :]
    @views @. X_tab = (p.mb[26] + p.comp * (X[6, :] + X[12, :])')

    r_comp = @views p.r[p.sel_en] .+ (p.comp .* sum(X[p.small_massives, :]))

    @views @. du[p.sel_en, :] = k_X_r[p.sel_en - 1, :] - kX_sel_en * r_comp - X_mb[p.sel_en, :]
    @views @. du[p.sel_unen, :] = kX_sel_en * r_comp + k_X_r[p.sel_unen, :] - X_mb[p.sel_unen, :]

    @views @. du[p.encrusting, :] = k_rec[p.enc, :] - k_X_r[p.encrusting, :] - X_mb[p.encrusting, :]

    @views @. du[p.small_massives, :] = k_X_r[p.small_massives - 1, :] - k_X_r[p.small_massives, :] - X[p.small_massives, :] * X_tab

    @views @. du[p.small, :] = k_rec[p.small_r, :] - k_X_r[p.small, :] - X_mb[p.small, :]
    @views @. du[p.mid, :] = k_X_r[p.mid - 1, :] - k_X_r[p.mid, :] - X_mb[p.mid, :]
    @views @. du[p.large, :] = k_X_r[p.large - 1 , :] + k_X_r[p.large, :] - X_mb[p.large, :]

    # Ensure no non-negative values
    du .= max.(du, 0.0)

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
    k = max(P - sum(X, 1), 0.0);

    # Total cover of small massives and encrusting
    X_sm = sum(X[26:28, :]);

    # Total cover of largest tabular Acropora
    X_tabular = (X[6, :] + X[12, :]); # this is for enhanced and unenhanced

    k_X_r = k .* X .* r;
    k_rec = k .* rec;
    X_mb = X .* mb;

    # Tabular Acropora Enhanced
    kX5 = k .* X[5, :]
    Y[5, :] = k_X_r[4, :] - kX5 .* (r[5] + comp .* X_sm) - X_mb[5, :];
    Y[6, :] = kX5 .* (r[5] + comp .* X_sm) + k_X_r[6, :] - X_mb[6, :];

    # Tabular Acropora Unenhanced
    kX11 = k .* X[11, :]
    Y[11, :] = k_X_r[10, :] - kX11 .* (r[11] + comp .* X_sm) - X_mb[11, :];
    Y[12, :] = kX11 .* (r[11] + comp .* X_sm) + k_X_r[12, :] - X_mb[12, :];

    # Corymbose Acropora Enhanced
    Y[13, :] = k_rec[3, :] - k .* X[13, :] .* r[13] - X_mb[13,:];

    # Small massives and encrusting Unenhanced
    Y[25, :] = k_rec[5, :] - k .* X[25, :] .* r[25] - X_mb[25,:];
    Y[26:28, :] = k_X_r[25:27, :] - k_X_r[26:28, :] - X[26:28, :] .* (mb[26] + comp .* X_tabular);

    # Small size classes
    Y[[1,7,19,31], :] = k_rec[[1,2,4,6], :] - k_X_r[[1,7,19,31], :] - X_mb[[1,7,19,31], :];

    # Mid size classes
    Y[[2:4,8:10,14:17,20:23,29,32:35], :] = k_X_r[[1:3,7:9,13:16,19:22,28,31:34], :] - k_X_r[[2:4,8:10,14:17,20:23,29,32:35], :] - X_mb[[2:4,8:10,14:17,20:23,29,32:35], :];

    # Larger size classes
    Y[[18,24,30,36], :] = k_X_r[[17,23,29,35], :] + k_X_r[[18,24,30,36], :] - X_mb[[18,24,30,36], :];

    # Ensure no non-negative values
    Y = max(Y, 0);
end


"""
    bleaching_mortality(tstep, n_p1, n_p2, a_adapt, n_adapt, bleach_resist, dhw)

Gompertz cumulative mortality function

Partial calibration using data by Hughes et al [1] (see Fig. 2C)

Parameters
----------
Y       : bleaching mortality for each coral species
tstep   : current time step
n_p1    : Gompertz distribution shape parameter 1
n_p2    : Gompertz distribution shape parameter 2
a_adapt : assisted adaptation
            where `sp` is the number of species considered
n_adapt : natural adaptation
            where `sp` is the number of species considered
dhw     : degree heating weeks for given time step for each site


References
----------
1. Hughes, T.P., Kerry, J.T., Baird, A.H., Connolly, S.R.,
    Dietzel, A., Eakin, C.M., Heron, S.F., Hoey, A.S.,
    Hoogenboom, M.O., Liu, G., McWilliam, M.J., Pears, R.J.,
    Pratchett, M.S., Skirving, W.J., Stella, J.S.
    and Torda, G. (2018)
    'Global warming transforms coral reef assemblages',
    Nature, 556(7702), pp. 492-496.
    doi:10.1038/s41586-018-0041-2.

2. Bozec, Y.-M. et. al. 2022 (in press). Cumulative impacts across
    Australia's Great Barrier Reef: A mechanistic evaluation.
    Ecological Monographs.
    https://doi.org/10.1101/2020.12.01.406413

3. Baird, A., Madin, J., Álvarez-Noriega, M., Fontoura, L.,
    Kerry, J., Kuo, C., Precoda, K., Torres-Pulliza, D., Woods, R.,
    Zawada, K., & Hughes, T. (2018).
    A decline in bleaching suggests that depth can provide a refuge
    from global warming in most coral taxa.
    Marine Ecology Progress Series, 603, 257-264.
    https://doi.org/10.3354/meps12732
"""
function bleaching_mortality!(Y::Array{Float64,2}, tstep::Int64, n_p1::Float64, n_p2::Float64,
    a_adapt::Vector{Float64}, n_adapt::Float64,
    bleach_resist::Vector{Float64}, dhw::Vector{Float64})::Nothing
    ad::Array{Float64} = a_adapt .+ bleach_resist .+ (tstep .* n_adapt)

    # Incorporate adaptation effect but maximum reduction is to 0
    capped::Array{Float64} = max.(0.0, dhw' .- ad)
    # Model 1: #Based on delta covers observed by Hughes et al. 2018 (Fig 2A)
    # and calibrated by Bozec et al. 2022
    # Halve bleaching mortality (see Baird et al., 2018)
    Y .= 1.0 .- exp.(n_p1 * (exp.(n_p2 * capped))) * 0.5

    return
end


"""
The scope that different coral groups and size classes have for
producing larvae without consideration of environment.

Coral fecundity per coral area of the different size classes.
When multiplied by the relative cover of each size class within taxa,
this produces an estimate of the relative fecundity of each coral group and size.
Total relative fecundity of a group is then calculated as the sum of
fecundities across size classes.

Parameters
----------
fec_groups : Matrix[n_classes, n_sites], memory cache to place results into
fec_all : Matrix[n_taxa, n_sites], temporary cache to place intermediate fecundity values into
fec_params : Vector, coral fecundity parameters
Y_pstep : Matrix[n_taxa, n_sites], of values in previous time step
site_area : Vector[n_sites], of site areas

Returns
-------
Matrix[n_classes, n_sites] : fecundity per m2 of coral
"""
function fecundity_scope!(fec_groups::Array{Float64, 2}, fec_all::Array{Float64, 2}, fec_params::Array{Float64},
                          Y_pstep::Array{Float64, 2}, site_area::Array{Float64})::Nothing
    ngroups::Int64 = size(fec_groups, 1)   # number of coral groups: 6
    nclasses::Int64 = size(fec_params, 1)  # number of coral size classes: 36

    fec_all .= fec_params .* Y_pstep .* site_area;
    for (i, (s, e)) in enumerate(zip(1:ngroups:nclasses, ngroups:ngroups:nclasses+1))
        @views fec_groups[i, :] = sum(fec_all[s:e, :], dims=1)
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


function larval_production(tstep, a_adapt, n_adapt, stresspast, LPdhwcoeff, DHWmaxtot, LPDprm2, n_groups)
    # Estimate how scope for larval production by each coral type changes as a
    # function of last year's heat stress. The function is theoretical and is
    # not yet verified by data. The rationale is that

    #
    # Inputs:
    #    tstep : int,
    #    assistadapt : array, DHW
    #    n_adapt : array, DHWs per year for all species
    #    stresspast : array, DHW at previous time step for each site
    #    LPdhwcoeff : float,
    #    DHWmaxtot : int, maximum DHW
    #    LPDprm2 : int, larval production parameter 2
    # Output:
    #    array of ngroups by nsites
    ad = @. a_adapt + tstep * n_adapt;

    # using half of DHWmaxtot as a placeholder
    # for the maximum capacity for thermal adaptation
    tmp_ad = @. (1 - ad / (DHWmaxtot/2));

    # One way around dimensional issue - tmp_ad for each class as the averaged
    # of the enhanced and unenhanced corals in that class
    # KA note: this works as it averages over size classes and not across groups.
    tmp_ad2 = mean(reshape(tmp_ad, Int(length(tmp_ad)/n_groups), n_groups), dims=1);

    return 1.0 .- exp.(-(exp.(-LPdhwcoeff .* (stresspast' .* tmp_ad2' .- LPDprm2))));
end
