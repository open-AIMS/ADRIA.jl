"""
    bleaching_mortality(tstep, n_p1, n_p2, a_adapt, n_adapt, bleach_resist, dhw)

Gompertz cumulative mortality function

Partial calibration using data by Hughes et al [1] (see Fig. 2C)

Parameters
----------
tstep   : int, current time step
n_p1    : float, Gompertz distribution shape parameter 1
n_p2    : float, Gompertz distribution shape parameter 2
a_adapt : array[sp*2, float], assisted adaptation
            where `sp` is the number of species considered
n_adapt : array[sp*2, float], natural adaptation
            where `sp` is the number of species considered
dhw     : float, degree heating weeks for given time step

Returns
-------
Y : Array[sp*2, float], bleaching mortality for each coral species

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
function bleaching_mortality(tstep, n_p1, n_p2, a_adapt, n_adapt, bleach_resist, dhw)
    ad = a_adapt + bleach_resist + (tstep .* n_adapt);

    # Incorporate adaptation effect but maximum reduction is to 0
    capped_dhw = max(0.0, dhw - ad);

    # Model 1: #Based on delta covers observed by Hughes et al. 2018 (Fig 2A)
    # and calibrated by Bozec et al. 2022
    Y = exp(n_p1 * (exp(n_p2 * capped_dhw)));

    # Halve bleaching mortality (see Baird et al., 2018)
    Y .= Y * 0.5;

end

