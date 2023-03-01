"""
    SimConstants

Struct of simulation constants for ADRIA

# References
1. Lough, J. M., Anderson, K. D., & Hughes, T. P. (2018).
   Increasing thermal stress for tropical coral reefs: 1871-2017.
   Scientific Reports, 8(1), 6079.
   https://doi.org/10.1038/s41598-018-24530-9

2. Hughes, T. P., Kerry, J. T., Baird, A. H., Connolly, S. R.,
     Dietzel, A., Eakin, C. M., Heron, S. F., Hoey, A. S.,
     Hoogenboom, M. O., Liu, G., McWilliam, M. J., Pears, R. J.,
     Pratchett, M. S., Skirving, W. J., Stella, J. S., & Torda, G. (2018).
   Global warming transforms coral reef assemblages.
   Nature, 556(7702), 492-496.
   https://doi.org/10.1038/s41586-018-0041-2

3. Bozec, Y.-M., Rowell, D., Harrison, L., Gaskell, J., Hock, K.,
     Callaghan, D., Gorton, R., Kovacs, E. M., Lyons, M., Mumby, P.,
     & Roelfsema, C. (2021).
   Baseline mapping to support reef restoration and
     resilience-based management in the Whitsundays.
   https://doi.org/10.13140/RG.2.2.26976.20482

4. Bozec, Y.-M., Hock, K., Mason, R. A. B., Baird, M. E., Castro-Sanguino, C.,
     Condie, S. A., Puotinen, M., Thompson, A., & Mumby, P. J. (2022).
   Cumulative impacts across Australia's Great Barrier Reef: A mechanistic evaluation. 
   Ecological Monographs, 92(1), e01494.
   https://doi.org/10.1002/ecm.1494
"""
Base.@kwdef mutable struct SimConstants
    n_site_int = 5 # max number of sites we intervene on in a given year.

    priority_sites = [] # sites to prioritize when seeding or shading
    priority_zones = ["DarkBlue", "Green"] # zones to prioritize when seeding or shading, in order of preference
    seed_criteria_names = ["heat_stress", "wave_stress", "in_conn", "out_conn", "cover_space", "zones", "predec"]
    shade_criteria_names = ["heat_stress", "wave_stress", "in_conn", "cover", "zones", "predec"]

    ## Environmental parameters
    # 50 DHW approximates the highest predicted value for the century for SSPs 3 and 5.
    DHWmaxtot = 50.0

    # Bleaching stress and coral fecundity parameters
    LPdhwcoeff = 0.4  # shape parameters relating dhw affecting cover to larval production
    LPDprm2 = 5.0  # parameter offsetting LPD curve

    # competition: probability that large Tabular Acropora overtop small massives
    comp = 0.3

    # True/False indicating Wwhether or not to mimic IPMF by loading only two coral types
    # currently unused
    mimic_IPMF = Int8(0)  # Use 0 or 1 as booleans cannot be stored in netCDF

    # Used by Bozec et al 2021 for Acropora, small massives and large massives
    # max_settler_density = [2.5, 2.5, 2.5, 2.5, 1.5, 1.5]

    # Modified for use in ADRIA
    max_settler_density = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0]

    max_larval_density = [5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0]
end
