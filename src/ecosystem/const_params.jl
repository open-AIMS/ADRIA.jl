"""
Struct of simulation constants for ADRIA


References:
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
"""
Base.@kwdef mutable struct SimConstants
    
    ## Base scenario parameters
    tf = 50; #number of years - e.g. year 2050 if we start deploying in year 2025 and run for 25 years.
    nsiteint = 5; # max number of sites we intervene on in a given year.
    
    # Default percent thresholds of max connectivity to filter out weak connections in connectivity network. 
    # Suggest we keep this low
    con_cutoff = 0.01;
    RCP = 45;  # RCP scenario to use
    prioritysites = []; # sites to prioritize when seeding or shading

    ## Environmental parameters
    beta = [1, 3]; # beta parameters for wave disturbance (distribution parameter)
    dhwmax25 = 5; # dhwmax at year 2025. NOTE: all warming simulations will change with new common DHW input for MDS team
    DHWmaxtot = 50; # max assumed DHW for all scenarios.  Will be obsolete when we move to new, shared inputs for DHW projections
    # wb1 = 0.55; # weibull parameter 1 for DHW distributions based on Lough et al 2018
    # wb2 = 2.24; # weibull parameter 2 for DHW distributions based on Lough et al 2018
    
    # Max total coral cover
    # used as a carrying capacity with 1-P representing space that is not
    # colonisable for corals
    max_coral_cover = 0.8;
    
    # Gompertz shape parameters 1 and 2 - for now applied to all coral species
    # equally. Based on Hughes et al 2018 and Bozec et al 2021.
    # Corrected to be consistent with zero bleaching mortality at DHW < 3.
    gompertz_p1 = 6.0;
    gompertz_p2 = 0.40;
    
    # Bleaching stress and coral fecundity parameters
    LPdhwcoeff = 0.4; # shape parameters relating dhw affecting cover to larval production
    LPDprm2 = 5; # parameter offsetting LPD curve
    
    # competition: probability that large tabular Acropora overtop small massives
    comp = 0.3;
    
    # True/False indicating Wwhether or not to mimic IPMF by loading only two coral types
    mimic_IPMF = Int8(0);  # Use 0 or 1 as booleans cannot be stored in netCDF

    max_settler_density = 2.5;                      # used by Bozec et al 2021 for Acropora
    density_ratio_of_settlers_to_larvae = 1 / 2000  # Bozec et al. 2021
    basal_area_per_settler = pi * ((1.25 / 100.0)^2)   # in m^2 assuming 1cm diameter
    # potential_settler_cover::Float64 = max_settler_density * basal_area_per_settler * density_ratio_of_settlers_to_larvae;
end
