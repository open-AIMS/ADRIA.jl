using ADRIA
using ADRIA: site_selection
using DataFrames

@info "Loading data package"
here = @__DIR__
ex_domain = ADRIA.load_domain(joinpath(here, "Example_domain"))

@info "Create criteria structure using scenario DataFrame."
criteria_df = ADRIA.sample(ex_domain, 5) # get scenario dataframe

area_to_seed = 1.5 * 10^-6 # area of seeded corals in km^2
nreps = 30 # number of dhw and wave replicates you want to use
ts = 5 # time step to perform site selection at
scen = 5
alg_ind = criteria_df.guided[scen]# MCDA algorithm to use (1-3)

ranks = site_selection(ex_domain, criteria_df, area_to_seed, ts, nreps, scen, alg_ind)
