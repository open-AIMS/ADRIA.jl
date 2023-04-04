using ADRIA
using ADRIA: run_site_selection
using DataFrames


using ADRIA

@info "Loading data package"
here = @__DIR__
dom = ADRIA.load_domain(joinpath(here, "Example_domain/"), "45")

@info "Sample criteria weights"
criteria_df = ADRIA.sample_site_selection(dom, 8) # get scenario dataframe

area_to_seed = 962.11  # area of seeded corals in m^2
ts = 5  # time step to perform site selection at

# define functions for tolerances
f_coral_cover(param) = area_to_seed * param

@info "Perform site selection"
# initial coral cover matching number of criteria samples (size = (no. criteria scens, no. of sites))
tolerances = (iv__coral_cover=(>, x -> f_coral_cover(x)),
    iv__heat_stress=(<, x -> x),
    iv__wave_stress=(<, x -> x))
ranks = run_site_selection(dom, criteria_df, tolerances, ts)
