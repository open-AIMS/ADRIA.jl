using ADRIA
using ADRIA: run_site_selection
using DataFrames


@info "Loading data package"
here = @__DIR__
dom = ADRIA.load_domain(joinpath(here, "Example_domain"))

criteria_df = ADRIA.sample(dom, 5) # get scenario dataframe

area_to_seed = 1.5 * 10^-6  # area of seeded corals in km^2
ts = 5  # time step to perform site selection at

# initial coral cover matching number of criteria samples (size = (no. criteria scens, no. of sites))
sum_cover = fill(0.1, nrow(criteria_df), nrow(dom.site_data))
ranks = run_site_selection(dom, criteria_df[criteria_df.guided.>0, :], sum_cover[criteria_df.guided.>0, :], area_to_seed, ts)
