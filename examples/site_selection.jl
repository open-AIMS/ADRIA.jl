using ADRIA
using ADRIA: run_site_selection
using DataFrames


@info "Loading data package"
here = @__DIR__
dom = ADRIA.load_domain(joinpath(here, "Example_domain"))

criteria_df = ADRIA.sample(ex_domain, 5) # get scenario dataframe

area_to_seed = 1.5 * 10^-6  # area of seeded corals in km^2
ts = 5  # time step to perform site selection at

sumcover = 0.1 .* ones(5, size(dom.site_data, 1))  # initial coral cover matching number of criteria samples
ranks = run_site_selection(dom, criteria_df[criteria_df.guided.>0, :], sumcover, area_to_seed, ts)
