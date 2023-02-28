using ADRIA
using ADRIA: run_site_selection
using DataFrames


@info "Loading data package"
here = @__DIR__
dom = ADRIA.load_domain(joinpath(here, "Example_domain"))

criteria_df = ADRIA.sample_site_selection(dom, 8) # get scenario dataframe

area_to_seed = 1.5 * 10^-6  # area of seeded corals in km^2
ts = 5  # time step to perform site selection at

# initial coral cover matching number of criteria samples (size = (no. criteria scens, no. of sites))
sum_cover = repeat(sum(dom.init_coral_cover, dims=1), size(criteria_df, 1))

ranks = run_site_selection(dom, criteria_df, sum_cover, area_to_seed, ts)
