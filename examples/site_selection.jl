using ADRIA
using ADRIA: run_site_selection, ranks_to_frequencies, ranks_to_frequencies_ts, ranks_to_location_order, location_selection_frequencies
using DataFrames
using Statistics, StatsBase

@info "Loading data package"
here = @__DIR__
dom = ADRIA.load_domain(joinpath(here, "Moore_2023-08-17"), "45")
scens = ADRIA.sample_site_selection(dom, 8) # Get scenario dataframe.

scens = ADRIA.sample_site_selection(dom, 8) # Get scenario dataframe.

area_to_seed = 962.11  # Area of seeded corals in m^2.

# Initial coral cover matching number of criteria samples (size = (no. criteria scens, no. of sites)).
sum_cover = repeat(sum(dom.init_coral_cover, dims=1), size(scens, 1))

@info "Run site selection"
# Use run_site_selection to get ranks
ranks = run_site_selection(dom, criteria_df, sum_cover, area_to_seed)
# Use an aggregation function to get location selection frequency.

@info "Calculate rank aggregations"
n_loc_int = 5 # number of sites selected at each step.
location_selection_frequency = location_selection_frequencies(rank, "seed", n_loc_int)

# Use aggregation function within run_site_selection to get direct output.
# To get rank frequencies:
rank_frequencies_seed = run_site_selection(dom, scens, sum_cover, area_to_seed, ranks_to_frequencies, "seed")
# To get location rank order as site ids:
location_order_seed = run_site_selection(dom, scens, sum_cover, area_to_seed, ranks_to_location_order, "seed")

# Example using ADRIA runs
scens = ADRIA.sample(dom, 8) # Get scenario dataframe.
rs = ADRIA.run_scenarios(scens, dom, "45") # Run scenarios.
rank_freq = ranks_to_frequencies(rs.ranks, "seed") # Get rank frequencies.
rank_freq_ts = ranks_to_frequencies_ts(rs.ranks, "seed") # Get rank frequencies over time.
