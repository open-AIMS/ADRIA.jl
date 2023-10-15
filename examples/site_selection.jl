using ADRIA
using ADRIA:
    run_site_selection,
    ranks_to_frequencies,
    location_selection_frequencies,
    summed_inverse_rank
using DataFrames
using Statistics, StatsBase

@info "Loading data package"
here = @__DIR__
dom = ADRIA.load_domain(joinpath(here, "Moore_2023-08-17"), "45")
scens = ADRIA.sample_site_selection(dom, 8) # Get scenario dataframe.

area_to_seed = 962.11  # Area of seeded corals in m^2.

# Initial coral cover matching number of criteria samples (size = (no. criteria scens, no. of sites)).
sum_cover = repeat(sum(dom.init_coral_cover; dims=1), size(scens, 1))

@info "Run site selection"
# Use run_site_selection to get ranks
ranks = run_site_selection(dom, scens, sum_cover, area_to_seed)

# Get frequencies with which each site is selected for each rank for set of stand alone location selections.
rank_freq = ranks_to_frequencies(ranks[intervention=1])

@info "Calculate rank aggregations"
# Get location selection freqencies for set of standalone location selections.
location_selection_frequency = location_selection_frequencies(ranks[intervention=1])
# Get summed inverse rank for set of standalone location selections.
# Measure of magnitude and frequency of high rank.
inv_summ_rank = summed_inverse_rank(ranks[intervention=1])

# Use aggregation function within run_site_selection to get direct output.
# To get rank frequencies:
rank_frequencies_seed = run_site_selection(
    dom, scens, sum_cover, area_to_seed, ranks_to_frequencies, 1
)
rank_frequencies_seed = run_site_selection(
    dom, scens, sum_cover, area_to_seed, location_selection_frequencies, 1
)
rank_frequencies_seed = run_site_selection(
    dom, scens, sum_cover, area_to_seed, summed_inverse_rank, 1
)

# Example using ADRIA runs
scens = ADRIA.sample(dom, 2^5) # Get scenario dataframe.
rs = ADRIA.run_scenarios(dom, scens, "45") # Run scenarios.

# Get frequencies with which each site is selected for each rank for set of runs.
rank_freq = ranks_to_frequencies(rs.ranks[intervention=1])

# Get selection frequencies for set of runs.
selection_freq = location_selection_frequencies(rs.ranks[intervention=1])

# Get selection frequencies over time for unguided runs only.
unguided_freq = location_selection_frequencies(
    rs.seed_log[scenarios=findall(scens.guided .== 1)]
)

# Get summed inverse rank for set of runs.
# Measure of magnitude and frequency of high rank.
inv_summ_rank = summed_inverse_rank(rs.ranks[intervention=1])
# Get summed inverse rank over time.
inv_summ_rank = summed_inverse_rank(rs.ranks[intervention=1]; dims=[:scenarios])
