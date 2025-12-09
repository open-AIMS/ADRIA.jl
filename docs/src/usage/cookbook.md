# Cookbook examples

## Extracting model details

Example showcasing how to extract model details, such as:

- the model parameter table
- model specification

and more specific information/data from the above.

```julia
using DataFrames
using ADRIA


# Loading a dataset for a study area (a "domain")
data_pkg = "./Example_domain"
dom = ADRIA.load_domain(data_pkg, 45)

# Get current parameter table (fieldnames and their values)
param_df = ADRIA.param_table(dom)

# Get model specification with lower/upper bounds separated
model_spec = ADRIA.model_spec(dom)

# Export model specification to CSV
ADRIA.model_spec(dom, "model_spec.csv")


# Get parameter details

## Parameter names
p_names = dom.model[:fieldname]

## Current values
p_vals = dom.model[:val]

## ADRIA parameter types
p_types = dom.model[:ptype]

## Parameter bounds (for e.g., to pass into a sampler or optimizer)
## Note: ADRIA integer parameter bounds are set such that ℓ ≤ x ≤ u+1,
## where ℓ is the lower bound and u is the upper bound.
## This is because `floor(x)` is assigned with `update_params!()`.
## Instances where ℓ := x := u indicate uncertain parameters that
## are nevertheless assumed to be constant.
p_bounds = dom.model[:bounds]

## Component groups
p_groups = dom.model[:component]

## All of above as a DataFrame
model_spec = DataFrame(dom.model)


# Get DataFrame of parameter information for a specific sub-component (Intervention, Criteria, Coral)
ADRIA.component_params(dom.model, Intervention)
```

## Generating and running scenarios

```julia
using ADRIA


# Loading data package
dom = ADRIA.load_domain("Example_domain", "<RCP>")

# Creating 128 scenarios based on parameter bounds using the Sobol' method
scens = ADRIA.sample(dom, 128)

# Can also use other samplers
# using Surrogates.QuasiMonteCarlo
# scens = ADRIA.sample(dom, 100, LatinHypercubeSample())

# Can also sample counterfactuals (scenarios with no interventions)
# or scenarios with guided interventions only
# s = ADRIA.sample_cf(dom, 32)
# s = ADRIA.sample_guided(dom, 32)

# Can also load previously generated scenarios
# p_df = ADRIA.load_scenarios(dom, joinpath(here, "example_scenarios.csv"))

# Batch run scenarios. Returns a ResultSet.
# Setting up and running scenarios
rs = ADRIA.run_scenarios(dom, p_df, "45")

# Multiple RCPs can be specified, so long as RCP-specific data is available.
# rs = ADRIA.run_scenarios(dom, p_df, ["45", "60"])

# Single scenario run (returns NamedTuple of results for a single environmental/intervention scenario).
# See documentation for more detail.
# scenario_id = 1
# result = ADRIA.run_scenario(domain::Domain, scenario_id, param_df::DataFrameRow)

# switch_RCPs!(domain, "45")
# res1 = ADRIA.run_scenario(domain, scens[1, :])
# res2 = ADRIA.run_scenario(domain, scens[2, :])
# res3 = ADRIA.run_scenario(domain, scens[3, :], "60")  # run for a different RCP

# The location of the outputs stored on disk
@info ADRIA.store_name(rs)
# "Example_domain__RCPs45__2022-10-19_12_01_26_965"

@info ADRIA.result_location(rs)
# "[some location]/Example_domain__RCPs45__2022-10-19_12_01_26_965"

# Can also load results using a path to the stored result set.
# rs = ADRIA.load_results("path to result set")

# Specific metrics found in the `metrics` submodule.
# tac = ADRIA.metrics.total_absolute_cover(rs)
```

## Intervention location selection

```julia
using ADRIA
using ADRIA: rank_locations


dom = ADRIA.load_domain("path to domain", "45")
scens = ADRIA.sample_site_selection(dom, 8)

# Area of seeded corals in m^2
area_to_seed = 962.11

# Initial coral cover matching number of criteria samples (size = (no. criteria scens, no. of sites)).
sum_cover = repeat(sum(dom.init_coral_cover; dims=1), size(scens, 1))

# Use rank_locations to get ranks
ranks = rank_locations(dom, scens, sum_cover, area_to_seed)
```

## Intervention location selection - summary functions

```julia
using ADRIA
using ADRIA:
    rank_locations,
    ranks_to_frequencies,
    location_selection_frequencies,
    selection_score
using DataFrames
using Statistics, StatsBase

# Load data package
dom = ADRIA.load_domain("path to Domain files", "RCP")

# Select locations for interventions without any model runs
scens = ADRIA.sample_site_selection(dom, 8)  # Get site selection scenario dataframe.

# Area of seeded corals in m^2
area_to_seed = 962.11

# Initial coral cover matching number of criteria samples
sum_cover = repeat(sum(dom.init_coral_cover; dims=1), size(scens, 1))

# Use rank_locations to get ranks
ranks = rank_locations(dom, scens, sum_cover, area_to_seed)

# Get frequencies with which each site is selected for each rank for set of stand alone
# location selections
rank_freq = ranks_to_frequencies(ranks[intervention=1])

# Calculate rank aggregations
# Get location selection freqencies for set of standalone location selections
location_selection_frequency = location_selection_frequencies(ranks[intervention=1])

# Get summed inverse rank for set of standalone location selections
# Measure of magnitude and frequency of high rank
sel_score = selection_score(ranks[intervention=1])

# Use aggregation function within rank_locations to get direct output
# To get rank frequencies:
rank_frequencies_seed = rank_locations(
    dom, scens, sum_cover, area_to_seed, ranks_to_frequencies, 1
)
rank_frequencies_seed = rank_locations(
    dom, scens, sum_cover, area_to_seed, location_selection_frequencies, 1
)
rank_frequencies_seed = rank_locations(
    dom, scens, sum_cover, area_to_seed, selection_score, 1
)

# Example using ADRIA runs
scens = ADRIA.sample(dom, 8)
rs = ADRIA.run_scenarios(dom, scens, "45")

# Get frequencies with which each site was selected for each rank for a set of runs
rank_freq = ranks_to_frequencies(rs.ranks[intervention=1])  # with timesteps not aggregated

# With timesteps aggregated
rank_freq = ranks_to_frequencies(
    rs.ranks[intervention=1];
    agg_func=x -> dropdims(sum(x; dims=:timesteps); dims=:timesteps),
)

# Get selection frequencies for set of runs
selection_freq = location_selection_frequencies(rs.ranks[intervention=1])

# Get selection frequencies over time for unguided runs only
unguided_freq = location_selection_frequencies(
    rs.seed_log[scenarios=findall(scens.guided .>= 1)]
)

# Get selection score for set of runs
# Measure of magnitude and frequency of high rank
sel_score = selection_score(rs.ranks[intervention=1])

# Get selection score for locations over time
sel_score = selection_score(rs.ranks[intervention=1]; dims=[:scenarios])
```
