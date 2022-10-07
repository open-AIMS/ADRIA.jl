@doc """
This example file assumes Julia is started in the examples folder.

The accompanying `config.toml` file specifies how many cores to use, 
and the output location to store results in.
"""

using ADRIA


@info "Loading data package"
here = @__DIR__
ex_domain = ADRIA.load_domain(joinpath(here, "Example_domain"))

@info "Loading example scenarios"
p_df = ADRIA.load_scenarios(ex_domain, joinpath(here, "example_scenarios.csv"))

# Batch run scenarios. Returns an updated domain object with the run ID used to gather results later.
@info "Setting up and running scenarios"
ex_domain = ADRIA.run_scenarios(p_df, ex_domain, "45")

# Multiple RCPs can be specified, so long as the data is available.
# ex_domain = ADRIA.run_scenarios(p_df, ex_domain, ["45", "60"])

# Single scenario run (returns NamedTuple of results for a single environmental/intervention scenario).
# See documentation for more detail.
# result = ADRIA.run_scenario(param_df::DataFrameRow, domain::Domain; rep_id=1)::NamedTuple

@info "Reloading results and saving figure"
rs = ADRIA.load_results(ex_domain)

# Can also load results using a string
# rs = ADRIA.load_results("path to result set")

# Specific metrics found in the `metrics` submodule.
# tac = ADRIA.metrics.total_absolute_cover(rs)
