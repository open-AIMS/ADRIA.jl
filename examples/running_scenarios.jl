@doc """
This example file assumes Julia is started in the examples folder.

The accompanying `config.toml` file specifies how many cores to use,
and the output location to store results in.
"""

using ADRIA


@info "Loading data package"
example_domain_path = joinpath(pkgdir(ADRIA), "examples", "Test_domain")
dom = ADRIA.load_domain(example_domain_path)

@info "Loading example scenarios"
p_df = ADRIA.load_scenarios(dom, joinpath(here, "example_scenarios.csv"))

# Batch run scenarios. Returns a ResultSet.
@info "Setting up and running scenarios"
rs = ADRIA.run_scenarios(dom, p_df, "45")

# Multiple RCPs can be specified, so long as the data is available.
# rs = ADRIA.run_scenarios(dom, p_df, ["45", "60"])

# Single scenario run (returns NamedTuple of results for a single environmental/intervention scenario).
# See documentation for more detail.
# scenario_id = 1
# result = ADRIA.run_scenario(domain::Domain, scenario_id::Int64, param_df::DataFrameRow)::NamedTuple


# The location of the outputs stored on disk
@info ADRIA.store_name(rs)
# "Example_domain__RCPs45__2022-10-19_12_01_26_965"

@info ADRIA.result_location(rs)
# "[some location]/Example_domain__RCPs45__2022-10-19_12_01_26_965"

# Can also load results using a path to the stored result set.
# rs = ADRIA.load_results("path to result set")

# Specific metrics found in the `metrics` submodule.
# tac = ADRIA.metrics.total_absolute_cover(rs)
