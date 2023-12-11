@doc """
This example file assumes Julia is started in the examples folder.

The accompanying `config.toml` file specifies how many cores to use,
and the output location to store results in.
"""

using ADRIA


@info "Loading data package"
here = @__DIR__
ex_domain = ADRIA.load_domain(joinpath(here, "Test_domain"))

@info "Creating 128 scenarios based on parameter bounds using the Sobol' method"
scens = ADRIA.sample(ex_domain, 128)

# Can also use other samplers
# using Surrogates.QuasiMonteCarlo
# s = ADRIA.sample(ex_domain, 100, LatinHypercubeSample())

# Can also sample counterfactuals (scenarios with no interventions)
# or scenarios with guided interventions only
# s = ADRIA.sample_cf(ex_domain, 20)
# s = ADRIA.sample_guided(ex_domain, 20)

# Batch run scenarios. Returns a ResultSet.
@info "Setting up and running scenarios for RCP 4.5"
rs = ADRIA.run_scenarios(ex_domain, scens, "45")

# Multiple RCPs can be specified, so long as the data is available.
# rs = ADRIA.run_scenarios(ex_domain, p_df, ["45", "60"])

# Single scenario run (returns NamedTuple of results for a single environmental/intervention scenario).
# See documentation for more detail.
# domain = ADRIA.switch_RCPs!(domain, "45")
# res1 = ADRIA.run_scenario(domain, scens[1, :])
# res2 = ADRIA.run_scenario(domain, scens[2, :])
# res3 = ADRIA.run_scenario(domain, scens[3, :], "60")  # run for a different RCP


# Name of result store
@info ADRIA.store_name(rs)
# "Example_domain__RCPs45__2022-10-19_12_01_26_965"

# The location of the outputs stored on disk
@info ADRIA.result_location(rs)
# "[some location]/Example_domain__RCPs45__2022-10-19_12_01_26_965"

# Can also load results using a path to the stored result set.
# rs = ADRIA.load_results("path to result set")

# Specific metrics found in the `metrics` submodule.
# tac = ADRIA.metrics.total_absolute_cover(rs)
