using TOML, Statistics, CSV, DataFrames
using ADRIA


ADRIA.setup()  # Load config and set up multiprocessing

@info "Loading data package"
data_pkg = "C:/development/ADRIA_brick/Brick"
brick = ADRIA.load_domain(data_pkg, 45)

@info "Loading example scenarios"
p_df = CSV.read("./example_scenarios.csv", DataFrame, comment="#")


# Batch run scenarios
@info "Setting up and running scenarios"
brick = ADRIA.run_scenarios(p_df, brick)

# Single scenario run (returns NamedTuple of results for a single environmental/intervention scenario)
# result = ADRIA.run_scenario(param_df::DataFrameRow, domain::Domain; rep_id=1)

@info "Reloading results and saving figure"
res = ADRIA.load_results(brick)

# Specific metrics found in the `metrics` submodule.
# Y_TC = ADRIA.metrics.coral_cover(res).total_cover
# Y_o = ADRIA.metrics.summarize_total_cover(res)

# Indicative result display for example only. This function to be removed.
ADRIA._indicative_result_display(res)
