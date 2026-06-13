using ADRIA

example_dir = joinpath(@__DIR__, "..", "test", "data", "Test_domain")
@assert isdir(example_dir) "Precompile fixture not found: $example_dir"

ENV["ADRIA_DEBUG"] = "false"
dom = ADRIA.load_domain(example_dir, "45")
p_df = ADRIA.sample(dom, 32)
rs = ADRIA.run_scenarios(dom, p_df, "45")
ADRIA.metrics.scenario_total_cover(rs)
delete!(ENV, "ADRIA_DEBUG")
