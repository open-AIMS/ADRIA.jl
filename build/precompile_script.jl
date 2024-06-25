using CSV, DataFrames, GeoDataFrames

precompile(CSV.read, (String, DataFrame))
precompile(GeoDataFrames.read, (String,))

# precompile(load_domain, (String, Int64))
# precompile(Domain, (String, Int64, String, String, String, String, String, String, String))
# precompile(EnvLayer, (String, String, String, String, String, String, String))

# let
#     here = @__DIR__
#     ex_dir = joinpath(here, "../examples")
#     @debug "Pre-running examples to reduce future spin-up time"

#     f() = begin
#         @showprogress 1 for _ in 1:10
#         end
#     end
#     b = redirect_stdout(f, devnull);

#     ex_domain = ADRIA.load_domain(joinpath(ex_dir, "Example_domain"), 45)
#     p_df = ADRIA.load_scenarios(ex_domain, joinpath(ex_dir, "example_scenarios.csv"))

#     ENV["ADRIA_THRESHOLD"] = 1e-6
#     ex_domain.sim_constants.tf = 3
#     ds = (raw=nothing, site_ranks=nothing, seed_log=nothing, fog_log=nothing, shade_log=nothing)
#     run_scenario(ex_domain, (1, p_df[1, :]), 1, ds)
#     run_scenario((ex_domain, 1, p_df[end, :]), 1, ds)
#     delete!(ENV, "ADRIA_THRESHOLD")
# end
