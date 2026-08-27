# # Scenario Discovery
#
# Scenario discovery identifies which scenarios meet a target outcome criterion,
# helping to understand what conditions or interventions lead to desired results.
#
# ## Setup
#
# ```julia
# using ADRIA, ADRIAanalysis, Statistics
# ```
#
# ### ResultSet
#
# The examples below assume an `ADRIAResultSet` `rs`:
#
# ```julia
# dom = ADRIA.load_domain("path to domain data", "<RCP>")
# scens = ADRIA.sample(dom, 4096)
# rs = ADRIA.run_scenarios(dom, scens, "45")
# ```
#
# See [Loading a Domain](@ref), [Generating scenarios](@ref) and
# [Running scenarios](@ref) for more detail.
#
# ## Screening scenarios
#
# `screen_scenarios` identifies scenarios where all outcomes meet a given rule after
# column-wise normalization. A scenario is selected only when the rule holds for
# every outcome column.
#
# ```julia
# tac = ADRIA.metrics.scenario_total_cover(rs)
# mean_tac = vec(mean(tac, dims=1))
#
# rsv = ADRIA.metrics.scenario_rsv(rs)
# mean_sv = vec(mean(rsv, dims=1))
#
# r_juves = ADRIA.metrics.scenario_relative_juveniles(rs)
# mean_juves = vec(mean(r_juves, dims=1))
#
# # Stack outcomes column-wise (one column per metric)
# y = hcat(mean_tac, mean_sv, mean_juves)
#
# # Select scenarios where all normalized outcomes are at or above the 30th percentile
# robust = screen_scenarios(y, x -> x >= 0.3)
# ```
#
# The returned `robust` is a vector of scenario indices. Use it to inspect the inputs
# that produced those outcomes:
#
# ```julia
# rs.inputs[robust, :]
# ```
#
# Or create a binary classification vector for further analysis (e.g., random forests):
#
# ```julia
# behave = zeros(size(y, 1))
# behave[robust] .= 1.0
# ```
