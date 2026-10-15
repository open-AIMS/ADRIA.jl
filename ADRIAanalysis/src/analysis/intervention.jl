"""
    intervention_frequency(rs::ResultSet, scen_indices::NamedTuple, log_type::Symbol)::YAXArray

Count number of times a location of selected for intervention
Count frequency of coral aquaculture sites for scenarios satisfying a condition.

# Arguments
- 'rs' : ResultSet
- `scen_indices` : rcp_id => scenario id that satisfy a condition of interest.
- 'log_type` : the intervention log to use in calculating frequencies (one of :CAq, :LvM or :Shd).
  `:Shd` uses the combined fog+shade log (`rs.iv_Shd_log`) — fog and shade are not tracked
  as separate logs.

# Returns
YAXArray(:locations, :rcps)

# Examples
```julia
using ADRIA, Statistics

rs = ADRIA.load_results("some result set")

tac = ADRIA.metrics.scenario_total_cover(rs)
rsv = ADRIA.metrics.scenario_rsv(rs)

# Create matrix of mean scenario outcomes
mean_tac = vec(mean(tac, dims=1))
mean_sv = vec(mean(rsv, dims=1))
y = hcat(mean_tac, mean_sv)

# Find all pareto optimal scenarios where all metrics >= 0.9
rule_func = x -> all(x .>= 0.9)
robust_scens = find_robust(rs, y, rule_func, [45, 60])

# Retrieve coral aquaculture intervention frequency for robust scenarios
robust_selection_frequencies = intervention_frequency(rs, robust_scens, :CAq)
"""
function intervention_frequency(
    rs::ResultSet, scen_indices::NamedTuple, log_type::Symbol
)::YAXArray
    log_type ∈ [:CAq, :LvM, :Shd] ||
        throw(ArgumentError("Unsupported log type: $log_type"))

    # CAq/LvM logs carry a `coral_id` axis (deployment is per functional group); Shd is
    # the combined fog+shade log and instead carries an `intervention` axis (fog and
    # shade summed into one signal, consistent with how ResultSet stores them as a
    # single `iv_Shd_log`).
    interv_log = getfield(rs, Symbol("iv_$(log_type)_log"))
    reduce_dim = log_type == :Shd ? :intervention : :coral_id

    rcps = collect(Symbol.(keys(scen_indices)))

    interv_freq = ZeroDataCube(; T=Float64, locations=rs.loc_ids, rcps=rcps)
    for rcp in rcps
        # Select scenarios satisfying condition and tally selection for each location
        logged_data = dropdims(
            sum(interv_log[scenarios = scen_indices[rcp]]; dims=reduce_dim); dims=reduce_dim
        )
        interv_freq[rcps = At(rcp)] .= vec(
            dropdims(sum(logged_data .> 0; dims=(:timesteps, :scenarios)); dims=:timesteps)
        )
    end

    return interv_freq
end
