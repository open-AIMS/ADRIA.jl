"""
    intervention_frequency(rs::ResultSet, scen_indices::NamedTuple, log_type::Symbol)::NamedDimsArray

Count number of times a location of selected for intervention
Count frequency of seeded locations for scenarios satisfying a condition.

# Arguments
- 'rs' : ResultSet
- `scen_indices` : rcp_id => scenario id that satisfy a condition of interest.
- 'log_type` : the intervention log to use in calculating frequencies (one of :seed, :shade or :fog).

# Returns
NamedDimsArray(:locations, :rcps)

# Example
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
robust_scens = ADRIA.analysis.find_robust(rs, y, rule_func, [45, 60])

# Retrieve seeding intervention frequency for robust scenarios
robust_selection_frequencies = ADRIA.analysis.intervention_frequency(rs, robust_scens, :seed)
"""
function intervention_frequency(rs::ResultSet, scen_indices::NamedTuple, log_type::Symbol)::NamedDimsArray
    log_type ∈ [:seed, :shade, :fog] || ValueError("Unsupported log")

    # Get requested log
    interv_log = getfield(rs, Symbol("$(log_type)_log"))
    rcps = collect(Symbol.(keys(scen_indices)))
    n_locs = n_locations(rs)

    interv_freq = NamedDimsArray(zeros(n_locs, length(rcps)), locations=rs.location_ids, rcps=rcps)
    for rcp in rcps
        # Select scenarios satisfying condition and tally selection for each location
        logged_data = dropdims(sum(interv_log[scenarios=scen_indices[rcp]], dims=:coral_id), dims=:coral_id)
        interv_freq(rcps=rcp) .= vec(dropdims(sum(logged_data .> 0, dims=[:timesteps, :scenarios]), dims=:timesteps))
    end

    return interv_freq
end
