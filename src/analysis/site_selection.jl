using NamedDims, AxisKeys
using ADRIA: ResultSet, n_locations


"""
    seeded_sites_frequency(rs::ResultSet,scens::NamedTuple, log_type::String)::NamedDimsArray

Count frequency of seeded sites for scenarios satisfying a condition.

# Arguments
- 'rs' : ResultSet
- `scens` : contains scenario ids for scenarios satisfying the condition of interest.
- 'log_type` : "seed", "shade" or "fog" indicating the intervention log to use in calculating frequencies.

# Returns
NamedDimsArray(:locations,:rcps)

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
robust = ADRIA.analysis.find_robust(rs, y, rule_func, [45, 60])

# find site selection frequencies for all robust scenarios
robust_selection_frequencies = intervention_frequency(rs, scen_indices)
"""
function intervention_frequency(rs::ResultSet, scen_indices::NamedTuple, log_type::String)::NamedDimsArray
    occursin(log_type, "seed shade fog") || ValueError("Unsupported log")

    interv_log = getfield(rs, Symbol("$(log_type)_log"))
    # retrieve RCPs
    rcps = Symbol.(keys(scen_indices))
    n_locs = n_locations(rs)

    # create frequencies storage container
    seeded_sites_store = NamedDimsArray(KeyedArray(zeros(n_locs, length(rcps)), (1:n_locs, [rcps[k] for k in eachindex(rcps)])), (:locations, :rcps))

    for rcp in rcps
        # select scenarios satisfying condition and sum up selection tally for each site
        seed_log = dropdims(sum(interv_log[scenarios=scen_indices[rcp]], dims=2), dims=2)
        seeded_sites_store(rcp) .= vec(dropdims(sum(seed_log .> 0, dims=[1, 3]), dims=3))
    end

    return seeded_sites_store
end
