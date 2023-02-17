using NamedDims
using ADRIA: ResultSet, n_locations


"""
    seeded_sites_frequency(rs::ResultSet,scens::NamedTuple)::NamedTuple

Count frequency of seeded sites for scenarios satisfying a condition.

# Arguments
- 'rs' : ResultSet
- `scens` : contains scenario ids for scenarios satisfying the condition of interest.

# Returns
NamedDimsArray, where each entry relates to an RCP of interest, e.g., 
[
    RCP45=[... frequency of selection for each site ...]; 
    RCP60=[ ... frequency of selection for each site ...]
]

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
function intervention_frequency(rs::ResultSet, scen_indices::NamedTuple)::NamedDimsArray

    # retrieve RCPs
    rcps = Symbol.(keys(scen_indices))

    # create frequencies storage container
    seeded_sites_store = NamedDimsArray(zeros(n_locations(rs), length(rcps)), (:locations, :rcps))
    for (indx, rcp) in enumerate(rcps)
        # select scenarios satisfying condition and sum up selection tally for each site
        seed_log = dropdims(sum(rs.seed_log[:, :, :, scen_indices[rcp]], dims=2), dims=2)
        seeded_sites_store[rcps=indx] .= vec(dropdims(sum(seed_log .> 0, dims=[1, 3]), dims=3))
    end

    return seeded_sites_store
end
