# Scenario Discovery


```julia
using ADRIA

dom = ADRIA.load_domain("...", "<RCP>")
scens = ADRIA.sample(dom, 4096)

rs = ADRIA.run_scenarios(dom, scens, "45")
# rs = ADRIA.load_results("...")

# Calculate representative statistic for all metrics of interest
# Here, total cover, shelter volume and juvenile population
tac = ADRIA.metrics.scenario_total_cover(rs)
mean_tac = vec(mean(tac, dims=1))

rsv = ADRIA.metrics.scenario_rsv(rs)
mean_sv = vec(mean(rsv, dims=1))

r_juves = ADRIA.metrics.scenario_relative_juveniles(rs)
mean_juves = vec(mean(r_juves, dims=1))

# Create matrix of all metrics
y = hcat(mean_tac, mean_sv, mean_juves)

# Define "robust scenario" as one where all metric outcomes >= 30th percentile.
rule_func = x -> all(x .>= 0.3)

# Identify robust scenarios for a specific RCP (4.5)
robust = ADRIA.analysis.find_robust(rs, y, rule_func, [45])

# Output robust scenario IDs
@info robust.RCP45

# Could qualitatively examine inputs that led to robust scenarios ...
rs.inputs[robust.RCP45, :]

# ... or mark behavioural and non-behavioural outcomes for further analysis
# e.g., with random forest.
behave = zeros(size(y, 1))
behave[robust.RCP45] .= 1.0

# Next step is to test these for robustness across environmental conditions ...
# [TODO]
```
