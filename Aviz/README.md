# Aviz 

The ADRIA Visualization extension package.

Uses Makie.jl to quickly produce indicative figures.

```
using ADRIA
using Aviz

# Load some results
rs = ADRIA.load_results("...")

# Scenario-level metrics

# Quick display of scenario results.
# Infers figure x/y labels from metric and ResultSet
Aviz.scenario(rs, ADRIA.metrics.scenario_total_cover)

# As above, but requires user to provide data and 
# specific axis options
s_tac = ADRIA.metrics.scenario_total_cover(rs)
Aviz.plot.scenario(rs, s_tac; axis_opts=Dict(:ylabel=>"Example Metric"))

# Can also compose subplots
tf = Figure(resolution=(1600, 600))  # resolution in pixels
Aviz.plot.scenario!(tf[1, 1], rs, ADRIA.metrics.scenario_total_cover; opts=Dict(:by_RCP => false), axis_opts=Dict(:title => "TAC"));
Aviz.plot.scenario!(tf[1, 2], rs, ADRIA.metrics.scenario_juveniles; opts=Dict(:by_RCP => false), axis_opts=Dict(:title => "Juveniles"));

tf  # show figure with subplots

# Sensitivity (to mean of scenario outcomes)
tac_Si = ADRIA.sensitivity.pawn(rs.inputs, vec(mean(s_tac, dims=1)))
Aviz.plot.pawn(tac_Si)
```

![Quick scenario plots](assets/imgs/aviz_scenario.png?raw=true "Quick scenario plots")


A GUI for quick visualization and analysis is also provided.
This can be launched programmatically, however, a standalone app will also be made available.

```
using ADRIA
using Aviz

# Load some results
rs = ADRIA.load_results("...")

Aviz.explore(rs)
```

![Standalone app for data exploration](assets/imgs/aviz_app.png?raw=true "Data Exploration App")