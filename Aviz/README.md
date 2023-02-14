# Aviz 

The ADRIA Visualization extension package.

Uses Makie.jl to quickly produce indicative figures.

```
using ADRIA
using Aviz

# Load some results
rs = ADRIA.load_results("...")

# Obtain data to plot (scenario-level metrics)
s_tac = ADRIA.metrics.scenario_total_cover(rs)

# Infers figure x/y labels from metric and ResultSet
# specific axis options
Aviz.plot.scenario(rs, s_tac; axis_opts=Dict(:ylabel=>"Example Metric"))

# Can also compose subplots
s_tac = ADRIA.metrics.scenario_total_cover(rs)
s_juves = ADRIA.metrics.scenario_relative_juveniles(rs)

# Compose figure
tf = Figure(resolution=(1600, 600))  # resolution in pixels
Aviz.plot.scenario!(tf[1, 1], rs, s_tac; opts=Dict(:by_RCP => false), axis_opts=Dict(:title => "TAC"));
Aviz.plot.scenario!(tf[1, 2], rs, s_juves; opts=Dict(:by_RCP => false), axis_opts=Dict(:title => "Juveniles"));

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