# Aviz 

The ADRIA Visualization extension package.

Uses Makie.jl to quickly produce indicative figures.

```
using ADRIA
using Aviz

# Load some results
rs = ADRIA.load_results("...")

Aviz.scenario(rs, ADRIA.metrics.scenario_total_cover)

tf = Figure(resolution=(1600, 600))  # resolution in pixels
Aviz.plot.scenario!(tf[1, 1], rs, ADRIA.metrics.scenario_total_cover; opts=Dict(:by_RCP => false), axis_opts=Dict(:title => "TAC"));
Aviz.plot.scenario!(tf[1, 2], rs, ADRIA.metrics.scenario_juveniles; opts=Dict(:by_RCP => false), axis_opts=Dict(:title => "Juveniles"));

tf  # show figure
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