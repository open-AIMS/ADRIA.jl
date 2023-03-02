# Aviz 

The ADRIA Visualization extension package.

Uses Makie.jl to quickly produce indicative figures.

# Scenario outcomes

```julia
using ADRIA
using Aviz


# First run some scenarios

# dom = ADRIA.load_domain("path to domain data")

# Create some scenarios
# num_samples = 4096
# scens = ADRIA.sample(dom, num_samples)

# Run the scenarios ...
# rs = ADRIA.run_scenarios(scens, dom, "45", remove_workers=true)

# ... or load pre-existing scenarios
rs = ADRIA.load_results("path to a result set")

# Obtain data to plot (here, scenario-level metrics)
s_tac = ADRIA.metrics.scenario_total_cover(rs)

# Plot a quick scenario overview
Aviz.scenario(rs, s_tac; axis_opts=Dict(:ylabel=>"Example Metric"))

# Can also compose a figure with subplots
s_juves = ADRIA.metrics.scenario_relative_juveniles(rs)

tf = Figure(resolution=(1600, 600))  # resolution in pixels

# Implicitly create a single figure with 2 columns
Aviz.scenario!(tf[1, 1], rs, s_tac; opts=Dict(:by_RCP => false), axis_opts=Dict(:title => "TAC [m²]"));
Aviz.scenario!(tf[1, 2], rs, s_juves; opts=Dict(:by_RCP => false), axis_opts=Dict(:title => "Juveniles [%]"));

tf  # display the figure
save("aviz_scenario.png", tf)
```

![Quick scenario plots](assets/imgs/aviz_scenario.png?raw=true "Quick scenario plots")

```julia
# Some shared options for the example plots below
fig_opts = Dict(:resolution => (1600, 800))

# Factors of Interest
opts = Dict(
    :factors => [:RCP, :dhw_scenario, :wave_scenario, :guided, :seed_TA, :seed_CA, :fogging, :SRM, :a_adapt, :n_adapt]
)
```

# PAWN sensitivity (heatmap overview)

```julia
# Sensitivity (of mean scenario outcomes to factors)
mean_s_tac = vec(mean(s_tac, dims=1))
tac_Si = ADRIA.sensitivity.pawn(rs, mean_s_tac)
pawn_fig = Aviz.pawn(
    tac_Si;
    opts,
    fig_opts
)
save("pawn_si.png", pawn_fig)
```

![PAWN sensitivity plots](assets/imgs/pawn_si.png?raw=true "PAWN sensitivity plots")

# Temporal Sensitivity Analysis

```julia
tsa_s = ADRIA.sensitivity.tsa(rs, s_tac)
tsa_fig = Aviz.tsa(
    rs,
    tsa_s;
    opts,
    fig_opts
)
save("tsa.png", tsa_fig)
```

![Plots of Temporal Sensitivities](assets/imgs/tsa.png?raw=true "Temporal Sensitivity Analysis")

# Regional Sensitivity Analysis

```julia
tac_rs = ADRIA.sensitivity.rsa(rs, mean_s_tac; S=10)
rsa_fig = Aviz.rsa(
    rs,
    tac_rs,
    ["n_adapt", "wave_scenario", "dhw_scenario", "seed_TA", "seed_CA", "fogging", "SRM"];
    opts,
    fig_opts
)
save("rsa.png", rsa_fig)
```

![Plots of Regional Sensitivities](assets/imgs/rsa.png?raw=true "Regional Sensitivity Analysis")

# Outcome mapping

```julia
tf = Figure(resolution=(1600, 1200))  # resolution in pixels

# Indicate factor values that are in the top 50 percentile
tac_om_50 = ADRIA.sensitivity.outcome_map(rs, mean_s_tac, x -> any(x .>= 0.5); S=20)
Aviz.outcome_map!(
    tf[1, 1],
    rs,
    tac_om_50,
    ["n_adapt", "wave_scenario", "dhw_scenario", "seed_TA", "seed_CA", "fogging", "SRM"];
    axis_opts=Dict(:title => "Regions which lead to Top 50th Percentile Outcomes", :ylabel => "TAC [m²]")
)

# Indicate factor values that are in the top 30 percentile
tac_om_70 = ADRIA.sensitivity.outcome_map(rs, mean_s_tac, x -> any(x .>= 0.7); S=20)
Aviz.outcome_map!(
    tf[2, 1],
    rs,
    tac_om_70,
    ["n_adapt", "wave_scenario", "dhw_scenario", "seed_TA", "seed_CA", "fogging", "SRM"];
    axis_opts=Dict(:title => "Regions which lead to Top 30th Percentile Outcomes", :ylabel => "TAC [m²]"))

save("outcome_map.png", tf)
```

![Outcome mapping](assets/imgs/outcome_map.png?raw=true "Outcome mapping")

# GUI

A GUI for quick visualization and analysis is also provided.
This can be launched programmatically from the REPL, however, a standalone app will also be made available.

```
using ADRIA, Aviz

# Load some results
rs = ADRIA.load_results("...")

Aviz.explore(rs)
```

![Standalone app for data exploration](assets/imgs/aviz_app.png?raw=true "Data Exploration App")