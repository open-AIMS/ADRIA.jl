# Analysis


## Available Metrics

A range of metrics are defined as part of the ADRIA framework.

See the [metrics](@ref) page for more details.

Here, we extract results for specific metrics for each timestep and sites
for all the scenarios run.

```julia
tac = ADRIA.metrics.total_absolute_cover(rs)
rsv = ADRIA.metrics.relative_shelter_volume(rs)
juves = ADRIA.metrics.juveniles(rs)
```

Some times we're more interested in the scenario-level performance:

```julia
s_tac = ADRIA.metrics.scenario_tac(rs)
s_rsv = ADRIA.metrics.scenario_rsv(rs)
s_juves = ADRIA.metrics.scenario_juveniles(rs)
```

## Visualizations

The Makie.jl ecosystem is used to produce figures.

```julia
# Install additional packages if necessary
# ]add GLMakie GeoMakie GraphMakie

# Import additional packages and the visualization extension will compile.
using GLMakie, GeoMakie, GraphMakie
# [ Info: Precompiling AvizExt [7cc86020-4844-5174-99c6-5a0a5943f024]

# If using VS Code, you may need to deactivate the inline plotting feature
# to make figures appear.
Makie.inline!(false)

using ADRIA
```


## Scenario outcomes

Plots scenario trajectories.

```julia
# First run some scenarios
# dom = ADRIA.load_domain("path to domain data")

# Create some scenarios
# num_samples = 4096
# scens = ADRIA.sample(dom, num_samples)

# Run the scenarios ...
rs = ADRIA.run_scenarios(scens, dom, "45")
# ... or load pre-existing scenarios
# rs = ADRIA.load_results("path to a result set")

# Obtain data to plot (here, scenario-level metrics)
s_tac = ADRIA.metrics.scenario_total_cover(rs)

# Plot a quick scenario overview
ADRIA.viz.scenario(rs, s_tac; axis_opts=Dict(:ylabel=>"Example Metric"))

# Can also compose a figure with subplots
s_juves = ADRIA.metrics.scenario_relative_juveniles(rs)

tf = Figure(resolution=(1600, 600))  # resolution in pixels

# Implicitly create a single figure with 2 columns
ADRIA.viz.scenario!(tf[1, 1], rs, s_tac; opts=Dict(:by_RCP => false), axis_opts=Dict(:title => "TAC [m²]"));
ADRIA.viz.scenario!(tf[1, 2], rs, s_juves; opts=Dict(:by_RCP => false), axis_opts=Dict(:title => "Juveniles [%]"));

tf  # display the figure
save("aviz_scenario.png", tf)  # save the figure to a file
```

![Quick scenario plots](assets/imgs/aviz_scenario.png?raw=true "Quick scenario plots")

# Other visualizations

The examples below are to illustrate usage. For further information on each method of analysis,
see the documentation for the given function.

Some options shared for the plots below are defined here.

```julia
# Some shared options for the example plots below
fig_opts = Dict(:resolution => (1600, 800))

# Factors of Interest
opts = Dict(
    :factors => [:RCP, :dhw_scenario, :wave_scenario, :guided, :seed_TA, :seed_CA, :fogging, :SRM, :a_adapt, :n_adapt]
)
```

## PAWN sensitivity (heatmap overview)

The PAWN sensitivity analysis method is a moment-independent approach to Global Sensitivity Analysis.
It is described as producing robust results at relatively low sample sizes, and is used to screen
factors (i.e., identification of important factors) and rank factors as well (ordering factors by
their relative contribution towards a given quantity of interest).

```julia
# Sensitivity (of mean scenario outcomes to factors)
mean_s_tac = vec(mean(s_tac, dims=1))
tac_Si = ADRIA.sensitivity.pawn(rs, mean_s_tac)
pawn_fig = ADRIA.viz.pawn(
    tac_Si;
    opts,
    fig_opts
)
save("pawn_si.png", pawn_fig)
```

![PAWN sensitivity plots](assets/imgs/pawn_si.png?raw=true "PAWN sensitivity plots")

## Temporal Sensitivity Analysis

Temporal (or Time-varying) Sensitivity Analysis applies sensitivity analysis to model outputs
over time. The relative importance of factors and their influence on outputs over time can 
then be examined through this analysis.

```julia
tsa_s = ADRIA.sensitivity.tsa(rs, s_tac)
tsa_fig = ADRIA.viz.tsa(
    rs,
    tsa_s;
    opts,
    fig_opts
)
save("tsa.png", tsa_fig)
```

![Plots of Temporal Sensitivities](assets/imgs/tsa.png?raw=true "Temporal Sensitivity Analysis")

## Regional Sensitivity Analysis

Regional Sensitivity Analysis is a monte-carlo filtering approach. The aim of RSA is to aid in
identifying which (group of) factors drive model outputs and their active areas of factor space.

```julia
tac_rs = ADRIA.sensitivity.rsa(rs, mean_s_tac; S=10)
rsa_fig = ADRIA.viz.rsa(
    rs,
    tac_rs,
    ["n_adapt", "wave_scenario", "dhw_scenario", "seed_TA", "seed_CA", "fogging", "SRM"];
    opts,
    fig_opts
)
save("rsa.png", rsa_fig)
```

![Plots of Regional Sensitivities](assets/imgs/rsa.png?raw=true "Regional Sensitivity Analysis")

## Outcome mapping

A monte-carlo filtering approach similar to Regional Sensitivity Analysis.

As the name implies, outcome mapping aids in identifying the relationship between model outputs
and the region of factor space that led to those outputs.

```julia
tf = Figure(resolution=(1600, 1200))  # resolution in pixels

# Indicate factor values that are in the top 50 percentile
tac_om_50 = ADRIA.sensitivity.outcome_map(rs, mean_s_tac, x -> any(x .>= 0.5); S=20)
ADRIA.viz.outcome_map!(
    tf[1, 1],
    rs,
    tac_om_50,
    ["n_adapt", "wave_scenario", "dhw_scenario", "seed_TA", "seed_CA", "fogging", "SRM"];
    axis_opts=Dict(:title => "Regions which lead to Top 50th Percentile Outcomes", :ylabel => "TAC [m²]")
)

# Indicate factor values that are in the top 30 percentile
tac_om_70 = ADRIA.sensitivity.outcome_map(rs, mean_s_tac, x -> any(x .>= 0.7); S=20)
ADRIA.viz.outcome_map!(
    tf[2, 1],
    rs,
    tac_om_70,
    ["n_adapt", "wave_scenario", "dhw_scenario", "seed_TA", "seed_CA", "fogging", "SRM"];
    axis_opts=Dict(:title => "Regions which lead to Top 30th Percentile Outcomes", :ylabel => "TAC [m²]"))

save("outcome_map.png", tf)
```

![Outcome mapping](assets/imgs/outcome_map.png?raw=true "Outcome mapping")

## GUI for high-level exploration (prototype only!)

```julia
# To explore results interactively
ADRIA.viz.explore("path to Result Set")

# or, if the result set is already loaded: 
# ADRIA.viz.explore(rs)
```

![Standalone app for data exploration](assets/imgs/aviz_app.png?raw=true "Data Exploration App")

