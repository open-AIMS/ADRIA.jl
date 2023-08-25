# Analysis

This section presents tools for analysing model generate data, including functions to
extract metrics and plot graphs.

## Setup
### Makie

The Makie.jl ecosystem is used to produce figures.

Install additional packages if necessary

```julia
]add GLMakie GeoMakie GraphMakie
```

Import additional packages and the visualization extension will compile.

```julia
using GLMakie, GeoMakie, GraphMakie
using ADRIA
```

If using VS Code, you may need to deactivate the inline plotting feature to make figures
appear.

```julia
Makie.inline!(false)
```

### Result Set

All metrics and visualization tools presented here can be used with data generated from
ADRIA's model or from any other model. Following, we show usage examples considering ADRIA
result set `rs`:

```julia
# Load domain data
dom = ADRIA.load_domain("path to domain data")

# Create some scenarios
num_samples = 4096
scens = ADRIA.sample(dom, num_samples)

# Run the model for generated scenarios
rcp_45 = "45"
rs = ADRIA.run_scenarios(scens, dom, rcp_45)
```

See the previous sections [Loading a Domain](@ref), [Generating scenarios](@ref) and
[Running scenarios](@ref) for more information.

## Metrics Extraction

A range of metrics are defined as part of the ADRIA framework. See the [Metrics](@ref)
page for more details.

Here, we extract results for specific metrics for each timestep and sites for all the
scenarios run. The result of each line above is a 3-dimensional Array of timesteps, sites
and scenarios:

```julia
tac = ADRIA.metrics.total_absolute_cover(rs)
rsv = ADRIA.metrics.relative_shelter_volume(rs)
juves = ADRIA.metrics.relative_juveniles(rs)
```

We can also look at scenario-level metrics. They aggregate the above metrics across the
`site` dimension. The result is a 2-dimensional Array of timesteps and scenarios:

```julia
s_tac = ADRIA.metrics.scenario_total_cover(rs)
s_rsv = ADRIA.metrics.scenario_rsv(rs)
s_juves = ADRIA.metrics.scenario_relative_juveniles(rs)
```

## Visualization

The examples below are to illustrate usage. For further information on each method of
analysis, see the documentation for the given function.

Some options shared for the plots below are defined here.

```julia
# Some shared options for the example plots below
fig_opts = Dict(:resolution => (1600, 800))

# Factors of Interest
opts = Dict(
    :factors => [
        :RCP,
        :dhw_scenario,
        :wave_scenario,
        :guided,
        :N_seed_TA,
        :N_seed_CA,
        :fogging,
        :SRM,
        :a_adapt
    ]
)
```

### Scenario outcomes

One can plot a quick scenario overview:

```julia
s_tac = ADRIA.metrics.scenario_total_cover(rs)
ADRIA.viz.scenario(rs, s_tac; axis_opts=Dict(:ylabel=>"Example Metric"))
```

And compose a figure with subplots:

```julia
s_tac = ADRIA.metrics.scenario_total_cover(rs)
s_juves = ADRIA.metrics.scenario_relative_juveniles(rs)

tf = Figure(resolution=(1600, 600))  # resolution in pixels

# Implicitly create a single figure with 2 columns
ADRIA.viz.scenario!(tf[1, 1], rs, s_tac; opts=Dict(:by_RCP => false, :legend=>false), axis_opts=Dict(:title => "TAC [m²]"));
ADRIA.viz.scenario!(tf[1, 2], rs, s_juves; axis_opts=Dict(:title => "Juveniles [%]"));

tf  # display the figure
save("aviz_scenario.png", tf)  # save the figure to a file
```

![Quick scenario plots](/ADRIA.jl/dev/assets/imgs/aviz_scenario.png?raw=true "Quick scenario plots")


### PAWN sensitivity (heatmap overview)

The PAWN sensitivity analysis method is a moment-independent approach to Global Sensitivity
Analysis. It is described as producing robust results at relatively low sample sizes, and
is used to screen factors (i.e., identification of important factors) and rank factors as
well (ordering factors by their relative contribution towards a given quantity of interest).

```julia
using Statistics

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

![PAWN sensitivity plots](/ADRIA.jl/dev/assets/imgs/pawn_si.png?raw=true "PAWN sensitivity plots")

### Temporal Sensitivity Analysis

Temporal (or Time-varying) Sensitivity Analysis applies sensitivity analysis to model
outputs over time. The relative importance of factors and their influence on outputs over
time can then be examined through this analysis.

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

![Plots of Temporal Sensitivities](/ADRIA.jl/dev/assets/imgs/tsa.png?raw=true "Temporal Sensitivity Analysis")


### Time Series Clustering

The Time Series Clustering algorithm clusters together series (typically time series)
with similar behavior. This is achieved by computing the Euclidian distance between each
pair of series weighted by a correlation factor that takes into account the quotient
between their complexities.

```julia
# Extract metric from scenarios
s_tac = ADRIA.metrics.scenario_total_cover(rs)

# Cluster scenarios
n_clusters = 6
clusters = ADRIA.analysis.time_series_clustering(s_tac, n_clusters)

axis_opts = Dict(
    :title => "Time Series Clustering with $n_clusters clusters",
    :ylabel => "TAC [m²]",
    :xlabel => "Timesteps [years]"
)
tsc_fig = ADRIA.viz.ts_cluster(
    s_tac,
    clusters;
    fig_opts=fig_opts,
    axis_opts=axis_opts
)

# Save final figure
save("tsc.png", tsc_fig)
```

![Plots of Temporal Sensitivities](/ADRIA.jl/dev/assets/imgs/tsc.png?raw=true "Time Series Cluster")

### Time Series Clustering Map

When using Time Series Clustering to cluster among multiple locations using some metric, it
is possible to visualize the result as a map.

```julia
using Statistics

# Extract metric from scenarios
tac = ADRIA.metrics.total_absolute_cover(rs)

# Summarize scenarios for each site
m_tac = ADRIA.metrics.loc_trajectory(median, tac)

# Cluster scenarios
n_clusters = 6
loc_clusters = ADRIA.analysis.time_series_clustering(rs, m_tac, n_clusters)

# Plot figure
tsc_map_fig = ADRIA.viz.map(rs, m_tac, loc_clusters)

# Save final figure
save("tsc_map.png", tsc_map_fig)
```

![Plots of Temporal Sensitivities](/ADRIA.jl/dev/assets/imgs/tsc_map.png?raw=true "Spatial Time Series Cluster")

### Regional Sensitivity Analysis

Regional Sensitivity Analysis is a monte-carlo filtering approach. The aim of RSA is to aid
in identifying which (group of) factors drive model outputs and their active areas of
factor space.

```julia
tac_rs = ADRIA.sensitivity.rsa(rs, mean_s_tac; S=10)
rsa_fig = ADRIA.viz.rsa(
    rs,
    tac_rs,
    ["dhw_scenario", "wave_scenario", "N_seed_TA", "N_seed_CA", "fogging", "SRM"];
    opts,
    fig_opts
)
save("rsa.png", rsa_fig)
```

![Plots of Regional Sensitivities](/ADRIA.jl/dev/assets/imgs/rsa.png?raw=true "Regional Sensitivity Analysis")

### Outcome mapping

A monte-carlo filtering approach similar to Regional Sensitivity Analysis.

As the name implies, outcome mapping aids in identifying the relationship between model
outputs and the region of factor space that led to those outputs.

```julia
tf = Figure(resolution=(1600, 1200))  # resolution in pixels

# Indicate factor values that are in the top 50 percentile
tac_om_50 = ADRIA.sensitivity.outcome_map(rs, mean_s_tac, x -> any(x .>= 0.5); S=20)
ADRIA.viz.outcome_map!(
    tf[1, 1],
    rs,
    tac_om_50,
    ["dhw_scenario", "wave_scenario", "N_seed_TA", "N_seed_CA", "fogging", "SRM"];
    axis_opts=Dict(:title => "Regions which lead to Top 50th Percentile Outcomes", :ylabel => "TAC [m²]")
)

# Indicate factor values that are in the top 30 percentile
tac_om_70 = ADRIA.sensitivity.outcome_map(rs, mean_s_tac, x -> any(x .>= 0.7); S=20)
ADRIA.viz.outcome_map!(
    tf[2, 1],
    rs,
    tac_om_70,
    ["dhw_scenario", "wave_scenario", "N_seed_TA", "N_seed_CA", "fogging", "SRM"];
    axis_opts=Dict(:title => "Regions which lead to Top 30th Percentile Outcomes", :ylabel => "TAC [m²]"))

save("outcome_map.png", tf)
```

![Outcome mapping](/ADRIA.jl/dev/assets/imgs/outcome_map.png?raw=true "Outcome mapping")

### GUI for high-level exploration (prototype only!)

```julia
# To explore results interactively
ADRIA.viz.explore("path to Result Set")

# or, if the result set is already loaded:
# ADRIA.viz.explore(rs)
```

![Standalone app for data exploration](/ADRIA.jl/dev/assets/imgs/aviz_app.png?raw=true "Data Exploration App")