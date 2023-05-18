# ADRIA.jl

ADRIA: Adaptive Dynamic Reef Intervention Algorithms.

[![Release](https://img.shields.io/github/v/release/open-AIMS/ADRIA.jl)](https://github.com/open-AIMS/ADRIA.jl/releases)  [![Documentation](https://img.shields.io/badge/docs-stable-blue)](https://open-aims.github.io/ADRIA.jl/stable/)  [![DOI](https://zenodo.org/badge/483052659.svg)](https://zenodo.org/badge/latestdoi/483052659)

ADRIA is a decision-support tool designed to help reef managers, modellers and decision-makers
address the challenges of adapting to climate change in coral reefs. It provides line of sight
to conservation solutions in complex settings where multiple objectives need to be considered,
and helps investors identify which options represent the highest likelihood of providing
returns on investment. ADRIA uses a set of dynamic Multi-Criteria Decision Analyses (dMCDA)
which simulates a reef decision maker to identify candidate locations for intervention
deployment which consider ecological, economic and social benefits.

ADRIA also includes a simplified coral ecosystem model to allow exploration of outcomes as a 
result of intervention decisions made across a wide range of possible future conditions.

Usage is demonstrated in the [documentation](https://open-aims.github.io/ADRIA.jl/stable/usage/getting_started/)


## Quick start

### Setup

To specify user-specific options, a `config.toml` file should be created with the following options (adjusted to suit your needs):

```toml
[operation]
num_cores = 2     # No. of cores to use. Values <= 0 will use all available cores.
threshold = 1e-6  # Result values below this will be set to 0.0
debug = false     # Disable multi-processing to allow error messages to be shown

[results]
output_dir = "./Outputs"  # Change this to point to where you want to store results
```

!!! tip "Performance"
    ADRIA uses an on-disk data store to hold results from model runs.
    Setting `output_dir` to a directory on an SSD (Solid State Drive)
    will maximize performance.


### Usage

```julia
# Import ADRIA package
using ADRIA


# Load data for a spatial domain
dom = ADRIA.load_domain("path to domain data package")

# Generate 128 scenarios based on available environmental data layers and model parameters
scens = ADRIA.sample(dom, 128)

# Run sampled scenarios for a given RCP
rs = ADRIA.run_scenarios(scens, dom, "45")

# ... or repeatedly run scenarios across several RCPs
rs = ADRIA.run_scenarios(scens, dom, ["45", "60", "85"])

# then extract metrics for analysis
tac = ADRIA.metrics.total_absolute_cover(rs)
```

# Analysis

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

## PAWN sensitivity (heatmap overview)

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
