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

## Workflow

TODO