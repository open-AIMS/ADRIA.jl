# ADRIAanalysis.jl

Post-simulation analysis tools for the [ADRIA](https://github.com/open-AIMS/ADRIA.jl) coral
reef decision support framework.

## Overview

ADRIAanalysis processes `ResultSet` outputs from ADRIA simulations to identify robust
intervention strategies, quantify factor sensitivities, and extract interpretable decision rules.

## Installation

```julia
using Pkg
Pkg.develop(path="path/to/ADRIAanalysis")
```

Requires Julia 1.11+ and ADRIA 0.17.

## Features

| Module | Functions | Purpose |
|--------|-----------|---------|
| Analysis | `find_pareto_optimal`, `find_robust` | Identify efficient/robust scenarios |
| Analysis | `data_envelopment_analysis` | DEA efficiency frontier |
| Analysis | `screen_scenarios`, `scenario_clusters`, `target_clusters` | Scenario filtering and clustering |
| Analysis | `feature_set` | Consolidate scenario inputs + DHW stats into a DataFrame |
| Analysis | `intervention_frequency` | Count seed/shade/fog deployments |
| Analysis | `cluster_rules`, `print_rules` | Extract interpretable decision rules (requires SIRUS, MLJ) |
| Sensitivity | `sensitivity.pawn` | PAWN sensitivity indices |
| Sensitivity | `sensitivity.tsa` | Two-Step Algorithm |
| Sensitivity | `sensitivity.rsa` | Regional Sensitivity Analysis |
| Sensitivity | `sensitivity.outcome_map` | Map outcomes to factor regions |

## Quick Start

```julia
using ADRIA, ADRIAanalysis, Statistics

rs = ADRIA.load_results("path/to/results")

# Compute outcomes
tac = ADRIA.metrics.scenario_total_cover(rs)
rsv = ADRIA.metrics.scenario_rsv(rs)

y = hcat(vec(mean(tac; dims=1)), vec(mean(rsv; dims=1)))

# Find Pareto-optimal scenarios across RCPs
optimal = find_pareto_optimal(rs, y, [45, 60])

# Filter to scenarios meeting a decision rule
robust = find_robust(rs, y, x -> all(x .>= 0.9), [45, 60])

# Count intervention deployments in robust scenarios
freqs = intervention_frequency(rs, robust, :seed)
```

### Sensitivity Analysis

```julia
X = feature_set(rs)

# PAWN sensitivity indices (10 bins)
si = sensitivity.pawn(X, vec(mean(tac; dims=1)), 10)

# Regional Sensitivity Analysis
sensitivity.rsa(X, vec(mean(tac; dims=1)), rs.model_spec)
```

### Data Envelopment Analysis

```julia
cost  = rs.inputs.cost
s_tac = dropdims(mean(tac; dims=:timesteps); dims=:timesteps)
s_sv  = dropdims(mean(ADRIA.metrics.scenario_shelter_volume(rs); dims=:timesteps); dims=:timesteps)

result = data_envelopment_analysis(cost, s_tac, s_sv)
```

## Optional Dependencies

Rule extraction (`cluster_rules`, `rules`) requires:

```julia
Pkg.add(["SIRUS", "MLJ"])
```

## License

See [LICENSE](LICENSE).
