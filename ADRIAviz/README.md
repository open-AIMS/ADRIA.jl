# ADRIAviz

Visualization subpackage for [ADRIA.jl](https://github.com/open-AIMS/ADRIA.jl).

ADRIAviz provides plotting functions for ADRIA model outputs via the `ADRIA.viz` namespace.
It is backend-agnostic: Makie (GeoMakie + GraphMakie) and PlotlyLight are supported through
Julia's package extension mechanism — only the backends you load are compiled.

> **v0.16.0 breaking change:** ADRIAviz is now a separate package.
> Code that previously relied on `using ADRIA, CairoMakie` to activate viz methods must
> now add `using ADRIAviz` (or `Pkg.add("ADRIAviz")` if outside the monorepo).

## Installation

From within the ADRIA.jl monorepo:

```julia
using Pkg
Pkg.develop(path="ADRIAviz")
```

As a standalone package (once registered):

```julia
] add ADRIAviz
```

## Usage

```julia
using ADRIA, ADRIAviz
ADRIAviz.activate()               # loads WGLMakie + GeoMakie + GraphMakie
ADRIAviz.activate("CairoMakie")   # non-interactive / CI

rs = ADRIA.load_results("path/to/results")

# Scenario overview
s_tc = ADRIA.metrics.scenario_total_cover(rs)
fig = ADRIA.viz.scenarios(rs, s_tc)

# Spatial map
fig = ADRIA.viz.map(rs, s_tc)

# Sensitivity analysis
Si = ADRIA.sensitivity.pawn(rs, s_tc)
fig = ADRIA.viz.pawn(Si)

# Interactive explorer (Makie only)
ADRIA.viz.explore(rs)
```

Alternatively, load backends explicitly for full control:

```julia
using ADRIA, ADRIAviz
using WGLMakie, GeoMakie, GraphMakie
```

## Backends

| Backend | Trigger packages | Supported functions |
|---------|-----------------|---------------------|
| Makie   | `Makie` + `GeoMakie` + `GraphMakie` | All `ADRIA.viz.*` functions, `explore()` GUI |
| PlotlyLight | `PlotlyLight` | `scenarios`, `taxonomy` (planned) |

Load any compatible Makie backend (`WGLMakie`, `GLMakie`, `CairoMakie`) alongside the
trigger packages to activate the full Makie extension.

## Running tests

```julia
# From ADRIAviz/
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.test()'

# Include Makie backend tests
$env:ADRIA_RUN_VIZ_TESTS = "1"
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.test()'
```
