# ADRIAviz

Visualization subpackage for [ADRIA.jl](https://github.com/open-AIMS/ADRIA.jl).

ADRIAviz provides plotting functions for ADRIA model outputs via the `ADRIA.viz` namespace.
It is backend-agnostic: Makie (GeoMakie + GraphMakie) and PlotlyBase are supported through
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
| PlotlyBase | `PlotlyBase` | `scenarios`, `clustered_scenarios`, `taxonomy`, `pawn`, `tsa`, `rsa`, `outcome_map`, `convergence`, `data_envelopment_analysis`, `rules_scatter` |
| PlotlyBase + PlotlyKaleido | `PlotlyBase` + `PlotlyKaleido` | All PlotlyBase functions + `savefig` to PNG/SVG/PDF |

Load any compatible Makie backend (`WGLMakie`, `GLMakie`, `CairoMakie`) alongside the
trigger packages to activate the full Makie extension.

### PlotlyBase backend

> **Display:** `PlotlyBase.Plot` auto-renders in Jupyter and Pluto. In VS Code with
> the Julia extension, `display(p)` opens the plot pane. In a plain terminal REPL,
> use `ADRIA.viz.show_in_browser(p)` to write a temp HTML file and open it in your
> default browser.

```julia
using ADRIA, ADRIAviz
ADRIAviz.activate("plotly")   # loads PlotlyBase + PlotlyKaleido (if installed)

rs = ADRIA.load_results("path/to/results")
s_tc = ADRIA.metrics.scenario_total_cover(rs)

# Returns a PlotlyBase.Plot
p = ADRIA.viz.scenarios(rs, s_tc)
p   # display inline in Jupyter / Pluto

# Sensitivity analysis
Si = ADRIA.sensitivity.pawn(rs, s_tc)
ADRIA.viz.pawn(Si)
ADRIA.viz.tsa(Si)
ADRIA.viz.rsa(Si, factor_values)

# Clustering
ADRIA.viz.clustered_scenarios(outcomes_matrix, cluster_labels)
ADRIA.viz.rules_scatter(scenarios_df, clusters, rules)

# Static export (requires PlotlyKaleido)
ADRIA.viz.savefig(p, "output.png")
ADRIA.viz.savefig(p, "output.svg")
ADRIA.viz.savefig(p, "output.pdf"; width=1200, height=800)
```

### Saving figures (PlotlyKaleido)

`savefig` is activated automatically by `ADRIAviz.activate("plotly")` when
`PlotlyKaleido` is installed. To install it:

```julia
] add PlotlyKaleido
```

Or load it explicitly alongside `PlotlyBase`:

```julia
using PlotlyBase, PlotlyKaleido

p = ADRIA.viz.scenarios(rs, s_tc)
ADRIA.viz.savefig(p, "output.png")                          # PNG 900×600 @ 3× scale
ADRIA.viz.savefig(p, "output.svg")                          # SVG (path extension auto-detected)
ADRIA.viz.savefig(p, "output.pdf"; width=1200, height=800)  # PDF, custom size
```

> **Note:** In-place (`!`) variants such as `scenarios!` are not supported by the
> Plotly backend. Use the non-mutating forms and compose layouts manually with
> `PlotlyBase.make_subplots`.

## Running tests

```julia
# From ADRIAviz/
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.test()'

# Include Makie backend tests (requires WGLMakie, GeoMakie, GraphMakie)
$env:ADRIA_RUN_VIZ_TESTS = "1"
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.test()'

# Include Plotly backend tests
$env:ADRIA_RUN_PLOTLY_TESTS = "1"
julia --project=. test/runtests.jl
```
