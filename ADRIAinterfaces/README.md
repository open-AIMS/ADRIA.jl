# ADRIAinterfaces

Contracts between the [ADRIA.jl](https://github.com/open-AIMS/ADRIA.jl) run loop and optional
submodels (models that run alongside the coral model).

This package contains contracts only: an abstract `Submodel`, `Phase` types, hook functions
with no-op defaults (`init_state`, `step!`, `record!`, `log_specs`, `validate`, `prepare`,
`active`, `seed_excluded`, `load_inputs`, `parameters`), and small concrete types
(`ParamSpec`, `LogSpec`, `SubmodelContext`, `CoralStepView`). It has no ADRIA types, no model
code and no simulation loop, so model packages can depend on it without depending on ADRIA.

ADRIA depends on this package and calls the hooks at fixed points of a scenario run. See the
docstrings of `ADRIAinterfaces` and of each hook for the call order and mutation rules.

## Implementing a submodel

```julia
using ADRIAinterfaces
import ADRIAinterfaces: phases, step!

struct MySubmodel <: Submodel end
phases(::MySubmodel) = (AfterDisturbances(),)

function step!(::MySubmodel, ::AfterDisturbances, state, ctx, view, tstep, year)
    cover(view) .*= 0.99   # decrease only
    return nothing
end
```

## Installation

From within the ADRIA.jl monorepo:

```julia
using Pkg
Pkg.develop(path="ADRIAinterfaces")
```

## Tests

```
julia --project=ADRIAinterfaces -e 'using Pkg; Pkg.test()'
```

`test/MockSubmodels.jl` holds mock submodels that ADRIA's own tests reuse.
