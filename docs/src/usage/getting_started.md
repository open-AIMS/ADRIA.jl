# Getting Started

ADRIA is a decision support platform which includes decision heuristics (MCDA) with a small coral ecosystem model.
The MCDA processes may be used separately.

This section outlines how ADRIA may be used to arrive at a select range of possible pathways that are robust to
possible future conditions.

Create a directory for your project, and start Julia inside that directory:

```bash
$ julia --project=.
```

ADRIA may be installed through the package manager or from the github repository (for the most recent development version).

```julia-repl
julia> ]add ADRIA

# OR, to install the latest development version:

julia> ]add https://github.com/open-AIMS/ADRIA.jl.git
```

Similarly, ADRIA can be updated as new releases are made:

```julia-repl
julia> ]up ADRIA
```

To setup ADRIA for development, see the [Development Setup](@ref) page.

## Visualizations

The Makie package ecosystem is used for plotting and need to be installed if visualizations are desired:

```julia-repl
julia> ]add GLMakie GeoMakie GraphMakie
```

## Running ADRIA simulations

Before anything can be done, a `config.toml` file is needed.
Create a `config.toml` file inside your project directory
with the following content below (adjusted for your needs).

This `config.toml` file is specific to your computer and project. It **should not be** committed to version control.

```toml
[operation]
num_cores = 2     # No. of cores to use. Values <= 0 will use all available cores.
threshold = 1e-8  # Result values below this will be set to 0.0 (to save disk space)
debug = false     # Disable multi-processing to allow error messages to be shown

[results]
output_dir = "./Outputs"  # Change this to point to where you want to store results
```

!!! tip "Performance"
    ADRIA uses an on-disk data store to hold results from model runs.
    Setting `output_dir` to a directory on an SSD (Solid State Drive)
    will maximize performance.

Start Julia from the project directory

```bash
$ julia --project=.
```


and import ADRIA to start.

```julia
# Import ADRIA package
using ADRIA
```

Below is an example workflow. Each line is explained in the next few pages.

```julia
using ADRIA

# Load a domain
dom = ADRIA.load_domain("path to some domain data package")

# Generate scenarios
scens = ADRIA.sample(128, dom)

# Run sampled scenarios for a given RCP
rs = ADRIA.run_scenarios(scens, dom, "45")

# ... or repeat scenario runs across multiple RCPs
rs = ADRIA.run_scenarios(scens, dom, ["45", "60", "85"])

# Obtain metrics
s_tac = ADRIA.metrics.scenario_tac(rs)
```

See [Analysis](@ref) for further examples.