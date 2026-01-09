# Getting Started

## Setup

This section outlines how ADRIA may be used to arrive at a select range of possible
pathways that are robust to possible future conditions.

Create a directory for your project, and start Julia inside that directory:

```bash
$ julia --project=.
```

ADRIA may be installed through the package manager or from the github repository (for the
most recent development version).

```julia-repl
julia> ]add ADRIA

# OR, to install the latest development version:

julia> ]add https://github.com/open-AIMS/ADRIA.jl.git
```

Similarly, ADRIA can be updated as new releases are made:

```julia-repl
julia> ]up ADRIA
```

If desired, you can create a `config.toml` file inside your project directory.
This is optional, and the assumed default values are shown below:

```toml
[operation]
num_cores = 2     # No. of cores to use. Values <= 0 will use all available cores.
threshold = 1e-8  # Result values below this will be set to 0.0 (to save disk space)
debug = false     # Disable multi-processing to allow error messages to be shown

[results]
output_dir = "./Outputs"  # Change this to point to where you want to store results
```

This `config.toml` file is specific to your computer and project. It **should not be**
committed to version control.


!!! tip "Performance"
    ADRIA uses an on-disk data store to hold results from model runs.
    Setting `output_dir` to a directory on an SSD (Solid State Drive)
    will maximize performance.

To setup ADRIA for development, see the [Development setup](@ref) page.

## Quick Start

A common workflow would be the following:

Start Julia from the project directory:

```bash
$ julia --project=.
```

Load data for a spatial domain. See [Loading a Domain](@ref) for more details:

```julia
using ADRIA

dom = ADRIA.load_domain("path to domain data package directory")
```

Generate scenarios based on available environmental data layers and model parameters. The
number of scenarios shoud be a power of two. See [Generating scenarios](@ref) for more
details:

```julia
num_scenarios = 128
scens = ADRIA.sample(dom, num_scenarios)
```

Run sampled scenarios for one or more RCPs. Be aware that this may take a while:

```julia
rcp_45 = "45"
rs = ADRIA.run_scenarios(dom, scens, rcp_45)
```

Or run scenarios across several RCPs:

```julia
rcps = ["45", "60", "85"]
rs = ADRIA.run_scenarios(dom, scens, rcps)
```

It is also possible to load previously run scenarios. See [Running scenarios](@ref) for
more details:

```julia
rs = ADRIA.load_results("path to results")
```

Extract some metric for analysis (e.g., the total absolute cover for each site and
timestep):

```julia
s_tc = ADRIA.metrics.scenario_total_cover(rs)
```

Use the visualization tools to plot the results. The Makie package ecosystem is used for
producing plots:

```julia
using GLMakie, GeoMakie, GraphMakie

# Plot a quick scenario overview
fig = ADRIA.viz.scenarios(rs, s_tc; axis_opts=Dict(:ylabel=>"Absolute Cover"))
save("path_to_save_figure", fig)
```

See [Analysis](@ref) for further examples of analysis and plots.

## Shared package depot paths

If parallel runs are to be conducted, it is recommended to set a shared `JULIA_DEPOT_PATH`.
This is so each individual worker does not race against each other to compile packages.

On Linux:

```shell
export JULIA_DEPOT_PATH="/some_shared_accessible_directory"
```

On Windows:

```shell
set JULIA_DEPOT_PATH="some_shared_accessible_directory"
```

https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH