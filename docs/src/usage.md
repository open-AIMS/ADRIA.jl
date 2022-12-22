# Using ADRIA

ADRIA is a decision support platform which includes decision heuristics (MCDA) with a small coral ecosystem model.
The MCDA processes may be used separately.

This section outlines how ADRIA may be used to arrive at a select range of possible pathways that are robust to
possible future conditions.

Julia may be installed through the package manager or from the github repository (for the most recent development version)

```julia-repl
julia> ]add ADRIA

# OR

julia> ]add https://github.com/open-AIMS/ADRIA.jl.git
```

Create a directory for your project. 

Before anything can be done, a `config.toml` file is needed.
Go to your project directory and create a `config.toml` file inside that directory
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

Start Julia from the project directory and import ADRIA to start.

```julia
# Import ADRIA package
using ADRIA
```

## Loading a Domain

ADRIA is designed to work with `Domain` data packages.
In short, these are pre-packaged data sets that hold the all necessary data to run
simulations for a given spatial domain.

See [Architecture](@ref) for more information.

A `Domain` may be loaded with the `load_domain` function.
By convention we assign the `Domain` to `dom`, although this variable can be named anything.

```julia
dom = ADRIA.load_domain("path to Input Set")
```

## Generating scenarios

Typical use of ADRIA is to generate a number of scenarios by sampling from combinations of
possible factors relating to environmental, intervention, and coral conditions.

A scenario is defined as a combination of all factors (i.e., all the model inputs).

```julia
# Generate 128 scenarios based on available environmental data layers and model parameters
scens = ADRIA.sample(dom, 128)
```

Here, the `scens` variable holds a DataFrame of scenarios of shape
``N`` by ``D``, where ``N`` is the number of scenarios (rows) and ``D``
is the number of factors.

!!! tip "Alternate samplers"
    The sampling method is compatible with any sampler supported by the 
    [Surrogates.jl](https://github.com/SciML/Surrogates.jl) package.

    The default sampler is the Sobol'. Below is an example of
    using Latin Hypercube sampling.

    ```julia
    import Surrogates.QuasiMonteCarlo: LatinHypercubeSample
    scens = ADRIA.sample(ex_domain, 100, LatinHypercubeSample())
    ```


### On model parameters and specifications

The current default values can be extracted with:

```julia
param_df = ADRIA.param_table(scenario_domain)
```

Again, `param_df` is a DataFrame of ``1`` by ``D``:
A single scenario with input values set to their default.

Running specific user-defined scenarios is as simple as modifying
the scenario DataFrame (referred to as the "scenario specification", 
or "scenario spec").

On a related note, details of the ADRIA model - parameter names, the 
default values, and their bounds - can be extracted as well.

```julia
# Get model specification
model_spec = ADRIA.model_spec(scenario_domain)

# Sometimes it is useful to export the model specification to CSV
ADRIA.model_spec(scenario_domain, "model_spec.csv")
```

## Running scenarios

```julia
# Run sampled scenarios for a given RCP
rs = ADRIA.run_scenarios(scens, dom, "45")

# ... or repeatedly run scenarios across multiple RCPs
rs = ADRIA.run_scenarios(scens, dom, ["45", "60", "85"])

# The location of the outputs stored on disk
@info ADRIA.store_name(rs)
# "Example_domain__RCPs45__2022-10-19_12_01_26_965"

@info ADRIA.store_location(rs)
# "[some location]/Example_domain__RCPs45__2022-10-19_12_01_26_965"
```

The `rs` variable is an `ResultSet` object which acts as an interface to the stored results.

The `ResultSet` provides:
- An overview of scenarios run
- Access to results from key ADRIA metrics
- Seeding/Shading/Fogging logs
- domain spatial data

```julia
print(rs)
```

!!! note "on-disk data store"
    ADRIA uses an on-disk data store (in Zarr format) to reduce memory use.
    The primary location for these is defined in the project's `config.toml` file
    (see instructions above).

!!! tip "Reloading results"
    Pre-existing results can also be reloaded by providing the path to the data store.

    ```julia
    rs = ADRIA.load_results("path to result set")
    ```

## Metrics

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

## Analysis

TODO