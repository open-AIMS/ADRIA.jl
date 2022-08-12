# ADRIA.jl

Julia implementation of ADRIA under development.

Scripts showcasing example usage are found in the `examples` directory.


## Setup

To specify user-specific options, a `config.toml` file should be created with the following options (adjusted to suit your needs):

```toml
[operation]
num_cores = 2     # No. of cores to use. Values <= 0 will use all available cores.
reps = 20         # No. of environmental scenarios to run
threshold = 1e-6  # Result values below this will be set to 0.0

[results]
output_dir = "./Outputs"  # Change this to point to where you want to store results
```

## Development setup

```bash
# Start julia specifying the current directory as the project
$ julia --project=.

# Instantiate project. Sets up project packages. Only need to do this once.
julia> ]instantiate
```

For development purposes, set up a sandbox environment **(setup only needs to be done once)**.
The steps below assumes you are in the ADRIA.jl project folder.

```bash
$ mkdir sandbox
$ cd sandbox
$ julia

# Switch to package management (`]`) and activate sandbox environment
julia> ]activate .

# Add ADRIA.jl as a local package under development
(sandbox) pkg> dev ../

# Add debugging tools to sandbox environment
(sandbox) pkg> add Revise Infiltrator ProfileView BenchmarkTools

# Press ctrl+c to exit the package manager
```

Development scripts/functions can then be worked on in the `sandbox` folder without these polluting the ADRIA project itself.


## Troubleshooting F.A.Q

ADRIA.jl is under active development and from time to time issues may arise.
Here are some answers to some issues encountered.

**Q. I get this warning when trying to load pre-existing results:**
  `Results were produced with an older version of ADRIA (v0.x.x). The installed version of ADRIA is: v0.y.y. Errors may occur when analyzing data.` 
  (where `x` and `y` are different numbers).

**A.** The result set being loaded were created by an older version of ADRIA, and stored in an older, possibly incompatible, format.
  Sometimes, results may still be produced/analyzed as normal. In other times, ADRIA.jl or the expected metadata in the result set may have changed
  leading to errors when conducting analyses.

  Either go back to the version indicated, or re-run the scenarios to obtain results in the updated format.

**Q. I get an error or warning about an `ENV` variable not being found or set.**

**A.** Double check the configuration settings in `config.toml` (see above).

**Q. How do I run my own scenarios?**

**A.** Scenarios are defined in a CSV file of parameter values (with values in columns, so that each row defines a scenario).

  - See `parameters.jl` file in the `examples` directory on how to extract the model specification and parameter table for a given domain.
  - See the `example_scenarios.csv` file in the `examples` directory for an idea of what this looks like.
  - See also the `running_scenarios.jl` example script which showcases how to run such a file for a given study area.


## Notes

The very first import of the ADRIA package will be very slow as it attempts to precompile common functions to reduce later start up time.
This slow initial precompilation has to be repeated if the package is modified, but will remain "fast" if no changes are made.

```julia
# First time importing ADRIA
julia> @time using ADRIA
 25.599935 seconds (40.47 M allocations: 2.990 GiB, 6.53% gc time, 22.64% compilation time)

# After precompilation
julia> @time using ADRIA
  0.002528 seconds (429 allocations: 27.859 KiB)
```

The same applies when running ADRIA for the first time:

```julia
# Running the first time
julia> @time include("running_scenarios.jl")
[ Info: Loading data package
[ Info: Loading example scenarios
[ Info: Setting up and running scenarios
[ Info: Reloading results and saving figure
441.297525 seconds (402.13 M allocations: 24.268 GiB, 2.05% gc time, 53.47% compilation time)

# Running a second time
julia> @time include("running_scenarios.jl")
[ Info: Loading data package
[ Info: Loading example scenarios
[ Info: Setting up and running scenarios
[ Info: Reloading results and saving figure
  5.066090 seconds (1.66 M allocations: 2.145 GiB, 6.70% gc time, 0.00% compilation time)
```

To ameliorate this start-up cost while developing, use the [Revise package](https://github.com/timholy/Revise.jl).

# App

[IN PROGRESS]

A command-line application can be created with PackageCompiler.jl

```julia
create_app("./", "app"; include_lazy_artifacts=true, force=true)
```

```bash
# Run scenarios and generate png of indicative results
$ app/bin/ADRIA.exe run [path to data package] [RCP scenario] [CSV of scenarios to run]

# Re-create plot of indicative results
$ app/bin/ADRIA.exe load [path to result set]

# Display options
$ app/bin/ADRIA.exe help
```
