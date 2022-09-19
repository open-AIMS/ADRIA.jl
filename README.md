# ADRIA.jl
![Tests](https://github.com/open-AIMS/ADRIA.jl/actions/workflows/ci.yml/badge.svg?branch=main)

Julia implementation of ADRIA: the Adaptive Dynamic Reef Intervention Algorithm.

Scripts showcasing example usage are found in the `examples` directory.


## Setup

To specify user-specific options, a `config.toml` file should be created with the following options (adjusted to suit your needs):

```toml
[operation]
num_cores = 2     # No. of cores to use. Values <= 0 will use all available cores.
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

**A.** Scenarios are defined in a CSV file (with parameter values in columns, so that each row defines a scenario).

  - See the `example_scenarios.csv` file in the `examples` directory for an idea of what this looks like.
  - See `parameters.jl` file in the `examples` directory on how to extract the model specification and parameter table for a given domain.
  - See also the `running_scenarios.jl` example script which showcases how to run such a file for a given study area.


## Testing

To run tests:

```bash
$ julia --project=.
julia> ]test 
```


## Notes

The very first import of the ADRIA package will be very slow as it attempts to precompile common functions to reduce later start up time.
The same applies when running ADRIA for the first time. This slow initial precompilation has to be repeated if the package is modified, but will remain "fast" if no changes are made.

To ameliorate this start-up cost while developing, use the [Revise package](https://github.com/timholy/Revise.jl).

A custom [sysimage](https://julialang.github.io/PackageCompiler.jl/dev/sysimages.html) can also be created to reduce start up times.
See the documentation [here](https://github.com/open-AIMS/ADRIA.jl/tree/main/build) for a quick how to.
Note: compilation time to create a sysimage can be upwards of 30mins, and has to be repeated if the included packages are to be updated.

```julia
# Timings here were taken with Julia v1.8

# Without custom sysimage
julia> @time using ADRIA
 97.883173 seconds (101.72 M allocations: 7.155 GiB, 4.22% gc time, 15.10% compilation time: 78% of which was recompilation)

# With custom sysimage
julia> @time using ADRIA
 0.012177 seconds (702 allocations: 40.648 KiB)
```
