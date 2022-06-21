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

For development purposes, set up a sandbox environment.
Steps below assumes you are in the ADRIA.jl project folder.

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


## Note

The very first import of the ADRIA package will be very slow as it attempts to precompile common functions to reduce later start up time.
This slow initial precompilation has to be repeated if the package is modified, but will remain "fast" if no changes are made.

To ameliorate this start-up cost while developing, use the [Revise package](https://github.com/timholy/Revise.jl).

```julia
julia> @time using ADRIA
[ Info: Precompiling ADRIA [7dc409a7-fbe5-4c9d-b3e2-b0c19a6ba600]
111.073238 seconds (63.24 M allocations: 4.103 GiB, 1.96% gc time, 9.91% compilation time)

julia> @time using ADRIA
  0.002154 seconds (430 allocations: 28.562 KiB)
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
