# ADRIA.jl

Julia implementation of ADRIA (in progress).

In the short-term, functionality may leverage the existing MATLAB version and so a working copy of MATLAB and ADRIA is recommended.
Interaction is handled via the MATLAB.jl interface (See intro here: https://github.com/JuliaInterop/MATLAB.jl).

The `example.ipynb` in the `examples` directory showcases an example use of MATLAB.jl


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
(sandbox) pkg> add Revise Infiltrator ProfileView BenchMarkTools

# Press ctrl+c to exit the package manager
```

Development scripts/functions can then be worked on in the `sandbox` folder without these polluting the repository.


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
