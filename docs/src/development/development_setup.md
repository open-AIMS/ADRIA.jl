# Development setup

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

# Add additional debugging tools to sandbox environment
(sandbox) pkg> add Revise Infiltrator ProfileView BenchmarkTools JET

# Press ctrl+c to exit the package manager
```

Development scripts/functions can then be worked on in the `sandbox` folder without these polluting the ADRIA project itself.


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

!!! note "VS Code"
    VS Code now has (experimental support) for generating a custom sysimage for its REPL.
    The same caveats as above apply (it has to be recreated if the project specification has changed for any reason).

    See: [This guide](https://www.julia-vscode.org/docs/dev/userguide/compilesysimage/)


```julia
# Timings here were taken with Julia v1.8.1 for ADRIA v0.4.0

# Without custom sysimage
julia> @time using ADRIA
 34.289738 seconds (84.97 M allocations: 5.771 GiB, 6.44% gc time, 16.08% compilation time: 74% of which was recompilation)

# With custom sysimage
julia> @time using ADRIA
 0.012177 seconds (702 allocations: 40.648 KiB)
```
