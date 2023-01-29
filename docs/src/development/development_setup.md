# Development setup

Install Julia if not already done so.

## Recommendations

- It is recommended that the [juliaup](https://github.com/JuliaLang/juliaup) tool be used to ease managing Julia versions.
- We recommend [VS Code](https://code.visualstudio.com/) with its Julia extension when developing ADRIA.
- We also recommend the built-in Julia REPL within VS Code be used (see the notes below).


## Initial Setup

Once installed, clone the ADRIA.jl repository, navigate to the project folder, and start Julia.

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
(sandbox) pkg> add Revise Infiltrator BenchmarkTools JET

# Press ctrl+c to exit the package manager
```

Development scripts/functions can then be worked on in the `sandbox` folder, and its sub-folders, without these polluting the ADRIA project itself.


## Testing

To run tests:

```bash
$ julia --project=.
julia> ]test 
```


## Code Style

Follow the standard Julia [style guide](https://docs.julialang.org/en/v1/manual/style-guide/) as much as possible.

In most cases, simply auto-formatting the code is enough.


## Notes

The very first import of the ADRIA package will be very slow as it attempts to precompile common functions to reduce later start up time.
The same applies when running ADRIA for the first time. This slow initial precompilation has to be repeated if the package is modified, but will remain "fast" if no changes are made.

To ameliorate this start-up cost while developing, use the [Revise package](https://github.com/timholy/Revise.jl).

A custom [sysimage](https://julialang.github.io/PackageCompiler.jl/dev/sysimages.html) can also be created to reduce start up times.

VS Code now has (experimental) support for generating a custom sysimage for its REPL.
The same caveats as above apply: the sysimage has to be recreated if the project specification has changed for any reason.

It is highly recommended that this sysimage be built and used.

See: [This guide](https://www.julia-vscode.org/docs/dev/userguide/compilesysimage/)

Otherwise, see the documentation [here](https://github.com/open-AIMS/ADRIA.jl/tree/main/build) for a quick how to.
Note: compilation time to create a sysimage can be upwards of 30mins, and has to be repeated if the project packages are updated.
