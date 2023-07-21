# Development setup

Install Julia if not already done so, preferably using [juliaup](https://github.com/JuliaLang/juliaup).

There may be issues installing juliaup from the Windows Store (or otherwise undesirable).
In such cases, use the alternate MSIX App installer (https://install.julialang.org/Julia.appinstaller)


## Recommendations

- It is recommended that the [juliaup](https://github.com/JuliaLang/juliaup) tool be used to ease managing Julia versions.
- We recommend [VS Code](https://code.visualstudio.com/) with its Julia extension when developing ADRIA.
- We also recommend the built-in Julia REPL within VS Code be used (see the notes below).


## Initial Setup

Once installed, clone the ADRIA.jl repository, navigate to the project folder, and start Julia.
**This only needs to be done once**.

```bash
# If not using the VS Code REPL, start julia specifying the current directory as the project environment
# Assumes you have started the terminal in the ADRIA directory.
$ julia --project=.

# Instantiate project. Sets up project packages. Only need to do this once.
julia> ]instantiate
```

For development purposes, set up a sandbox environment **(setup only needs to be done once)**.
The steps below assumes you are in the ADRIA.jl project folder.

```bash
$ mkdir sandbox
$ cd sandbox
$ julia --project=.

# Switch to the package manager (`]`)
julia> ]

# Add ADRIA.jl as a local package under development
(sandbox) pkg> dev ../

# Install additional packages for visualizations
(sandbox) pkg> add GLMakie GeoMakie GraphMakie

# Add additional debugging tools to sandbox environment
(sandbox) pkg> add Revise Infiltrator BenchmarkTools JET

# Press ctrl+c to exit the package manager
```

Development scripts/functions can then be worked on in the `sandbox` folder, and its sub-folders, without these polluting the ADRIA project itself.


## Testing

To run the full test suite, rebuilding the environment as necessary:

```bash
$ julia --project=.
julia> ]test 
```

Rebuilding the environment can be unnecessary for every test run during development.
It such cases, `include()` the `runtests.jl` file directly.

```bash
# Assuming the current working directory is the project root.
# Adjust the filepath as necessary if this is not the case.
include("test/runtests.jl")
```

If a specific test case is being run, write the test file to be a standalone script
(importing all necessary packages, including `Test`) and run it directly.

Doing so allows use of debugging packages if necessary.

```bash
include("test/some_test_file.jl")
```

Once the test is complete, put the tests in a `testset` as and if required.
If a new file is added to the test suite, `include()` it in `test/runtests.jl`

See [the Test documentation](https://docs.julialang.org/en/v1/stdlib/Test/#Basic-Unit-Tests)
for further details.


## Code Style

Follow the standard Julia [style guide](https://docs.julialang.org/en/v1/manual/style-guide/) as much as possible.

In most cases, simply auto-formatting the code is enough.


## Notes

The very first import of the ADRIA package will be very slow as it attempts to precompile common functions to reduce later start up time.
The same applies when running ADRIA for the first time. This slow initial precompilation has to be repeated if the package is modified, but will remain "fast" if no changes are made.

To ameliorate this avoid having to repeatedly restart the REPL to incorporate code changes, use the [Revise package](https://github.com/timholy/Revise.jl).
By default, the VS Code REPL will auto-load this package.

A custom [sysimage](https://julialang.github.io/PackageCompiler.jl/dev/sysimages.html) can also be created to reduce start up times.

VS Code now has (experimental) support for generating a custom sysimage for its REPL.
The same caveats as above apply: the sysimage has to be recreated if the project specification has changed for any reason.

It is highly recommended that this sysimage be built and used.

See: [This guide](https://www.julia-vscode.org/docs/dev/userguide/compilesysimage/)

Otherwise, see the documentation [here](https://github.com/open-AIMS/ADRIA.jl/tree/main/build) for a quick how to.
Note: compilation time to create a sysimage can be upwards of 30mins and, again, has to be repeated if the project packages are updated.
