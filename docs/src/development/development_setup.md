# Development setup

## Install Julia

Install Julia if not already done so, preferably using [juliaup](https://github.com/JuliaLang/juliaup).

### Windows

There may be issues installing juliaup from the Windows Store (or otherwise undesirable).
In such cases, use the alternate MSIX App installer (https://install.julialang.org/Julia.appinstaller)


## Recommendations

- It is recommended that the [juliaup](https://github.com/JuliaLang/juliaup) tool be used to ease managing Julia versions.
- We recommend [VS Code](https://code.visualstudio.com/) with its Julia extension when developing ADRIA.
- We also recommend the built-in Julia REPL within VS Code be used (see the notes below).


## Initial Setup

### Using the console

Once Julia is installed, clone the ADRIA.jl repository and navigate to the project folder:

```bash
$ git clone git@github.com:open-AIMS/ADRIA.jl.git
$ cd ./ADRIA.jl
```

Start Julia specifying the current directory as the project environment:

```bash
$ julia --project=.
```

Switch to the package manager (`]`) and instantiate the project. **This only needs to be done once**.

```julia-REPL
julia> ]
(ADRIA.jl) pkg> instantiate
```

This will sets up the project packages.

## Sandbox

For development purposes, set up a sandbox environment **(setup only needs to be done once)**.
This environment will function as a project apart from ADRIA, where you can install any
packages, including ADRIA.jl, and run your code. When installing ADRIA.jl at the sandbox,
use the `dev` command instead of `add`. For more information, please refer to
[Pkg.jl documentation](https://pkgdocs.julialang.org/v1/managing-packages/#developing).

Once you are inside ADRIA.jl project folder, create a folder named `sandbox` and start
julia inside it:

```bash
$ mkdir sandbox
$ cd sandbox
$ julia --project=.
```

Switch to the package manager (`]`) and add ADRIA.jl as a local package under development

```julia-REPL
julia> ]
(sandbox) pkg> dev ../
```

You may also install additional packages for visualizations and debugging tools

```julia-REPL
(sandbox) pkg> add GLMakie GeoMakie GraphMakie
(sandbox) pkg> add Revise Infiltrator BenchmarkTools JET
```

Press backspace or Ctrl+C to leave the package manager.

Development scripts/functions can now be worked on in the `sandbox` folder, and its
sub-folders, without these polluting the ADRIA project itself.


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

We are currently moving to follow [Blue Style Guide](https://github.com/invenia/BlueStyle).
All PRs should follow this style guide.

## Notes

The very first import of the ADRIA package will be very slow as it attempts to precompile common functions to reduce later start up time.
The same applies when running ADRIA for the first time. This slow initial precompilation has to be repeated if the package is modified, but will remain "fast" if no changes are made.

Use the [Revise package](https://github.com/timholy/Revise.jl) to avoid having to repeatedly restart the REPL to incorporate code changes.
By default, the VS Code REPL will auto-load this package.

A custom [sysimage](https://julialang.github.io/PackageCompiler.jl/dev/sysimages.html) can also be created to reduce start up times.

VS Code now has (experimental) support for generating a custom sysimage for its REPL.
Prior to Julia v1.9, a custom sysimage for the development/sandbox environment was highly recommended.
Julia v1.9 introduced an improved precompilation process and the concept of extension packages.
As many packages are still in the process of taking advantage of these changes, the sysimage may not
successfully build. Given precompilation is now much faster than previously, the sysimage can be
considered to be a "nice to have".

The same caveats as above apply: the sysimage has to be recreated if the project specification (e.g., expected package dependencies) changes.

See: [This guide](https://www.julia-vscode.org/docs/dev/userguide/compilesysimage/)

Otherwise, if the VS Code build task cannot be used, see the documentation [here](https://github.com/open-AIMS/ADRIA.jl/tree/main/build) for a quick how to.
Note: compilation time to create a sysimage can be upwards of 15mins and, again, has to be repeated if the project packages are updated.