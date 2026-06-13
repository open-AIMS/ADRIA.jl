# Development setup

## Install Julia

Install Julia if not already done so, preferably using [juliaup](https://github.com/JuliaLang/juliaup).

### Windows

There may be issues installing juliaup from the Windows Store (or otherwise undesirable).
In such cases, use the alternate MSIX App installer (https://install.julialang.org/Julia.appinstaller)


## Recommendations

- It is recommended that the [juliaup](https://github.com/JuliaLang/juliaup) tool be used to ease managing Julia versions.
- We recommend [VS Code](https://code.visualstudio.com/) with its Julia extension when developing ADRIA.
- Install the VS Code Julia Formatter extension (note: **not** the JuliaFormatter.jl package).
- We also recommend the built-in Julia REPL within VS Code be used (see the notes below).
- If parallel runs will be done, it is recommended to set a shared [`JULIA_DEPOT_PATH`](https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH)

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

```julia
julia> ]
(ADRIA.jl) pkg> instantiate
```

This will sets up the project packages.

## Sandbox

### Concept

The `sandbox` is a personal Julia project environment that lives inside the ADRIA.jl
repository folder but is excluded from version control (it appears in `.gitignore`).
Think of it as a workbench: a place to experiment with new features, prototype algorithms,
or run study-specific analyses without those scripts ever touching ADRIA's own git history.

Using the `dev` command (rather than `add`) to install ADRIA into the sandbox means that
any edits you make to ADRIA's source files are reflected immediately in your sandbox
scripts -- no reinstall required.

A typical sandbox layout looks like this:

```
ADRIA.jl/
|-- ADRIA/               <- ADRIA package source (version-controlled)
|-- sandbox/             <- your personal workspace (git-ignored)
    |-- Project.toml
    |-- Manifest.toml
    |-- scratch.jl       <- ad-hoc scripts
    |-- study_A/         <- a self-contained study
    |   |-- .git/        <- optionally its own git repository
    |   |-- ...
    |-- study_B/
        |-- ...
```

Sub-directories such as `study_A/` can be individually version-controlled with their
own `git init`, letting you track and share the code for a specific study while keeping
it completely separate from the main ADRIA repository.

### Alternative: external project

An equally valid approach is to create a project directory *outside* the ADRIA.jl
repository and `dev`-install ADRIA using its absolute path:

```julia
(my_study) pkg> dev /path/to/ADRIA.jl/ADRIA
(my_study) pkg> dev /path/to/ADRIA.jl/ADRIAanalysis
```

This is a good option if you prefer to keep study code entirely separate from the ADRIA
clone. The trade-off is that you maintain separate directories (and potentially separate
repositories) alongside the ADRIA clone, which can become harder to manage as the number
of studies grows.

### Setup (one-time)

From the ADRIA.jl repository root, create the `sandbox` folder and start Julia inside it:

```bash
$ mkdir sandbox
$ cd sandbox
$ julia --project=. --threads=auto
```

Switch to the package manager (`]`) and add ADRIA.jl as a local package under development:

```julia
julia> ]
(sandbox) pkg> dev ../ADRIA
(sandbox) pkg> dev ../ADRIAanalysis
```

You may also install additional packages for visualizations and debugging tools:

```julia
(sandbox) pkg> add ADRIAviz GLMakie GeoMakie GraphMakie
(sandbox) pkg> add Revise Infiltrator BenchmarkTools JET
```

Press backspace or Ctrl+C to leave the package manager.

Development scripts can now be worked on inside `sandbox/` and its sub-folders without
affecting the ADRIA project itself.


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
All PRs should follow this style guide. The [julia-format](https://github.com/julia-actions/julia-format)
GitHub workflow will check that your PR's code is formatted. Note that this check requires
all code in the repo to be formatted, not only the files modified by your PR.

To set things up:

1. Do **not** install the VS Code Julia Formatter extension. \
   If you have installed it, remove it.

2. Open VS Code settings and search for `default formatter`

3. Set it to `Julia` (julialang.language-julia)

Use the VSCode `Format Document` or `Format Selection` actions to format your code
(`shift+alt+f` is the shortcut, or `ctrl+shift+p` and search for `format`).

The applied formatting is defined in the project `.JuliaFormatter.toml` file.

### VSCode Settings

Open Settings (`Ctrl+,`). Search for `trim` and ensure the following options are all ticked/enabled:

- Files: Trim Final Newlines
- Files: Trim Trailing Whitespace
- Editor: Trim Auto Whitespace

An optional, but recommended, step would be to add a ruler guide to indicate where the
character limit/width is.

Search for `rulers` and click on "Edit in settings.json" under "Editor: Rulers"

Add "92" to the list of ruler lengths, such that the `editor.rulers` entry looks like this:

```json
"editor.rulers": [
        92
    ]
```

Adding multiple values adds more guide lines at the indicated widths.

*Important:* if you installed the *oh7z Julia Formatter* VSCode extension, uninstall or disable it for this workspace.
That extension always uses its formatter settings and does not support `.JuliaFormatter.toml` whereas the main Julia extension does.
The only reason to use the oh7z extension is for Julia projects that do not have a `.JuliaFormatter.toml` file.

### With the JuliaFormatter package

To reformat the entire project:

```julia
using JuliaFormatter
format(".")
```

*If this returns `false`, call `format()` again.*

Formatter configuration is defined in `.JuliaFormatter.toml`, see
[JuliaFormatter docs](https://domluna.github.io/JuliaFormatter.jl/stable/).

## Git blame ignore revs

If you have GitLens (or similar extension), it will show the author of a line of code using git blame.
To ignore commits like "re-format entire project", run this once:

`git config --local blame.ignoreRevsFile .git-blame-ignore-revs`

If you reformat code and want to ignore that commit, simply add the commit SHA to the `.git-blame-ignore-revs` file.

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
