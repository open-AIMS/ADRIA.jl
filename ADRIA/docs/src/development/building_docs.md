# Building Documentation

ADRIA documentation is built using [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl)
with [DocumenterVitepress.jl](https://github.com/LuxDL/DocumenterVitepress.jl) as the site
theme and [Literate.jl](https://github.com/fredrikekre/Literate.jl) for tutorial pages.

Tutorial source files are `.jl` scripts in `ADRIA/docs/src/usage/`. Running `make.jl`
converts them to markdown before Documenter processes the full site. The generated `.md`
files are not tracked in version control.

## Refreshing documentation figures

Example figures shown in the documentation are PNG files committed under
`ADRIA/docs/src/assets/imgs/`. They are generated offline using the scripts in
`ADRIA/docs/scripts/` and should be re-run whenever the relevant visualizations change.

The scripts require a separate Julia environment with heavy dependencies
(CairoMakie, ADRIAanalysis, MLJ, SIRUS) that are not part of the main docs build.
Set it up once:

```julia
julia> ]activate ADRIA/docs/scripts
(scripts) pkg> dev ADRIA ADRIAanalysis ADRIAviz
(scripts) pkg> add CairoMakie GeoMakie GraphMakie MLJ SIRUS DataFrames
```

Set `ADRIA_TEST_DOMAIN` to the path of a local domain directory. You can export it
once per shell session to avoid repeating it:

**Bash / Linux / macOS:**
```bash
$ export ADRIA_TEST_DOMAIN=/path/to/domain
```

**PowerShell (Windows):**
```powershell
$env:ADRIA_TEST_DOMAIN = 'C:\path\to\domain'
```

**Command Prompt (Windows):**
```cmd
set ADRIA_TEST_DOMAIN=C:\path\to\domain
```

Then regenerate figures:

```bash
$ julia --project=ADRIA/docs/scripts ADRIA/docs/scripts/makie_viz_check.jl
```

PNGs are written to `ADRIA/docs/scripts/makie_viz_output/` by default. To write
directly to the committed docs assets path, override `ADRIA_FIGURE_DIR`:

```bash
$ ADRIA_FIGURE_DIR=ADRIA/docs/src/assets/imgs/analysis \
  julia --project=ADRIA/docs/scripts ADRIA/docs/scripts/makie_viz_check.jl
```

Updated PNGs are written directly to `ADRIA/docs/src/assets/imgs/analysis/` and should
be committed. Model outputs (Zarr result sets written to `Outputs/`) are gitignored.

## Setup

Add the required documentation dependencies if they are not already in the docs environment:

```julia
julia> ]activate ADRIA/docs
(docs) pkg> dev ADRIA ADRIAanalysis
(docs) pkg> add Literate DocumenterVitepress
```

## Building locally

From the repository root:

```bash
$ julia --project=ADRIA/docs/ --threads=auto
(docs) pkg> resolve
(docs) pkg> up
```

Then build:

```bash
$ julia --project=ADRIA/docs/ --threads=auto ADRIA/docs/make.jl
```

This generates markdown from the Literate sources and runs `npm install` automatically
using the Node.js bundled with DocumenterVitepress (`NodeJS_20_jll`).

On **Linux/macOS**, `makedocs` also runs the full Vitepress build and writes the final
static HTML to `ADRIA/docs/build/final_site/`.

On **Windows**, the Vitepress build step is skipped by `makedocs` (a current limitation
of the JLL runtime on Windows). To complete the build or start a live-reload dev server,
install [Node.js](https://nodejs.org/en) (v22.11.0 or later) and run from `ADRIA/docs/`:

```bash
$ npm run docs:build   # produce static HTML in build/final_site/
$ npm run docs:dev     # start live-reload dev server at http://localhost:5173
```

## Documentation deployment

Documentation is hosted on [GitHub Pages](https://pages.github.com/) via
[GitHub Actions](https://github.com/features/actions).

Configuration is in [`.github/workflows/documentation.yml`](https://github.com/open-AIMS/ADRIA.jl/blob/main/.github/workflows/documentation.yml).

Documentation is built and deployed automatically on commit to `main`. When a PR
targeting `main` is submitted, a preview URL is generated (e.g. `previews/PR###`).
