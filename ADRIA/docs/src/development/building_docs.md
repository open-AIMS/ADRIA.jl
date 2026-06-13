# Building Documentation

ADRIA documentation is built using [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl)
with [DocumenterVitepress.jl](https://github.com/LuxDL/DocumenterVitepress.jl) as the site
theme and [Literate.jl](https://github.com/fredrikekre/Literate.jl) for tutorial pages.

Tutorial source files are `.jl` scripts in `ADRIA/docs/src/usage/`. Running `make.jl`
converts them to markdown before Documenter processes the full site. The generated `.md`
files are not tracked in version control.

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
