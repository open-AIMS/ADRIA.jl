# ADRIA.jl

ADRIA: Adaptive Dynamic Reef Intervention Algorithms.

[![Release](https://img.shields.io/github/v/release/open-AIMS/ADRIA.jl)](https://github.com/open-AIMS/ADRIA.jl/releases) [![DOI](https://zenodo.org/badge/483052659.svg)](https://zenodo.org/badge/latestdoi/483052659)

[![Documentation](https://img.shields.io/badge/docs-stable-blue)](https://open-aims.github.io/ADRIA.jl/stable/) [![Documentation](https://img.shields.io/badge/docs-dev-blue)](https://open-aims.github.io/ADRIA.jl/dev/)

<div >
  <a href="https://github.com/invenia/BlueStyle">
    <img valign="bottom" src="https://img.shields.io/badge/code%20style-blue-4495d1.svg" alt="Code Style Blue">
    </a>
    <p>This project is currently moving to use <a href="https://github.com/invenia/BlueStyle">Blue Style Guide</a>.
Follow the style guide when submiting a PR.</p>
</div>

ADRIA is a decision-support tool designed to help reef managers, modellers and decision-makers
address the challenges of adapting to climate change in coral reefs. It provides line of sight
to conservation solutions in complex settings where multiple objectives need to be considered,
and helps investors identify which options represent the highest likelihood of providing
returns on investment. ADRIA uses a set of dynamic Multi-Criteria Decision Analyses (dMCDA)
which simulates a reef decision maker to identify candidate locations for intervention
deployment which consider ecological, economic and social benefits.

ADRIA also includes a simplified coral ecosystem model to allow exploration of outcomes as a
result of intervention decisions made across a wide range of possible future conditions.

ADRIA requires Julia v1.9 and above.

Setup and usage is demonstrated in the
[Documentation](https://open-aims.github.io/ADRIA.jl/stable/usage/getting_started/).

For developers, refer to the
[Developers setup guide](https://open-aims.github.io/ADRIA.jl/stable/development/development_setup/).

# ADRIA Docker Setup

This repo contains a multistage Dockerfile for generating containerised applications based on the ADRIA source code in this repository.

- [docker configuration](#docker-configuration)
- [adria-base](#adria-base)
  - [Building `adria-base`](#building-adria-base)
  - [Running `adria-base` with a non-interactive Julia command](#running-adria-base-with-a-non-interactive-julia-command)
  - [Running `adria-base` as an interactive Julia shell](#running-adria-base-as-an-interactive-julia-shell)
  - [Running `adria-base` with a non-Julia entrypoint](#running-adria-base-with-a-non-julia-entrypoint)
  - [Deriving an image from `adria-base`](#deriving-an-image-from-adria-base)
- [adria-dev](#adria-dev)
  - [Working with the `ADRIA.jl` submodule](#working-with-the-adriajl-submodule)
  - [Building `adria-dev`](#building-adria-dev)
  - [Running `adria-dev` as an interactive Julia shell](#running-adria-dev-as-an-interactive-julia-shell)
  - [Running ADRIA tests with `adria-dev`](#running-adria-tests-with-adria-dev)
- [adria-sandbox](#adria-sandbox)
  - [Building `adria-sandbox`](#building-adria-sandbox)
  - [Running `adria-sandbox`](#running-adria-sandbox)
  - [Interacting with `adria-sandbox`](#interacting-with-adria-sandbox)

---

## Docker configuration

The following build args and defaults are available to configure the build behaviour.

- ARG SANDBOX_FROM="adria-dev": What base image should be used for the sandbox

- ARG ADRIA_VERSION="0.11.0": What version of ADRIA from package registry to install in adria-base

- ARG JULIA_VERSION="1.10.4": See https://hub.docker.com/\*/julia for valid versions.

## adria-base

The `adria-base` image variant is a Julia image with [ADRIA.jl](https://github.com/open-AIMS/ADRIA.jl) installed to a shared environment that is included in the Julia LOAD_PATH directly from it's github origin.

The Docker entrypoint for the `adria-base` image is the `julia` binary, so you can run just this base container if you want to invoke julia commands, including commands that depend on the ADRIA package being installed.

It can also be used as a base-image for any Julia application that needs to use ADRIA.

### Building `adria-base`

```bash
# EITHER using docker compose:
docker compose build adria-base

# OR using just `docker build`:
docker build --target "adria-base" --tag ADRIA.jl/adria-base:latest .
```

You can also opt to specify some custom [build arguments](https://docs.docker.com/reference/cli/docker/image/build/#build-arg) to change the versions of Julia or ADRIA.jl that get installed. Supported arguments are:

- `ADRIA_REPO`: URL for the repository that ADRIA.jl should be cloned from. Defaults to <https://github.com/open-AIMS/ADRIA.jl.git>
- `ADRIA_REFSPEC`: the branch-name or tag of the `ADRIA_REPO` that you want to install. Defaults to `main`.
- `JULIA_VERSION`: The version of the Julia platform you want to install ADRIA.jl into. This must be one of the versions available for the official [Julia base image](https://hub.docker.com/_/julia). Defaults to `1.10.1`.

See the [docker-compose.yaml](./docker-compose.yaml) file for an example of how to specify build arguments in docker compose.

&nbsp;

### Running `adria-base` with a non-interactive Julia command

e.g. to list the packages installed in the `@adria` shared environment:

```bash
# EITHER using docker compose:
docker compose run --rm adria-base --project=@adria -e 'using Pkg; Pkg.status()'

# OR using just docker run:
docker run --rm ADRIA.jl/adria-base:latest --project=@adria -e 'using Pkg; Pkg.status()'
```

&nbsp;

### Running `adria-base` as an interactive Julia shell

To launch an interactive Julia shell:

```bash
# EITHER using docker compose:
docker compose run --rm adria-base

# OR using just docker run:
docker run --rm --interactive --tty ADRIA.jl/adria-base:latest
```

In both cases, type `CTRL-d` to exit the shell and stop the container.

&nbsp;

### Running `adria-base` with a non-Julia entrypoint

If you want to use this image to run something _other_ than a Julia command, you can specify an alternate entrypoint at runtime as well as an alternate command. e.g. to launch an interactive `bash` shell in the container for checking filesystem permissions or similar:

```bash
# EITHER using docker compose:
docker compose run --rm --entrypoint /bin/bash adria-base

# OR using just docker run:
docker run --rm --interactive --tty --entrypoint /bin/bash ADRIA.jl/adria-base:latest
```

&nbsp;

### Deriving an image from `adria-base`

To make a derived ADRIA application:

- Use `FROM ADRIA.jl/adria-base:latest` (or a related tag)
- Include a `CMD` line in your Dockerfile that provides appropriate arguments to the `julia` command line to invoke your application.

The section of the [Dockerfile](./Dockerfile) that defines the `adria-sandbox` target described below might be useful inspiration.

&nbsp;

---

## adria-dev

The `adria-dev` image variant is an _alternative_ to the `adria-base` image, not a derived application.
Instead of installing [ADRIA.jl](https://github.com/open-AIMS/ADRIA.jl) as a normal package, it looks for
the ADRIA.jl source code in a local subdirectory, and installs that as a [Julia development package](https://pkgdocs.julialang.org/v1/api/#Pkg.develop).

This allows you to use the `adria-dev` container as an `ADRIA.jl` development environment: you can run tests,
bind-mount and edit the code, re-resolve dependencies and all sorts of other useful things without needing
to a native installation of Julia.

### Working with the `ADRIA.jl` submodule

To support this usage, the ADRIA.jl repository has been configured as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules)
of this repository. This ensures the source code is always available in the proper location.

If you are not familiar with git submodules, some useful commands are:

- revert to the last used version of ADRIA.jl: `git submodule update --init`
- to switch to an existing branch of ADRIA.jl: `cd ADRIA.jl && git checkout <branch-name>`
- to pull recent commits on your current branch of ADRIA.jl: `cd ADRIA.jl && git pull`
- To commit and push changes you have made to ADRIA.jl source:
  - First, commit and push just as usual from _inside_ the submodule repository. As per good development practice, work in a dedicated branch:
    ```bash
    cd ADRIA.jl
    git checkout -b my-dev-branch
    git add .
    git commit -m"explain your changes here"
    git push --set-upstream origin my-dev-branch
    ```
  - At this point, your changes will be safely pushed to a new branch of the
    [ADRIA.jl origin repository](https://github.com/open-AIMS/ADRIA.jl), and you can make pull requests just like normal.
  - Then from the root of _this_ repository:
    ```bash
    git add ADRIA.jl
    git commit -m"explain your changes again"
    ```

&nbsp;

### Building `adria-dev`

Once you have ensured the ADRIA.jl source code is available, you can build the image like so:

```bash
# EITHER using docker compose:
docker compose build adria-dev

# OR using just `docker build`:
docker build --target "adria-dev" --tag ADRIA.jl/adria-dev:latest .
```

The same `JULIA_VERSION` build argument that works with `adria-base` will also work with `adria-dev`.

&nbsp;

### Running `adria-dev` as an interactive Julia shell

Very useful for running commands to update Package manifests and similar!
The command will activate the shared `@adria` environment by default, but
you can switch to any other environment as you need to.

```bash
# EITHER using docker compose:
docker compose run --rm adria-dev

# OR using just docker run:
docker run --rm --interactive --tty \
  --mount type=bind,source="$(pwd)"/ADRIA.jl/,target=/usr/local/src/adria/
  ADRIA.jl/adria-dev:latest
```

In both cases, type `CTRL-d` to exit the shell and stop the container.

&nbsp;

---

### Running ADRIA tests with `adria-dev`

ADRIA tests can be run like any other Julia command, but the GOTCHA is that
they need to run in the local project environment, NOT in the shared `@adria`
environment:

```bash
# EITHER using docker compose:
docker compose run --rm adria-dev --project=. -e 'using Pkg; Pkg.test();'

# OR using just docker run:
docker run --rm ADRIA.jl/adria-dev:latest --project=. -e 'using Pkg; Pkg.test();'
```

This method of running tests is suitable to use in a containerised continuous integration pipeline.

&nbsp;

---

## adria-sandbox

The `adria-sandbox` image variant is set up to run the `sandbox` Julia application which has its source code in this repository.

This application uses the pre-installed ADRIA package, and can be built on _either_ of `adria-base` or `adria-dev`,
depending on whether you want the development version of the package or not.

It depends on input and output data files which must be provided at runtime.

### Building `adria-sandbox`

The sandbox build supports an optional `SANDBOX_FROM` build argument which is used to specify a base image.
Supported values are `adria-dev` (the default) or `adria-base`.

```bash
# EITHER using docker compose (edit the compose file to specify which base image to use):
docker compose build adria-sandbox

# OR using just `docker build`:
docker build --build-arg SANDBOX_FROM=adria-dev --target "adria-sandbox" --tag ADRIA.jl/adria-sandbox:latest .
```

&nbsp;

### Running `adria-sandbox`

The `adria-sandbox` image is configured to automatically run the `dev.jl` script when a container
made from that image runs.

For that script to be useful, though, you must first mount two filesystem locations to the container:

- A directory containing your input data files should be mounted at `/data/input`
- A directory where the sandbox application can create output files should be mounted at `/data/output`

This documentation and the [docker-compose.yaml](./docker-compose.yaml) file in the repository demonstrate
using [bind mounts](https://docs.docker.com/storage/bind-mounts/) to data-directories which are local to
your current working directory for this purpose, but [docker volumes](https://docs.docker.com/storage/volumes/)
with any supported storage driver should also work fine.

```bash
# EITHER using docker compose (which has the bind-mounts predefined):
docker compose up adria-sandbox

# OR using docker run:
docker run --rm \
    --mount type=bind,source="$(pwd)"/input,target=/data/input \
    --mount type=bind,source="$(pwd)"/output,target=/data/output \
    ADRIA.jl/adria-sandbox:latest
```

The `dev.jl` script should run with any files you have provided in your input volume, and will
create and output files in your output volume. The container will be removed once the script completes.

&nbsp;

### Interacting with `adria-sandbox`

If you prefer to work with the sandbox code from an interactive shell, then you will need to override the default entrypoint and command combination to make the container launch your preferred shell at startup.

In this case, you may also like to bind-mount your sandbox and/or ADRIA.jl source code into the container,
so that you can edit it without having to re-build the container each time.

e.g. to launch the sandbox container as a development environment with an interactive bash shell and use:

```bash
# EITHER using docker compose:
docker compose run --rm --entrypoint /bin/bash adria-sandbox

# OR using just docker run:
docker run --rm --interactive --tty --entrypoint /bin/bash \
    --mount type=bind,source="$(pwd)"/input,target=/data/input \
    --mount type=bind,source="$(pwd)"/output,target=/data/output \
    --mount type=bind,source="$(pwd)"/sandbox,target=/opt/adria-sandbox/src \
    --mount type=bind,source="$(pwd)"/ADRIA.jl,target=/usr/local/src/adria \
    ADRIA.jl/adria-sandbox:latest
```

In both cases, the initial working directory will be the installation-location for the sandbox application's source code,
and you can invoke any command you like from that shell prompt, e.g:

```bash
julia --project=@. dev.jl
```

**Warning:** For `julia` commands, you will probably need to use the `--project=@.` argument. This tells Julia that it's working environment is based in a parent directory of the one the sandbox source code is installed to, which is where the build
configured all the precompiled dependencies. If you omit this, the pre-installed packages may not all be available.
