# ADRIA Docker Setup

This repo contains a multistage Dockerfile for generating containerised applications based on the ADRIA source code in this repository.

- [docker configuration](#docker-configuration)
- [adria-base](#adria-base)
  - [A note about platform build targets](#a-note-about-platform-build-targets)
  - [Published base images](#published-base-images)
  - [Building `adria-base`](#building-adria-base)
  - [Running `adria-base` with a non-interactive Julia command](#running-adria-base-with-a-non-interactive-julia-command)
  - [Running `adria-base` as an interactive Julia shell](#running-adria-base-as-an-interactive-julia-shell)
  - [Running `adria-base` with a non-Julia entrypoint](#running-adria-base-with-a-non-julia-entrypoint)
  - [Deriving an image from `adria-base`](#deriving-an-image-from-adria-base)
- [adria-dev](#adria-dev)
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

### A note about platform build targets

This line in the Dockerfile prompts Julia to precompile the ADRIA package and it's predependencies for a larger range of target platforms. Given that the building of the image could occur on any host, we need to specify a range of platforms so that precompilation (which is CPU architecture dependent) does not need to reoccur when the Docker image is run on the user's system.

The below list came from [here](https://discourse.julialang.org/t/seemingly-unnecessary-precompilation-with-julia-1-9-in-a-docker-container/99892/3) and is confirmed to work for a `x86_64` based architecture.

```Dockerfile
# Try to coerce Julia to build across multiple targets
ENV JULIA_CPU_TARGET=x86_64;haswell;skylake;skylake-avx512;tigerlake
```

If more platforms are needing to be targeted, please lodge a pull request with the output of `uname -m` and other details, and we can expand this list of target platforms.

### Published base images

This repository uses GitHub actions to publish, upon version release, the `adria-base` image - see `.github/workflows/PublishDockerImage.yml`. The below commands such as running and using the ADRIA base image apply to these published images, replacing your locally tagged image with the GitHub Container Registry label e.g. the following will run in an interactive Julia shell, a precompiled terminal with `ADRIA` installed, ready to be used.

```
docker run --tty --interactive ghcr.io/open-aims/adria.jl/adria-base:latest
```

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

See the `docker-compose.yaml` file for an example of how to specify build arguments in docker compose.

### Running `adria-base` with a non-interactive Julia command

e.g. to list the packages installed in the `@adria` shared environment:

```bash
# EITHER using docker compose:
docker compose run --rm adria-base --project=@adria -e 'using Pkg; Pkg.status()'

# OR using just docker run:
docker run --rm ADRIA.jl/adria-base:latest --project=@adria -e 'using Pkg; Pkg.status()'
```

### Running `adria-base` as an interactive Julia shell

To launch an interactive Julia shell:

```bash
# EITHER using docker compose:
docker compose run --rm adria-base

# OR using just docker run:
docker run --rm --interactive --tty ADRIA.jl/adria-base:latest
```

In both cases, type `CTRL-d` to exit the shell and stop the container.

### Running `adria-base` with a non-Julia entrypoint

If you want to use this image to run something _other_ than a Julia command, you can specify an alternate entrypoint at runtime as well as an alternate command. e.g. to launch an interactive `bash` shell in the container for checking filesystem permissions or similar:

```bash
# EITHER using docker compose:
docker compose run --rm --entrypoint /bin/bash adria-base

# OR using just docker run:
docker run --rm --interactive --tty --entrypoint /bin/bash ADRIA.jl/adria-base:latest
```

### Deriving an image from `adria-base`

To make a derived ADRIA application:

- Use `FROM ADRIA.jl/adria-base:latest` (or a related tag)
- Include a `CMD` line in your Dockerfile that provides appropriate arguments to the `julia` command line to invoke your application.

The section of the `Dockerfile` that defines the `adria-sandbox` target described below might be useful inspiration.

---

## adria-dev

The `adria-dev` image variant is an _alternative_ to the `adria-base` image, not a derived application.
Instead of installing [ADRIA.jl](https://github.com/open-AIMS/ADRIA.jl) as a normal package, it looks for
the ADRIA.jl source code in a local subdirectory, and installs that as a [Julia development package](https://pkgdocs.julialang.org/v1/api/#Pkg.develop).

This allows you to use the `adria-dev` container as an `ADRIA.jl` development environment: you can run tests,
bind-mount and edit the code, re-resolve dependencies and all sorts of other useful things without needing
to a native installation of Julia.

### Building `adria-dev`

From within the root of the `ADRIA.jl` project, you can build the `adria-dev` image like so:

```bash
# EITHER using docker compose:
docker compose build adria-dev

# OR using just `docker build`:
docker build --target "adria-dev" --tag ADRIA.jl/adria-dev:latest .
```

The same `JULIA_VERSION` build argument that works with `adria-base` will also work with `adria-dev`.

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

### Running `adria-sandbox`

The `adria-sandbox` image is configured to automatically run the `dev.jl` script when a container
made from that image runs.

For that script to be useful, though, you must first mount two filesystem locations to the container:

- A directory containing your input data files should be mounted at `/data/input`
- A directory where the sandbox application can create output files should be mounted at `/data/output`

This documentation and the `docker-compose.yaml` file in the repository demonstrate
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
