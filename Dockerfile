ARG SANDBOX_FROM="adria-dev"

#------------------------------------------------------------------------------
# internal-base build target: julia with OS updates and an empty @adria
# Julia environment prepared for use. NOT intended for standalone use.
#------------------------------------------------------------------------------
# Work with an apt-based Debian 12 (bookworm) OS, but allow
# the Julia platform version to be overridden at build-time
# See https://hub.docker.com/_/julia for valid versions.
ARG JULIA_VERSION="1.10.4"
FROM julia:${JULIA_VERSION}-bookworm AS internal-base

# Record the actual base image used from the FROM command as label in the compiled image
ARG BASE_IMAGE=$BASE_IMAGE
LABEL org.opencontainers.image.base.name=${BASE_IMAGE}

# Update all pre-installed OS packages (to get security updates)
# and add a few extra utilities
RUN --mount=target=/var/lib/apt/lists,type=cache,sharing=locked \
    --mount=target=/var/cache/apt,type=cache,sharing=locked \
    apt-get update \
    && apt-get -y upgrade \
    && apt-get install --no-install-recommends -y \
        git \
        less \
        nano \
    && apt-get clean \
    && apt-get autoremove --purge \
    && rm -rf /var/lib/apt/lists/*

# Tweak the JULIA_DEPOT_PATH setting so that our shared environments will end up
# in a user-agnostic location, not in ~/.julia => /root/.julia which is the default.
# See https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH
# This allows apps derived from this image to drop privileges and run as non-root
# user accounts, but still activate environments configured by this dockerfile.
ENV JULIA_DEPOT_PATH="/usr/local/share/julia"

# Prepare an empty @adria Julia environment for derived images to use
RUN mkdir -p "${JULIA_DEPOT_PATH}" && \
    chmod 0755 "${JULIA_DEPOT_PATH}" && \
    julia -e 'using Pkg; Pkg.activate("adria", shared=true)'

# Ensure the @adria environment is in the load path for Julia, so that apps derived
# from this image can access any packages installed to there.
# (See https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_LOAD_PATH)
ENV JULIA_LOAD_PATH="@:@adria:@v#.#:@stdlib"

# Run Julia commands by default as the container launches.
# Derived applications should override the command.
ENTRYPOINT ["julia"]


#------------------------------------------------------------------------------
# adria-base build target: ADRIA.jl preinstalled as a non-dev package.
# Use `--target=adria-base` on your `docker build command to build *just* this.
#------------------------------------------------------------------------------
FROM internal-base AS adria-base

# Install ADRIA.jl into the @adria shared environment as an unregistered package.
# - Allow the package source and version to be overridden at build-time.
# - Include citation information for ADRIA.jl in the image labels.
ARG ADRIA_REPO="https://github.com/open-AIMS/ADRIA.jl.git" \
    ADRIA_REFSPEC="main"
RUN mkdir -p "${JULIA_DEPOT_PATH}" && \
    chmod 0755 "${JULIA_DEPOT_PATH}" && \
    julia --project=@adria -e "using Pkg; Pkg.add(url=\"${ADRIA_REPO}\", rev=\"${ADRIA_REFSPEC}\"); using ADRIA"
LABEL au.gov.aims.adria.source="${ADRIA_REPO}" \
      au.gov.aims.adria.branch="${ADRIA_REFSPEC}" \
      au.gov.aims.adria.vendor="Australian Institute of Marine Science" \
      au.gov.aims.adria.licenses=MIT

#------------------------------------------------------------------------------
# adria-dev build target: Assumes you have the ADRIA.jl source code
# available in a subdirectory of your docker context named `ADRIA.jl`
# (recommended method is to make it a submodule), and configures that
# as a Julia development package in the @adria shared environment.
#------------------------------------------------------------------------------
FROM internal-base AS adria-dev

ENV ADRIA_ENV_DIR="${JULIA_DEPOT_PATH}/environments/adria" \
    ADRIA_SRC_DIR="/usr/local/src/adria"

# Install the versioned .toml file(s) into the shared adria environment and use
# those to set up the ADRIA source code as a development package in the
# shared @adria environment, pre-installing and precompiling dependencies.
# This should *hugely* speeds up the build when other ADRIA files have changed,
# as this very slow step won't need to be run again.

# NOTE: the Manifest.toml file is NOT currently versioned by the ADRIA.jl
#       respository, but is a pre-requisite for this to work.  Putting it
#       in this adria-docker repository instead is a horrible hack, but works.
WORKDIR "${ADRIA_SRC_DIR}"
COPY ./Project.toml ./Project.toml
# TODO bring this back but for now let it build the manifest
#COPY ./Manifest.toml ./Manifest.toml
RUN julia --project=@adria -e 'using Pkg;  Pkg.instantiate(verbose=true);'

# Install the ADRIA source code and configure it as a development
# package in the @adria shared environment.
# Should be v speedy if the .toml file is unchanged, because all the
# dependencies *should* already be installed.
COPY . .
RUN julia --project=@adria \
          -e  'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.precompile(); using ADRIA;'


#------------------------------------------------------------------------------
# adria-sandbox build target: Test harness for the ADRIA.jl package.
# Use `--target=adria-sandbox` on your `docker build command, or omit the
# `--target` argument entirely to build this one.
# This can use EITHER adria-base OR adria-dev as a base image.
#------------------------------------------------------------------------------
FROM "${SANDBOX_FROM}" AS adria-sandbox

# Ensure the input and output data directories have mountpoints
# prepared and will always be external persistent volumes, even
# if the caller forgets to mount something.
VOLUME /data/input
VOLUME /data/output

# Prepare a working base/project-environment directory for the sandbox code,
# pre-install any additional packages required by the sandbox application
# which are not already dependencies of ADRIA, and force an initial precompile
# of the package cache for all the critical dependencies. May be very slow!
# (This will create a Project.toml and Manifest.toml file in the directory).
WORKDIR /opt/adria-sandbox
RUN julia --project=. -e 'using Pkg; Pkg.add(["WGLMakie", "Bonito", "GeoMakie", "GraphMakie"])' && \
    julia --project=. -e 'using Bonito, WGLMakie, GeoMakie, GraphMakie, ADRIA'

# Install the actual sandbox code into a *subdirectory* of the sandbox project.
# This allows you to safely bind-mount the same files when developing without
# overriding the .toml files that identify package dependencies.
WORKDIR /opt/adria-sandbox/src
COPY sandbox .

# Run the sandbox application at startup, ensuring that Julia picks up
# the parent directory and precompiled packages as the Julia project base.
CMD [ "--project=@.", "dev.jl" ]
