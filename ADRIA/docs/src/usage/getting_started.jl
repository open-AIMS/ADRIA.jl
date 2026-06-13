# # Getting Started
#
# ## Setup
#
# This section outlines how ADRIA may be used to arrive at a select range of possible
# pathways that are robust to possible future conditions.
#
# Create a directory for your project, and start Julia inside that directory:
#
# ```bash
# $ julia --project=. --threads=auto
# ```
#
# The ADRIA ecosystem consists of three packages. Install them all from the package manager:
#
# ```julia
# julia> ]add ADRIA ADRIAviz ADRIAanalysis
# ```
#
# To install the latest development version of all three packages:
#
# ```julia
# using Pkg
# Pkg.add([
#     Pkg.PackageSpec(url="https://github.com/open-AIMS/ADRIA.jl", subdir="ADRIA"),
#     Pkg.PackageSpec(url="https://github.com/open-AIMS/ADRIA.jl", subdir="ADRIAviz"),
#     Pkg.PackageSpec(url="https://github.com/open-AIMS/ADRIA.jl", subdir="ADRIAanalysis"),
# ])
# ```
#
# Update all three packages as new releases are made:
#
# ```julia
# julia> ]up ADRIA ADRIAviz ADRIAanalysis
# ```
#
# Package roles:
#
# - **ADRIA** provides the core simulation engine, MCDA, and base metrics.
# - **ADRIAviz** provides visualization tools. The Plotly backend (recommended) requires
#   `PlotlyBase`; Makie backends (`"WGLMakie"`, `"GLMakie"`, `"CairoMakie"`) require
#   GeoMakie and GraphMakie.
# - **ADRIAanalysis** provides extended analysis functions including `target_clusters`,
#   `screen_scenarios`, `pawn`, `data_envelopment_analysis`, `rules`, and others.
#
# Install the Plotly backend dependencies:
#
# ```julia
# julia> ]add PlotlyBase            # required for Plotly backend
# julia> ]add PlotlyKaleido         # optional, enables static image export
# ```
#
# If desired, create a `config.toml` file inside your project directory.
# The assumed default values are shown below:
#
# ```toml
# [operation]
# threshold = 1e-8      # Result values below this will be set to 0.0 (to save disk space).
# debug = false         # Disable multi-threading to allow error messages to be shown.
# log_dhw_tols = false  # Log per-location coral DHW tolerance trajectories (increases disk usage).
# log_cover = false     # Log raw results (coral cover) for all timesteps, locations, functional groups, size classes and scenarios.
# rng_seed = false      # Set to an integer to be used as RNG seed for all runs.
#
# [results]
# output_dir = "./Outputs"  # Change this to point to where you want to store results.
# ```
#
# This `config.toml` file is specific to your computer and project. It **should not be**
# committed to version control.
#
# !!! tip "Performance"
#     ADRIA uses an on-disk data store to hold results from model runs.
#     Setting `output_dir` to a directory on an SSD (Solid State Drive) will maximize
#     performance.
#
# To setup ADRIA for development, see the [Development setup](@ref) page.
#
# ## Quick Start
#
# A common workflow would be the following.
#
# Start Julia from the project directory with multi-threading enabled:
#
# ```bash
# $ julia --project=. --threads=auto
# ```
#
# Load data for a spatial domain. See [Loading a Domain](@ref) for more details:
#
# ```julia
# using ADRIA
#
# dom = ADRIA.load_domain("path to domain data package directory", "<RCP>")
# ```
#
# Generate scenarios based on available environmental data layers and model parameters. The
# number of scenarios should be a power of two. See [Generating scenarios](@ref) for more
# details:
#
# ```julia
# num_scenarios = 128
# scens = ADRIA.sample(dom, num_scenarios)
# ```
#
# Run sampled scenarios for one or more RCPs. This may take a while:
#
# ```julia
# rcp_45 = "45"
# rs = ADRIA.run_scenarios(dom, scens, rcp_45)
# ```
#
# Or run scenarios across several RCPs:
#
# ```julia
# rcps = ["45", "60", "85"]
# rs = ADRIA.run_scenarios(dom, scens, rcps)
# ```
#
# It is also possible to load previously run scenarios. See [Running scenarios](@ref) for
# more details:
#
# ```julia
# rs = ADRIA.load_results("path to results")
# ```
#
# Extract some metric for analysis (e.g., the total absolute cover for each site and
# timestep):
#
# ```julia
# s_tc = ADRIA.metrics.scenario_total_cover(rs)
# ```
#
# Use ADRIAviz to plot the results. Load the package and activate the Plotly backend before
# calling any `ADRIA.viz.*` function:
#
# ```julia
# using ADRIA, ADRIAviz, PlotlyBase
# ADRIAviz.activate(:plotly)
#
# fig = ADRIA.viz.scenarios(rs, s_tc; axis_opts=Dict(:ylabel => "Absolute Cover"))
# ADRIA.viz.savefig(fig, "scenarios.html")
# ```
#
# For extended analysis, load ADRIAanalysis alongside ADRIA:
#
# ```julia
# using ADRIA, ADRIAanalysis
#
# # Cluster scenarios by temporal behaviour
# tac = ADRIA.metrics.scenario_total_cover(rs)
# clusters = scenario_clusters(tac)
#
# # Sensitivity analysis (PAWN method)
# scens = ADRIA.param_table(rs)
# Si = pawn(scens, vec(mean(tac; dims=(:timesteps, :locations))), ADRIA.component_params(dom))
# ```
#
# See [Analysis](@ref) for further examples of analysis and plots.
#
# ## Shared package depot paths
#
# If multiple Julia processes are used (e.g. running several independent ADRIA instances),
# it is recommended to set a shared `JULIA_DEPOT_PATH`
# so each process does not race against the others to compile packages.
#
# This is typically defined in your `.bashrc` (or equivalent) on Linux, or configured via
# the user environment variable control panel on Windows.
#
# To set a temporary environment variable for a session:
#
# On Linux:
#
# ```shell
# export JULIA_DEPOT_PATH="some_shared_accessible_directory"
# ```
#
# On Windows (Command Prompt):
#
# ```shell
# set JULIA_DEPOT_PATH="some_shared_accessible_directory"
# ```
#
# On Windows (Powershell):
#
# ```powershell
# $Env:JULIA_DEPOT_PATH="some_shared_accessible_directory"
# ```
#
# For VS Code, open settings (Ctrl+,), search for `terminal.integrated.env.windows`, and
# add `"JULIA_DEPOT_PATH": "<path to depot>"`.
#
# See also:
# - [Julia documentation on `JULIA_DEPOT_PATH`](https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH)
# - [Pkg depot documentation](https://pkgdocs.julialang.org/dev/depots/#Platform-specific-configuration)
