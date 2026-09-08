# # Loading Results
#
# ## Loading ADRIA Results
#
# Results produced by `ADRIA.run_scenarios` are stored on disk and can be reloaded at any
# time using their path:
#
# ```julia
# rs = ADRIA.load_results("path/to/result_set")
# ```
#
# The returned `rs` is an `ADRIAResultSet` that gives access to everything needed for
# analysis and visualization.
#
# ### Key fields
#
# | Field | Description |
# |-------|-------------|
# | `rs.name` | Name of the result set (derived from the domain name and run time) |
# | `rs.RCP` | RCP scenario string (e.g. `"45"`) |
# | `rs.inputs` | `DataFrame` of scenario inputs used for the run |
# | `rs.model_spec` | `DataFrame` describing all model parameters and their bounds |
# | `rs.outcomes` | `Dict` of named outcome arrays (coral cover, shelter volume, etc.) |
# | `rs.ranks` | Location ranking log; dims `(timesteps, locations, intervention, scenarios)` |
# | `rs.seed_log` | Seeding deployment log; dims `(timesteps, locations, species, scenarios)` |
# | `rs.shading_log` | Fogging/shading log; dims `(timesteps, locations, intervention, scenarios)` where `intervention` is `["fog", "shade"]` |
# | `rs.mc_log` | Moving-coral (assisted migration) log |
# | `rs.coral_dhw_tol_log` | Per-location DHW tolerance trajectories (only populated when `log_dhw_tols = true` in `config.toml`) |
# | `rs.coral_cover_log` | Raw coral cover for all size classes (only populated when `log_cover = true` in `config.toml`) |
# | `rs.loc_ids` | Location identifiers |
# | `rs.loc_area` | Location areas (m2) |
# | `rs.loc_centroids` | Location centroid coordinates |
# | `rs.loc_data` | `DataFrame` with spatial attributes for each location |
# | `rs.dhw_stats` | Summary statistics for the DHW projections used |
# | `rs.wave_stats` | Summary statistics for the wave stress projections used |
# | `rs.connectivity_data` | Connectivity matrix data |
#
# ### Accessing outcomes
#
# Individual outcome arrays can be extracted directly from `rs.outcomes` or via the
# `ADRIA.metrics.*` functions (preferred):
#
# ```julia
# # Via the metrics API (recommended)
# tac = ADRIA.metrics.total_absolute_cover(rs)
# rsv = ADRIA.metrics.relative_shelter_volume(rs)
#
# # Directly from the outcomes dict
# rs.outcomes[:relative_cover]
# ```
#
# See [Running scenarios](@ref) and the [Metrics](@ref) page for more detail on
# available metrics and result set properties.
#
# ## Loading ReefModEngine Results
#
# Results from ReefModEngine.jl can be loaded with the `load_results` function.
#
# ```julia
# rs = ADRIA.load_results(RMEResultSet, "<path to data dir>")
# ```
#
# Expected data directory structure:
#
# ```bash
# data_dir
# |
# +---con_bin
# |       CONNECT_ACRO_2010_11.bin
# |       CONNECT_ACRO_2011_12.bin
# |       ...
# |
# +---id
# |       id_list_2023_03_30.csv
# |
# +---region
# |       reefmod_gbr.gpkg
# |
# +---results
#         results.nc
#         scenarios.csv
# ```
#
# To reduce duplication of geospatial and connectivity data, the data directory and results
# directory can be supplied separately to avoid keeping copies for each result set analysed.
#
# ```julia
# rs = ADRIA.load_results(RMEResultSet, "<path to data dir>", "<path to results dir>")
# ```
#
# ## Loading C~scape Results
#
# Results from C~scape can be loaded with the `load_results` function.
#
# The first argument is always the **C~scape data package directory** — the folder holding
# `ScenarioID.csv` plus the `connectivity/`, `site_data/` and `initial_cover/`
# subdirectories (see the tree below).
#
# ```julia
# rs = ADRIA.load_results(
#     CScapeResultSet, "<path to C~scape data package>";
#     result_dir="<path to result NetCDF directory>",
#     result_files=["NetCDF_Scn_140001.nc", "NetCDF_Scn_142162.nc"],
#     show_progress=true
# )
# ```
#
# All keyword arguments are optional:
#
# - Omit `result_dir` and the NetCDFs are read from the data package's own `results/`
#   subdirectory. Set it to point at NetCDFs kept outside the data package.
# - Omit `result_files` and every `NetCDF_Scn_*` file in `result_dir` is loaded. Set it to a
#   list of NetCDF paths to load only those; `result_dir` is then ignored.
# - `show_progress` (default `true`) toggles the progress bar shown while outcomes are
#   computed.
#
# Expected C~scape data package structure (the directory passed as the first argument):
#
# ```bash
# cscape_data_package
# |   ScenarioID.csv
# |
# +---connectivity
# |       connectivity.csv
# |
# +---site_data
# |       geospatial_data.gpkg
# |
# +---initial_cover
# |       initial_cover.csv
# |
# +---results (optional)
#         NetCDF_Scn_140001.nc
#         NetCDF_Scn_140002.nc
#         ...
# ```
#
# ### C~scape data package
#
# Most of these files can be sourced from the RRAP data store (published dataset names in
# italics below). They are usually distributed as R `.Rdata` objects or `write.table` text
# and need light reformatting into the CSV / GeoPackage layout shown above.
#
# | File | Contents | Data store source |
# |------|----------|-------------------|
# | `ScenarioID.csv` | One row per scenario: input parameters, intervention settings and the datasets each run used. The `ID` column matches the `NetCDF_Scn_<ID>` result files. | not yet published |
# | `connectivity/connectivity.csv` | Larval connectivity matrix between locations, with `reef_siteid` row and column labels. | *Spatial inputs - Moore cluster 2022 v2* (`MEAN_all_Connectivity_MooreReef_cluster_221019.Rdata`) |
# | `site_data/*.gpkg` | Location polygons and their spatial attributes (`reef_siteid`, `k`, `area`, depth, ...). The first `.gpkg` found in the folder is used. | *Spatial inputs - Moore cluster 2022 v2* (`MooreReefCluster_Polygon_Geometry.Rdata`) |
# | `initial_cover/initial_cover.csv` | Initial coral cover per location and functional group. | *Coral Cover Initialisation data inputs - C~scape - Counterfactuals Mar 2024* |
#
# ### C~scape model outputs
#
# The results directory holds one NetCDF per scenario. Files must contain the
# `NetCDF_Scn_<ID>` prefix to be discovered automatically.
#
# The full model output set is large (~100 GB). RRAP M&DS publishes instructions for
# downloading it via the AWS CLI on the *Model Outputs* data store page (download tab).
#
# ### Accessing C~scape outcomes
#
# Only relative cover is loaded automatically. All other outcomes are computed on demand
# via the `ADRIA.metrics.*` functions and cached in `rs.outcomes`:
#
# ```julia
# settlers = ADRIA.metrics.total_settlers(rs)
# rs.outcomes[:total_settlers]  # now cached
# ```
