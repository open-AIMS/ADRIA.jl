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
