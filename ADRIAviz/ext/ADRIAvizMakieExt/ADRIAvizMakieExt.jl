module ADRIAvizMakieExt

using Base.Iterators
using Reexport

@reexport using GeoMakie

using Statistics, Distributions, Random
using DataFrames

using ImageIO, GeoInterface

import GeoMakie.GeoJSON.FeatureCollection as FC

import ADRIA.FileIO, ADRIA.GFT
using ADRIA
using ADRIA:
    ResultSet

using ADRIAviz
using ADRIAviz:
    COLORS, labels,
    _get_scenario_groups,
    _scenario_types, _scenario_rcps, _scenario_clusters,
    relative_sensitivities, outcome_probability,
    outcome_title, outcome_label, set_plot_opts!,
    OPT_TYPE, DEFAULT_OPT_TYPE,
    _time_labels, _calc_gridsize, timesteps

import ADRIA: timesteps as AD_timesteps

Random.seed!(101)

include("./plotting.jl")
include("./theme.jl")
include("./viz/viz.jl")

end  # module
