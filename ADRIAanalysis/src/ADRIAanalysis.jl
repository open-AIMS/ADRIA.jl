module ADRIAanalysis

using ADRIA
using ADRIA: ResultSet, ZeroDataCube, DataCube, n_locations, model_spec
using ADRIA.metrics: nds

using
    Bootstrap,
    Clustering,
    DataEnvelopmentAnalysis,
    DataFrames,
    Distances,
    Distributions,
    HypothesisTests,
    Random,
    StaticArrays,
    Statistics,
    StatsBase,
    YAXArrays

import ADRIA.analysis: col_normalize, scenario_types, scenario_rcps

include("analysis/clustering.jl")
include("analysis/rule_extraction.jl")
include("analysis/pareto.jl")
include("analysis/screening.jl")
include("analysis/scenario.jl")

end
