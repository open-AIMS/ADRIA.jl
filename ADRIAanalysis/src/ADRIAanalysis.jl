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
    MLJ,
    Random,
    SIRUS,
    StableRNGs,
    StaticArrays,
    Statistics,
    StatsBase,
    YAXArrays

import ADRIA.analysis: col_normalize, scenario_types, scenario_rcps

include("analysis/clustering.jl")
include("analysis/rule_extraction.jl")
include("analysis/data_envelopment.jl")
include("analysis/feature_set.jl")
include("analysis/intervention.jl")
include("analysis/pareto.jl")
include("analysis/screening.jl")
include("analysis/scenario.jl")
include("sensitivity/sensitivity.jl")

end
