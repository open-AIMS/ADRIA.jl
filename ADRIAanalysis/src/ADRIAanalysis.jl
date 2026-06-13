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
    PrecompileTools,
    Random,
    StaticArrays,
    Statistics,
    StatsBase,
    YAXArrays

import ADRIA.analysis: col_normalize, scenario_types, scenario_rcps

include("analysis/clustering.jl")
include("analysis/data_envelopment.jl")
include("analysis/feature_set.jl")
include("analysis/rule_extraction.jl")
include("analysis/pareto.jl")
include("analysis/screening.jl")
include("analysis/scenario.jl")
include("sensitivity/sensitivity.jl")

@compile_workload begin
    _n = 12
    _X = rand(Float64, _n, 3)
    _y = rand(Float64, _n)
    _Y2 = rand(Float64, _n, 4)
    _df = DataFrame(_X, [:a, :b, :c])
    _ci = [1, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2]
    _bv = BitVector(isodd.(_ci))
    _ms = DataFrame(;
        fieldname=[:a, :b, :c],
        ptype=fill("continuous", 3),
        lower_bound=zeros(3),
        upper_bound=ones(3)
    )
    _df_rcp = insertcols(copy(_df), :RCP => repeat([45, 60], _n ÷ 2))

    scenario_clusters(_bv)
    scenario_clusters(_ci)

    screen_scenarios(_y, x -> x .> 0.5)
    target_clusters(_ci, _X)

    pawn(_X, _y, ["a", "b", "c"]; S=2)
    pawn(_df, _y; S=2)
    tsa(_df, _Y2)
    rsa(_df, _y, _ms; S=2)
    outcome_map(_df, _y, _bv, [:a, :b, :c], _ms; S=2, n_boot=2)

    find_pareto_optimal(_df_rcp, hcat(_y, rand(_n)), [45, 60])
end

end
