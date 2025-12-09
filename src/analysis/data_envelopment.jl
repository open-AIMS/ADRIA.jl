using DataEnvelopmentAnalysis: DataEnvelopmentAnalysis as DEA
using ADRIA: ResultSet
using DataFrames, YAXArrays

struct DEAResult{V,V2,V3}
    crs_vals::V # Inverse efficiencies using constant returns to scale.
    vrs_vals::V # Inverse efficiencies using variable returns to scale.
    fdh_vals::V # Inverse efficiencies using free disposability hull (non-convexity assumption).
    crs_peers::V2 # Scenarios on the efficiency frontier.
    vrs_peers::V2 # Scenarios on the efficiency frontier.
    fdh_peers::V2 # Scenarios on the efficiency frontier.
    X::V # Inputs
    Y::V3 # Outputs
end

"""
    data_envelopment_analysis(rs::ResultSet, input_function::Function, metrics...; orient::Symbol=:Output, dea_model::Function=DEA.deabigdata)::DEAResult
    data_envelopment_analysis(rs::ResultSet, Y::AbstractArray, input_function::Function; orient::Symbol=:Output, dea_model::Function=DEA.deabigdata)::DEAResult
    data_envelopment_analysis(X::YAXArray, metrics...; orient::Symbol=:Output, dea_model::Function=DEA.deabigdata)::DEAResult
    data_envelopment_analysis(X::YAXArray, Y::AbstractArray; orient::Symbol=:Output, dea_model::Function=DEA.deabigdata)::DEAResult
    data_envelopment_analysis(X::Union{Vector{Float64}, Matrix{Float64}}, Y::Matrix{Float64}; dea_model::Function=DEA.deabigdata)::DEAResult

Performs output-oriented (default) Data Envelopment Analysis (DEA) given inputs X and output
metrics Y. DEA is used to measure the performance of entities (scenarios), where inputs are
converted to outputs via some process. Each scenario's "efficiency score" is calculated
relative to an "efficiency fromtier", a region representing scenarios for which outputs
cannot be further increased by changing inputs (scenario settings). Scenarios on the
frontier serve as "benchmarks" or "peers", associated with best practice restoration
scenarios. Scenarios with efficiencies not equal to 1 can be improved to be more efficient.

# Arguments
- `rs` : ADRIA ResultSet
- `input_function` : function which calculates an input for each scenario (e.g. cost, effort) given the scenario
    dataframe as input.
- `X` : Model inputs for each scenario (usually costs).
- `Y` : Model outputs for each scenario (metrics such as tac, rci etc.).
- `orient` : Orientation of the analysis. Can be output oriented (`orient=:Output`), which
    seeks to maximise outputs for a given level of input, or input oriented (`orient=:Input`),
    which seeks to minimise an input for a given level of output.
- `dea_model` : model to use to calculate DEA frontier (see https://javierbarbero.github.io/DataEnvelopmentAnalysis.jl/stable/)

# Returns
DEAResult, which summarizes inputs, outputs, efficiencies and peers for each scenario.

# Examples
```julia
dom = ADRIA.load_domain("example_domain", "<RCP>")
scens = ADRIA.sample(dom, 128)
rs = ADRIA.run_scenarios(dom, scens, "45")

# Get cost of deploying corals in each scenario, with user-specified function
cost = cost_function(scens)

# Get mean coral cover and shelter volume for each scenario
s_tac = dropdims(
    mean(ADRIA.metrics.scenario_total_cover(rs); dims=:timesteps); dims=:timesteps
)
s_sv = dropdims(
    mean(ADRIA.metrics.scenario_shelter_volume(rs); dims=:timesteps); dims=:timesteps
    )

# Do output oriented DEA analysis seeking to maximise cover and shelter volume for minimum
# deployment cost.analysis.data_envelopment_analysis(cost, s_tac, s_sv)

```
# References
1. Huguenin, J-M., 2012.
   Data Envelopment Analysis (DEA): A pedagogical guide for decision makers in the public
   sector. https://api.semanticscholar.org/CorpusID:188267263
2. Pascoe, S., Cannard, T., Dowling, N.A., et. al., 2023.
   Use of Data Envelopment Analysis (DEA) to assess management alternatives in the presence
   of multiple objectives.
   Marine Policy, 148, 105444.
   https://doi.org/10.1016/j.marpol.2022.105444
3. Pascoe, S., 2024.
   On the use of Data Envelopment Analysis for Multi-Criteria Decision Analysis.
   Algorithms, 17:89.
   https://doi.org/10.3390/a17030089
"""
function data_envelopment_analysis(
    rs::ResultSet,
    input_function::Function,
    metrics...;
    orient::Symbol=:Output,
    dea_model::Function=DEA.deabigdata
)::DEAResult
    return data_envelopment_analysis(
        input_function(rs.inputs), metrics; orient=orient, dea_model=dea_model
    )
end
function data_envelopment_analysis(
    rs::ResultSet,
    Y::AbstractArray,
    input_function::Function;
    orient::Symbol=:Output,
    dea_model::Function=DEA.deabigdata
)::DEAResult
    return data_envelopment_analysis(
        input_function(rs.inputs), Y; orient=orient, dea_model=dea_model
    )
end
function data_envelopment_analysis(
    X::YAXArray,
    metrics...;
    orient::Symbol=:Output,
    dea_model::Function=DEA.deabigdata
)::DEAResult
    Y = Array(hcat(metrics...))
    return data_envelopment_analysis(
        X.data, Y; orient=orient, dea_model=dea_model
    )
end
function data_envelopment_analysis(
    X::YAXArray,
    Y::AbstractArray;
    orient::Symbol=:Output,
    dea_model::Function=DEA.deabigdata
)::DEAResult
    return data_envelopment_analysis(
        X.data, Array(Y); orient=orient, dea_model=dea_model
    )
end
function data_envelopment_analysis(
    X::Union{Vector{Float64},Matrix{Float64}},
    Y::Matrix{Float64};
    orient::Symbol=:Output,
    dea_model::Function=DEA.deabigdata
)::DEAResult
    result_CRS = dea_model(X, Y; orient=orient, rts=:CRS)
    result_VRS = dea_model(X, Y; orient=orient, rts=:VRS)
    result_FDH = dea_model(X, Y; orient=orient, rts=:FDH)

    return DEAResult(
        1 ./ result_CRS.eff,
        1 ./ result_VRS.eff,
        1 ./ result_FDH.eff,
        DEA.peers(result_CRS),
        DEA.peers(result_VRS),
        DEA.peers(result_FDH),
        X,
        Y
    )
end
