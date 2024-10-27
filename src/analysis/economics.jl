module economics

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
    DEAResult(CRS_eff::Vector{Float64}, VRS_eff::Vector{Float64}, FDH_eff::Vector{Float64},
        CRS_peers::DEA.DEAPeers, VRS_peers::DEA.DEAPeers, FDH_peers::DEA.DEAPeers,
        X::Matrix{Float64}, Y::Vector{Float64})::DEAResult

Constructor for DEAResult type.

# Arguments
- `CRS_eff` : efficiencies from CRS DEA analysis.
- `VRS_eff` : efficiencies from VRS DEA analysis.
- `FDH_eff` : efficiencies from FDH DEA analysis.
- `CRS_peers` : peers indices from CRS DEA analysis.
- `VRS_peers` : peers indices from VRS DEA analysis.
- `FDH_peers` : peers indices from FDH DEA analysis.
- `X` : inputs.
- `Y` : outputs.
"""
function DEAResult(CRS_eff::Vector{Float64}, VRS_eff::Vector{Float64},
    FDH_eff::Vector{Float64}, CRS_peers::DEA.DEAPeers, VRS_peers::DEA.DEAPeers,
    FDH_peers::DEA.DEAPeers, X::Matrix{Float64}, Y::Vector{Float64}
)::DEAResult
    return DEAResult(1 ./ CRS_eff,
        1 ./ VRS_eff,
        1 ./ FDH_eff,
        CRS_peers,
        VRS_peers,
        FDH_peers,
        X,
        Y)
end

"""
    data_envelopment_analysis(X::YAXArray, Y::YAXArray; orient::Symbol=:Output,
        dea_model::Function=DEA.deabigdata)::DEAResult
    data_envelopment_analysis(X::YAXArray, metrics...; orient::Symbol=:Output,
        dea_model::Function=DEA.deabigdata)::DEAResult

Performs output-oriented (default) Data Envelopment Analysis (DEA) given inputs X and output
metrics Y. DEA is used to measure the performance of entities (scenarios), where inputs are
converted to outputs via some process. Each scenario's "efficiency score" is calculated
relative to an "efficiency fromtier", a region representing scenarios for which outputs
cannot be further increased by changing inputs (scenario settings). Scenarios on the
frontier serve as "benchmarks" or "peers", associated with best practice restoration
scenarios. Scenarios with efficiencies not equal to 1 can be improved to be more efficient.


# Arguments
- `X` : Model inputs for each scenario (usually costs).
- `Y` : Model outputs for each scenario (metrics such as tac, rci etc.).
- `orient` : Orientation of the analysis. Can be output oriented (`orient=:Output`), which
    seeks to maximise outputs for a given level of input, or input oriented (`orient=:Input`),
    which seeks to minimise an input for a given level of output.

# Returns
DEAResult, which summarizes inputs, outputs, efficiencies and peers for each scenario.

# Examples
```julia
dom = ADRIA.load_domain("example_domain")
scens = ADRIA.sample(dom, 128)
rs = ADRIA.run_scenarios(dom, scens, "45")

# Get cost of deploying corals in each scenario, with user-specified function
cost = cost_function(scens)

# Get mean coral cover and shelter volume for each scenario
s_tac = dropdims(
    mean(ADRIA.metrics.scenario_total_cover(rs); dims=:timesteps); dims=:timesteps
)
s_sv::Vector{Float64} =
    dropdims(
        mean(mean(ADRIA.metrics.absolute_shelter_volume(rs); dims=:timesteps); dims=:sites);
        dims=(:timesteps, :sites)
    )

# Do output oriented DEA analysis seeking to maximise cover and shelter volume for minimum
# deployment cost.
DEA_scens = ADRIA.economics.data_envelopment_analysis(cost, s_tac, s_sv)

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
function data_envelopment_analysis(rs::ResultSet, input_function::Function,
    metrics...; orient::Symbol=:Output, dea_model::Function=DEA.deabigdata
)::DEAResult
    X = input_function(rs.inputs)
    return data_envelopment_analysis(
        X, metrics; orient=orient, dea_model=dea_model
    )
end
function data_envelopment_analysis(rs::ResultSet, Y::YAXArray, input_function::Function;
    orient::Symbol=:Output, dea_model::Function=DEA.deabigdata
)::DEAResult
    X = input_function(rs.inputs)
    return data_envelopment_analysis(
        X, Y; orient=orient, dea_model=dea_model
    )
end
function data_envelopment_analysis(
    X::YAXArray, metrics...; orient::Symbol=:Output, dea_model::Function=DEA.deabigdata
)::DEAResult
    Y = Array(hcat(metrics...))
    return data_envelopment_analysis(
        Array(X), Y; orient=orient, dea_model=dea_model
    )
end
function data_envelopment_analysis(
    X::YAXArray, Y::YAXArray; orient::Symbol=:Output, dea_model::Function=DEA.deabigdata
)::DEAResult
    return data_envelopment_analysis(
        Array(X), Array(Y); orient=orient, dea_model=dea_model
    )
end
function data_envelopment_analysis(
    X::Vector{Float64}, Y::Matrix{Float64}; orient::Symbol=:Output,
    dea_model::Function=DEA.deabigdata
)::DEAResult
    result_CRS = dea_model(X, Y; orient=orient, rts=:CRS)
    result_VRS = dea_model(X, Y; orient=orient, rts=:VRS)
    result_FDH = dea_model(X, Y; orient=orient, rts=:FDH)

    CRS_peers = DEA.peers(result_CRS)
    VRS_peers = DEA.peers(result_VRS)
    FDH_peers = DEA.peers(result_FDH)

    return DEAResult(
        result_CRS.eff,
        result_VRS.eff,
        result_FDH.eff,
        CRS_peers,
        VRS_peers,
        FDH_peers,
        X,
        Y
    )
end

end
