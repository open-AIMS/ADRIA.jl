module economics

using CSV
using DataEnvelopmentAnalysis: DataEnvelopmentAnalysis as DEA
using BasicInterpolators
using ADRIA: ResultSet
using DataFrames, YAXArrays

"""
    data_envelopment_analysis(X::YAXArray, Y::YAXArray; rts::Symbol=:VRS,
        orient::Symbol=:Output, dea_model::Function=DEA.deabigdata)
        ::Tuple{DEA.AbstractDEAModel,AbstractArray}
    data_envelopment_analysis(X::YAXArray, metrics...; rts::Symbol=:VRS,
        orient::Symbol=:Output, dea_model::Function=DEA.deabigdata)
        ::Tuple{DEA.AbstractDEAModel,AbstractArray}

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
- `rts` : Returns to scale, can be constant returns to scale (assumes all scenarios are at
    their optimal scale, ie. a multiple would not be more optimal, `rts=:CRS`), or variable
    returns to scale (doesn't assume the optimal scale has been acheived, `rts=:VRS`)
- `orient` : Orientation of the analysis. Can be output oriented (`orient=:Output`), which
    seeks to maximise outputs for a given level of input, or input oriented (`orient=:Input`),
    which seeks to minimise an input for a given level of output.

# Returns
DEA.AbstractDEAModel, which summarizes efficiencies and slack variables (or multipliers
depending on the DEA method) for each scenario. Slack variables or multipliers give
information about the effort required to move an inefficient scenario towards it's peer
efficient scenario. X, a matrix of the outputs which can be used to plot the efficiency
frontier.


# Examples
```julia
dom = ADRIA.load_domain("example_domain")
scens = ADRIA.sample(dom, 128)
rs = ADRIA.run_scenarios(dom, scens, "45")

# Get cost of deploying corals in each scenario
CAD_cost = ADRIA.economics.CAD_cost(scens)

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
DEA_scens = ADRIA.economics.data_envelopment_analysis(CAD_cost, s_tac, s_sv)

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
   On the use of Data Envelopement Analysis for Multi-Criteria Decision Analysis.
   Algorithms, 17:89.
   https://doi.org/10.3390/a17030089

"""
function data_envelopment_analysis(
    X::YAXArray, metrics...; rts::Symbol=:VRS,
    orient::Symbol=:Output, dea_model::Function=DEA.deabigdata
)::Tuple{DEA.AbstractDEAModel,AbstractArray}
    Y = Array(hcat(metrics...))
    return data_envelopment_analysis(
        X, Y; rts=rts, orient=orient, dea_model=dea_model
    )
end
function data_envelopment_analysis(
    rs::ResultSet, cost::Vector{Float64}, metrics...; rts::Symbol=:VRS,
    orient::Symbol=:Output, dea_model::Function=DEA.deabigdata
)::Tuple{DEA.AbstractDEAModel,AbstractArray}
    X = Array(hcat(metrics...))
    return data_envelopment_analysis(
        rs, cost, X; rts=rts, orient=orient, dea_model=dea_model
    )
end
function data_envelopment_analysis(
    rs::ResultSet, cost::Vector{Float64}, X::Matrix{Float64}; rts::Symbol=:VRS,
    orient::Symbol=:Output, dea_model::Function=DEA.deabigdata
)::Tuple{DEA.AbstractDEAModel,AbstractArray}
    result = dea_model(Array(cost), X; orient=orient, rts=rts)
    return result, X
end

"""
    CAD_cost(rs; Reef::String="Moore")::YAXArray
    CAD_cost(scenarios::DataFrame; Reef::String="Moore")::YAXArray

Calculate the cost of coral deployments for a set of scenarios. Based on piecewise linear
    interpolations of cost data collected from `3.5.1 CA Deployment Model.xls`. Assumes the
    ship Cape Ferguson is used and a 28 day deployment window (effects set-up cost only).

# Arguments
- `rs` : ResultSet
- `scenarios` : sampled input scenario specifications.
- `Reef` : Reef to travel to (impacts cost significantly for <10000 devices deployed).
    Currently cost data only available for ["Moore", "Davies", "Swains", "Keppel"], but
    could be used as an approximation for nearby clusters.

# Returns
YAXArray, containing estimated cost for each scenario in `scenarios`.

"""
function CAD_cost(rs; Reef::String="Moore")::YAXArray
    return CAD_cost(rs.inputs; Reef=Reef)
end
function CAD_cost(scenarios::DataFrame; Reef::String="Moore")::YAXArray
    # No. of corals deployed in each scenario
    scen_no_corals =
        scenarios[:, :N_seed_CA] .+ scenarios[:, :N_seed_SM] .+ scenarios[:, :N_seed_TA]

    scen_no_years = scenarios[:, :seed_years]
    deploy_op_cost = CSV.read("deploy_op_cost.csv", DataFrame; header=false)
    deploy_cap_cost = CSV.read("deploy_cap_cost.csv", DataFrame; header=false)
    reef_ind = findfirst(deploy_op_cost[2:end, 1] .== Reef)

    OP_lin = LinearInterpolator(
        Array(deploy_op_cost[1, 2:end]), Array(deploy_op_cost[reef_ind, 2:end]),
        NoBoundaries()
    )
    CAP_lin = LinearInterpolator(
        Array(deploy_cap_cost[1, :]), Array(deploy_cap_cost[2, :]), NoBoundaries()
    )

    return (OP_lin.(scen_no_corals) .* scen_no_years .+ CAP_lin.(scen_no_corals)) ./ (10^6)
end

end
