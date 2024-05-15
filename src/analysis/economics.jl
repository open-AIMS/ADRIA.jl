module economics

using CSV
using DataEnvelopmentAnalysis: DataEnvelopmentAnalysis as DEA
using BasicInterpolators
using ADRIA: ResultSet
using DataFrames, YAXArrays

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

function CAD_cost(scenarios::DataFrame; Reef::String="Moore")
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
