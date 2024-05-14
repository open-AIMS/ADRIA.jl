using DataEnvelopmentAnalysis: DataEnvelopmentAnalysis as DEA
using DataFrames, YAXARrays
using StableRNGs

function deployment_cost(rs::ResultSet)::YAXArray
    scenarios = rs.inputs
    site_data = rs.site_data

    # Get total cost for no. of corals deployed in each scenario
    prod_cost = prodcution_cost(scenarios, site_data)
    # Get logistics cost for corals deployed in each scenarios
    logistics_cost = logistics_cost(scenarios, site_data)

    return YAXArray((Dim{scenarios}(1:size(scenarios, 1)),), prod_cost .+ logistics_cost)
end
function data_envelopment_analysis(
    rs::ResultSet, metrics...; dea_function::Function=dea, rts::Symbol=:VRS,
    orient::Symbol=:Output
)::Tuple{AbstractDEAModel,AbstractArray}
    cost = deployment_cost(rs)
    X = metrics[keys(metrics)[1]]
    for metric in keys(metrics)[2:end]
        vcat!(X, Array(metrics[metric]))
    end

    result = dea_function(Array(cost), X; orient=orient, rts=rts)
    return result, X
end
