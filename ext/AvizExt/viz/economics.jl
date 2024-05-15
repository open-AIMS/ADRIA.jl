
function ADRIA.viz.data_envelopment_analysis(
    rs::ResultSet, X::Array{Float64}, efficiencies::Array{Float64}; axis_opts=Dict(),
    opts=Dict()
)
    f = Figure()
    g = Axis(f[1, 1]; axis_opts...)
    return ADRIA.viz.data_envelopment_analysis(g, rs, X, efficiencies; opts=opts)
end
function ADRIA.viz.data_envelopment_analysis(g::Union{GridLayout,GridPosition},
    rs::ResultSet, X::Array{Float64}, efficiencies::Array{Float64}; opts=Dict()
)
    return ADRIA.viz.data_envelopment_analysis(
        g, X[1, :], X[2, :], efficiencies; opts=opts
    )
end
function ADRIA.viz.data_envelopment_analysis(g::Union{GridLayout,GridPosition},
    metric_1::AbstractArray, metric_2::AbstractArray,
    efficiencies::Array; opts=Dict())
    line_color = get(opts, :line_color, :red)
    data_color = get(opts, :data_color, :black)
    frontier_name = get(opts, :frontier_name, "Best practice frontier")
    data_name = get(opts, :data_name, "Scenario data cloud")

    # Find points on best practice frontier
    best_practice_scens = findall(efficiencies .== 1.0)

    # Plot frontier and data cloud
    frontier = lines!(
        g, metric_1[best_practice_scens], metric_2[best_practice_scens]; color=line_color
    )
    data = scatter!(g, metric_1, metric_2; color=data_color)
    Legend(g, [frontier, data], [frontier_name, data_name])

    return g
end
