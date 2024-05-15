
function ADRIA.viz.data_envelopment_analysis(
    rs::ResultSet, costs::Vector{Float64}, X::Vector{Float64}, efficiencies::Array{Float64};
    axis_opts=Dict(),
    opts=Dict()
)
    f = Figure()
    g = Axis(f[1, 1]; axis_opts...)
    return ADRIA.viz.data_envelopment_analysis(g, rs, costs, X, efficiencies; opts=opts)
end
function ADRIA.viz.data_envelopment_analysis(g::Union{GridLayout,GridPosition},
    rs::ResultSet, costs::Vector{Float64}, X::Vector{Float64}, efficiencies::Array{Float64};
    opts=Dict()
)
    return ADRIA.viz.data_envelopment_analysis(
        g, costs, X, efficiencies; opts=opts
    )
end
function ADRIA.viz.data_envelopment_analysis(g::Union{GridLayout,GridPosition},
    Y::Vector{Float64}, X::Vector{Float64},
    efficiencies::Vector{Float64}; opts=Dict())
    line_color = get(opts, :line_color, :red)
    data_color = get(opts, :data_color, :black)
    frontier_name = get(opts, :frontier_name, "Best practice frontier")
    data_name = get(opts, :data_name, "Scenario data cloud")

    # Find points on best practice frontier
    best_practice_scens = findall(efficiencies .== 1.0)

    # Plot frontier and data cloud
    frontier = lines!(
        g, Y[best_practice_scens], X[best_practice_scens]; color=line_color
    )
    data = scatter!(g, Y, X; color=data_color)
    Legend(g, [frontier, data], [frontier_name, data_name])

    return g
end
