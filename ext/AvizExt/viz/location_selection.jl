using YAXArrays
using ADRIA: ResultSet

"""
    ADRIA.viz.ranks_to_frequencies!(g::Union{GridLayout,GridPosition},rs::ResultSet,
        frequencies::YAXArray,rank_ids::Vector{Int64};opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())
    ADRIA.viz.ranks_to_frequencies!(g::Union{GridLayout,GridPosition},rs::ResultSet,
        frequencies::YAXArray,rank_id::Int64;opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(:color_map => :CMRmap),axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())
    ADRIA.viz.ranks_to_frequencies(rs::ResultSet,frequencies::YAXArray,rank_ids::Union{Int64,Vector{Int64}};
        opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

Plot a spatial map of location selection frequencies.

# Arguments
- `g` : Figure GridPosition or GridLayout.
- `rs` : Result set.
- `frequencies` : Set of frequencies for each rank over a set of scenarios and/or timesteps.
    As calculated using`ranks_to_frequencies`.
- `rank_id`/`rank_ids` : Rank or set of ranks to plot frequency maps for. E.g. 1, [1,2,3].
- `opts` : Aviz options
    - `colorbar_label`, label for colorbar. Defaults to "Relative Cover".
    - `color_map`, preferred colormap for plotting heatmaps.
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `fig_opts` : Additional options to pass to adjust Figure creation
  See: https://docs.makie.org/v0.19/api/index.html#Figure

# Returns
Figure
"""
function ADRIA.viz.ranks_to_frequencies!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    frequencies::YAXArray,
    rank_ids::Vector{Int64};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    sym_rank_ids = Symbol.(rank_ids)
    rank_groups = Dict(rank_grp => rank_grp .== sym_rank_ids for rank_grp in sym_rank_ids)

    if :colormap in keys(opts)
        @assert opts[:color_map] isa Dict
        all_colormaps = opts[:color_map]
    else
        alpha_vals = alphas(rank_groups)
        all_colormaps = _default_colormap(rank_groups, alpha_vals)
    end

    opts[:color_map] = all_colormaps[sym_rank_ids[1]]
    geodata = _get_geoms(rs)
    legend_els = Vector{PolyElement}(undef, length(rank_ids))
    legend_labels = Vector{String}(undef, length(rank_ids))
    opts[:show_colorbar] = get(opts, :show_colorbar, false)

    ADRIA.viz.map!(
        g,
        rs,
        frequencies[ranks=rank_ids[1]];
        opts=opts,
        axis_opts=axis_opts
    )
    legend_els[1] = PolyElement(;
        color=all_colormaps[Symbol(rank_ids[1])][2], strokecolor=:grey, strokewidth=1
    )
    legend_labels[1] = string("Rank ", string(rank_ids[1]))
    ax = content(g[1, 1])  # get GeoAxis

    for rr in rank_ids[2:end]
        poly!(
            ax,
            geodata;
            color=collect(frequencies[ranks=rr]),
            colormap=all_colormaps[Symbol(rr)],
            strokecolor=:grey,
            strokewidth=0.5,
            linestyle=:solid,
            overdraw=true
        )
        legend_els[rr] = PolyElement(;
            color=all_colormaps[Symbol(rr)][2], strokecolor=:grey, strokewidth=1
        )
        legend_labels[rr] = string("Rank ", string(rr))
    end
    Legend(g[1, 2], legend_els, legend_labels; patchsize=(35, 35), rowgap=10)
    return g
end
function ADRIA.viz.ranks_to_frequencies!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    frequencies::YAXArray,
    rank_id::Int64;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    opts[:colorbar_label] = get(opts, :colorbar_label, "Selection frequency")
    opts[:color_map] = get(
        opts,
        :color_map,
        [RGBA{Float32}(1.0, 1.0, 1.0, 1.0), RGBA{Float32}(0.00784314, 0.243137, 1.0, 1.0)]
    )

    return ADRIA.viz.map!(
        g,
        rs,
        frequencies[ranks=rank_id].data;
        opts=opts,
        axis_opts=axis_opts
    )
end
function ADRIA.viz.ranks_to_frequencies(
    rs::ResultSet,
    frequencies::YAXArray,
    rank_ids::Union{Int64,Vector{Int64}};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.ranks_to_frequencies!(
        g,
        rs,
        frequencies,
        rank_ids;
        opts=opts,
        axis_opts=axis_opts
    )

    return f
end

"""
    ADRIA.viz.selection_criteria_map(rs::Union{Domain,ResultSet}, decision_matrix::YAXArray, scores::Vector{Float64}; criteria::Vector{Symbol}=Array(decision_matrix.criteria), opts::OPT_TYPE=DEFAULT_OPT_TYPE(), fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()), axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE()))
    ADRIA.viz.selection_criteria_map!(g::Union{GridLayout,GridPosition}, rs::Union{Domain,ResultSet}, decision_matrix::YAXArray, scores::Vector{Float64}; criteria::Vector{Symbol}=Array(decision_matrix.criteria), opts::OPT_TYPE=DEFAULT_OPT_TYPE(), axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE()))

Plot several spatial maps showing spatial distribution of location selection criteria, with the final map being the
aggregate location selection score.

# Arguments
- `g` : Figure GridPosition or GridLayout.
- `rs` : Result set or Domain.
- `decision_matrix` : Decision matrix of size (n_locs*n_criteria)
- `scores` : Aggregate criteria scores calculated from `decision_matrix`.
- `criteria` : Names of criteria used, default uses decision matrix YAXArray dimension names.
- `opts` : Aviz options
    - `colorbar_label`, label for colorbar. Defaults to "Relative Cover".
    - `color_map`, preferred colormap for plotting heatmaps.
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `fig_opts` : Additional options to pass to adjust Figure creation
  See: https://docs.makie.org/v0.19/api/index.html#Figure

# Returns
Figure
"""
function ADRIA.viz.selection_criteria_map(
    rs::Union{Domain,ResultSet},
    decision_matrix::YAXArray,
    scores::Vector{Float64};
    criteria::Vector{Symbol}=Array(decision_matrix.criteria),
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.selection_criteria_map!(
        g, rs, decision_matrix, scores; criteria=criteria, opts=opts, axis_opts=axis_opts
    )
    return f
end
function ADRIA.viz.selection_criteria_map!(
    g::Union{GridLayout,GridPosition},
    rs::Union{Domain,ResultSet},
    decision_matrix::YAXArray,
    scores::Vector{Float64};
    criteria::Vector{Symbol}=Array(decision_matrix.criteria),
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)
    if length(rs.loc_data.site_id) != size(decision_matrix, 1)
        error("Only unfiltered decision matrices can be plotted.")
    end

    outputs_matrix = hcat(Matrix(decision_matrix), scores)
    titles = human_readable_name(string.(vcat(criteria, :agregate_score)); title_case=true)

    return ADRIA.viz.map!(
        g,
        rs,
        outputs_matrix,
        titles;
        opts=opts,
        axis_opts=axis_opts
    )
end

"""
    ADRIA.viz.mcda_options(data::DataFrame; fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(), axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE())::Figure

Plot heatmap with weights given to the MCDA criteria for different options. This method is used in the pathway diversity
analysis where options are defined based on different set of weights given to the MCDA method.

# Example
```julia
options = ADRIA.analysis.option_seed_preference(include_weights=true)
fig = ADRIA.viz.mcda_options(options)
save("option_weight_matrix.png", fig)
```

# Arguments
- `data` : DataFrame with the weights given to different options.
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `fig_opts` : Additional options to pass to adjust Figure creation
  See: https://docs.makie.org/v0.19/api/index.html#Figure

# Returns
Figure
"""
function ADRIA.viz.mcda_options(
    data::DataFrame;
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    option_names = string.(data.option_name)
    criteria_names = string.(names(data))
    criteria_names = criteria_names[startswith.(criteria_names, "seed_")]
    number_options = size(data, 1)
    data = Matrix(data[:, criteria_names])

    fig = Figure(; fig_opts...)
    ax = Axis(
        fig[1, 1];
        xlabel="MCDA criteria",
        ylabel="Options",
        xticks=(1:length(criteria_names), criteria_names),
        yticks=(1:number_options, option_names),
        xticklabelrotation=2.0 / π,
        axis_opts...
    )
    limits = (0.0, 1.0)

    heatmap!(ax, transpose(data); colorrange=limits)
    Colorbar(fig[1, 2]; limits=limits)
    text!(ax,
        string.(round.(vec(data); digits=1));
        position=[
            Point2f(x, y) for x in 1:length(criteria_names) for y in 1:number_options
        ],
        align=(:center, :center),
        color=:white,
        fontsize=14
    )
    return fig
end

"""
    _default_colormap(rank_groups::Dict{Symbol,BitVector}, alpha_vals::Dict{Symbol,Float64})

Retrieve set of colormaps for plotting overlayed colormaps.

# Arguments
- `rank_groups` : Maps identifying key to be plotted to Boolean vector indicating members of group.
- `alpha_vals` : Maps identifying key to alpha values for colormap of each group (as greated by `alphas()`).

# Returns
Maps for each key in rank_groups to a unique colormap.
"""
function _default_colormap(
    rank_groups::Dict{Symbol,BitVector}, alpha_vals::Dict{Symbol,Float64}
)
    rank_colors = colors(rank_groups, alpha_vals)
    rank_ids = keys(rank_groups)
    return Dict(
        rank_grp =>
            [RGBA{Float32}(1.0, 1.0, 1.0, 0.1), rank_colors[rank_grp]] for
        rank_grp in rank_ids
    )
end
