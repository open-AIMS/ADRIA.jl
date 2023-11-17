using AxisKeys, NamedDims
using ADRIA: ResultSet

"""
    ADRIA.viz.ranks_to_frequencies!(g::Union{GridLayout,GridPosition},rs::ResultSet,
        frequencies::NamedDimsArray,rank_ids::Vector{Int64};opts::Dict=Dict(),axis_opts::Dict=Dict(),)
    ADRIA.viz.ranks_to_frequencies!(g::Union{GridLayout,GridPosition},rs::ResultSet,
        frequencies::NamedDimsArray,rank_id::Int64;opts::Dict=Dict(:color_map => :CMRmap),axis_opts::Dict=Dict())
    ADRIA.viz.ranks_to_frequencies(rs::ResultSet,frequencies::NamedDimsArray,rank_ids::Union{Int64,Vector{Int64}};
        opts::Dict=Dict(),fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
        
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
    frequencies::NamedDimsArray,
    rank_ids::Vector{Int64};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
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
    geodata = get_geojson_copy(rs)
    legend_els = Vector{Any}(undef, length(rank_ids))
    legend_labels = Vector{String}(undef, length(rank_ids))
    opts[:show_colorbar] = get(opts, :show_colorbar, false)

    ADRIA.viz.map!(
        g,
        rs,
        frequencies[ranks=rank_ids[1]];
        opts=opts,
        axis_opts=axis_opts,
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
            overdraw=true,
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
    frequencies::NamedDimsArray,
    rank_id::Int64;
    opts::Dict=Dict(),
    axis_opts::Dict=Dict())
    opts[:colorbar_label] = get(opts, :colorbar_label, "Selection frequency")
    opts[:color_map] = get(
        opts,
        :color_map,
        [RGBA{Float32}(1.0, 1.0, 1.0, 1.0), RGBA{Float32}(0.00784314, 0.243137, 1.0, 1.0)],
    )

    return ADRIA.viz.map!(
        g,
        rs,
        AxisKeys.keyless(NamedDims.unname(frequencies[ranks=rank_id]));
        opts=opts,
        axis_opts=axis_opts,
    )
end
function ADRIA.viz.ranks_to_frequencies(
    rs::ResultSet,
    frequencies::NamedDimsArray,
    rank_ids::Union{Int64,Vector{Int64}};
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(), axis_opts::Dict=Dict())

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.ranks_to_frequencies!(
        g,
        rs,
        frequencies,
        rank_ids;
        opts=opts,
        axis_opts=axis_opts,
    )

    return f
end

"""
    ADRIA.viz.decision_matrices(rs::ResultSet, S::NamedDimsArray, criteria::Vector{Symbol}; opts::Dict=Dict(),
        axis_opts::Dict=Dict(), fig_opts::Dict=Dict())
    ADRIA.viz.decision_matrices!(g::Union{GridLayout,GridPosition}, rs::ResultSet, S::NamedDimsArray, criteria::Vector{Symbol};
        opts::Dict=Dict(), axis_opts::Dict=Dict())
        
Plot a grid of spatial maps for a selection of location selection criteria.

# Arguments
- `g` : Figure GridPosition or GridLayout
- `rs` : Result set
- `S` : A decision matrix calculated using decison.decision_matrices
- `criteria` : Names of criteria to be plotted, e.g [:heat_stress, :wave_stress]
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
function ADRIA.viz.decision_matrices(
    rs::ResultSet,
    S::NamedDimsArray,
    criteria::Vector{Symbol};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.decision_matrices!(g, rs, S, criteria; opts=opts, axis_opts=axis_opts)
    return f
end
function ADRIA.viz.decision_matrices!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    S::NamedDimsArray,
    criteria::Vector{Symbol};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    n_criteria::Int64 = length(criteria)
    opts[:color_map] = get(opts, :color_map, :viridis)
    n_rows, n_cols = _calc_gridsize(n_criteria)
    criteria_names = String.(criteria)
    step::Int64 = 1

    for row in 1:n_rows, col in 1:n_cols
        axis_opts_temp = Dict(:title => criteria_names[step]; axis_opts...)
        ADRIA.viz.map!(
            g[row, col],
            rs,
            vec(S(criteria[step]));
            opts=opts,
            axis_opts=axis_opts_temp,
        )

        step += 1
        if step > n_criteria
            break
        end
    end

    try
        # Clear empty figures
        trim!(g)
    catch err
        if !(err isa MethodError)
            # GridPosition plots a single figure so does
            # not need empty figures to be cleared
            # If any other error is encountered, something else happened.
            rethrow(err)
        end
    end
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
