"""
    ADRIA.viz.map(rs::ADRIA.ResultSet; opts=Dict(by_RCP => false), fig_opts=Dict(), axis_opts=Dict(), series_opts=Dict())
    ADRIA.viz.map(rs::ADRIA.ResultSet, y::NamedDimsArray; opts=Dict(by_RCP => false), fig_opts=Dict(), axis_opts=Dict(), series_opts=Dict())
    ADRIA.viz.map!(f::Union{GridLayout,GridPosition}, rs::ADRIA.ResultSet, y::NamedDimsArray; opts=Dict(by_RCP => false), axis_opts=Dict(), series_opts=Dict())

Plot spatial outcomes.

# Arguments
- `rs` : ResultSet
- `y` : results of scenario metric
- `opts` : Aviz options 
    - `by_RCP`, color by RCP otherwise color by scenario type. Defaults to false.
- `axis_opts` : Additional options to pass to adjust Axis attributes  
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes  
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
GridPosition
"""
function ADRIA.viz.map!(g::Union{GridLayout,GridPosition}, rs::ResultSet, y::NamedDimsArray;
    opts::Dict=Dict(:by_RCP => false), axis_opts::Dict=Dict(), series_opts::Dict=Dict())

    # Ensure last year is always shown in x-axis
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(y)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)

    ax = Axis(
        g[1, 1],
        xticks=xtick_vals,
        xticklabelrotation=xtick_rot;
        axis_opts...
    )

    min_step = (1 / 0.05)
    color_weight = min((1.0 / (size(y, 1) / min_step)), 0.6)

    if :color ∉ keys(series_opts)
        hide_idx = get(series_opts, :hide_series, BitVector())

        if get(opts, :by_RCP, false) == false
            series_opts = merge(series_opts, Dict(:color => scenario_colors(rs, color_weight, hide_idx)))

            cf = LineElement(color=COLORS[:counterfactual], linestyle=nothing)
            ug = LineElement(color=COLORS[:unguided], linestyle=nothing)
            gu = LineElement(color=COLORS[:guided], linestyle=nothing)

            eles = [cf, ug, gu]
            labels = ["No Intervention", "Unguided", "Guided"]
        else
            rcp_ids = sort(Int.(unique(rs.inputs[:, :RCP])))
            c = [COLORS[_r] for _r in [Symbol("RCP$(r_id)") for r_id in rcp_ids]]
            r_s = Symbol[Symbol("RCP$(r_id)") for r_id in rcp_ids]

            eles = LineElement[
                LineElement(color=_c, linestyle=nothing)
                for _c in [COLORS[_r] for _r in r_s]
            ]

            series_opts = merge(series_opts, Dict(:color => map(x -> (COLORS[Symbol("RCP$(Int(x))")], color_weight), rs.inputs[:, :RCP])))
            labels = String.(r_s)
        end

        # Add legend
        Legend(g[1, 2], eles, labels, halign=:left, valign=:top, margin=(10, 10, 10, 10))
    end

    ls = series!(ax, y'; series_opts...)
    # ax.ylabel = metric_label(metric)
    ax.xlabel = "Year"

    return g
end
function ADRIA.viz.map(rs::ResultSet, y::NamedDimsArray; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict(), series_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.map!(g, rs, y; opts, axis_opts, series_opts)

    return f
end
function ADRIA.viz.map(rs::Union{Domain,ResultSet}; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict(), series_opts::Dict=Dict())
    f = Figure(; fig_opts...)

    geo_fn = make_geojson_copy(rs)
    geodata = GeoMakie.GeoJSON.read(read(geo_fn))
    data = Observable(rs.site_data.k)

    create_map!(f, geodata, data, nothing, ADRIA.centroids(rs))

    return f
end


"""
    make_geojson_copy(ds::Union{ResultSet,Domain})::String

# Arguments
ds : Domain or ResultSet containing spatial data

# Returns
Path to temporary copy of GeoJSON file.
"""
function make_geojson_copy(ds::Union{ResultSet,Domain})::String
    # Make temporary copy of GeoPackage as GeoJSON
    tmpdir = ADRIA.viz.tmpdir
    local geo_fn = joinpath(tmpdir, "Aviz_$(ds.name).geojson")
    if !isfile(geo_fn)
        try
            GDF.write(geo_fn, ds.site_data; driver="geojson")
        catch
            GDF.write(geo_fn, ds.site_data; geom_columns=(:geom,), driver="geojson")
        end
    end

    return geo_fn
end
