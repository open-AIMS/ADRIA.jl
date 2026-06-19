"""
    _confint_traces(x, lo, mid, hi; name, color, alpha=0.2)

Returns `[band_trace, line_trace]` for a confidence-interval plot.
`color` must be a `"#RRGGBB"` hex string. Both traces have `type` set explicitly.
"""
function _confint_traces(x, lo, mid, hi; name::String, color::String, alpha::Real=0.2)
    x_vec = vec(x)
    fill_x = vcat(x_vec, reverse(x_vec))
    fill_y = vcat(vec(hi), reverse(vec(lo)))

    band = scatter(;
        x=fill_x,
        y=fill_y,
        fill="toself",
        fillcolor=_hex_to_rgba(color, alpha),
        line_color="rgba(0,0,0,0)",
        showlegend=false,
        hoverinfo="skip",
        name=name * "_band",
        type="scatter"
    )
    line = scatter(;
        x=x_vec, y=vec(mid), mode="lines",
        line_color=_hex_to_rgb(color), line_width=2,
        name=name, type="scatter"
    )

    return [band, line]
end

"""
    _group_colors(groups; palette=PLOTLY_COLORS)

Returns `Dict` mapping each group key to a Plotly RGB string.
Falls back to a generated color if the key is not in `palette`.
"""
function _group_colors(groups; palette=PLOTLY_COLORS)
    fallback = ["#377eb8", "#ff7f00", "#4daf4a", "#984ea3", "#e41a1c",
        "#a65628", "#f781bf", "#999999"]
    result = Dict{Symbol,String}()
    for (i, k) in enumerate(keys(groups))
        hex = get(ADRIAviz.COLORS, k, nothing)
        result[k] = isnothing(hex) ? fallback[mod1(i, length(fallback))] : hex
    end
    return result
end

"""
    _year_ticks(x_vals)

Returns `(tickvals, ticktext)` for a Plotly x-axis. Labels timestep integers
as calendar years starting from 2025.
"""
function _year_ticks(x_vals)
    tickvals = collect(x_vals)
    if first(tickvals) < 1900
        @warn "Timesteps do not appear to be calendar years. Defaulting to 2025 base year offset."
        base_year = 2025
        ticktext = string.(base_year .+ tickvals .- first(tickvals))
    else
        ticktext = string.(tickvals)
    end
    return tickvals, ticktext
end
