function _hex_to_rgb(hex::String)
    h = lstrip(hex, '#')
    r, g, b = parse.(Int, (h[1:2], h[3:4], h[5:6]); base=16)
    return "rgb($r,$g,$b)"
end

function _hex_to_rgba(hex::String, alpha::Real)
    h = lstrip(hex, '#')
    r, g, b = parse.(Int, (h[1:2], h[3:4], h[5:6]); base=16)
    return "rgba($r,$g,$b,$alpha)"
end

const PLOTLY_COLORS = Dict{Symbol,String}(
    k => _hex_to_rgb(v) for (k, v) in ADRIAviz.COLORS
)

const ADRIA_LAYOUT_DEFAULTS = Layout(;
    font=attr(; family="Open Sans, sans-serif", size=12),
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=attr(; showgrid=false, linecolor="black", linewidth=1),
    yaxis=attr(; gridcolor="#e5e5e5", linecolor="black", linewidth=1)
)

# Tier values mirror set_typography_defaults! in src/viz/viz.jl.
# DPI note: Plotly renders at ~96 DPI (browser/screen); Makie targets 300 DPI print.
# The same numeric pt values are intentionally shared - the perceptual difference is
# acceptable because the two backends serve different audiences (dashboard vs. publication).
# Users exporting Plotly to PNG for print should override via the layout keyword.
function _plotly_font_sizes(n_panels::Int=1)
    if n_panels <= 1
        return (title=16, label=12, tick=10)
    elseif n_panels <= 4
        return (title=12, label=10, tick=8)
    else
        return (title=10, label=9, tick=8)
    end
end
