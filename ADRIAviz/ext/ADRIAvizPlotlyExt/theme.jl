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
