module ADRIAvizPlotlyKaleidoExt

using ADRIA, ADRIAviz, PlotlyBase, PlotlyKaleido

"""
    ADRIA.viz.savefig(p::PlotlyBase.Plot, path; width=900, height=600, scale=3)

Export a PlotlyBase figure to a static file (`.png`, `.svg`, `.pdf`).
`scale=3` maps 96-dpi screen resolution to ~300 dpi for print.
Requires `using PlotlyKaleido` to be loaded.
"""
function ADRIA.viz.savefig(p::PlotlyBase.Plot, path::AbstractString;
    width::Int=900, height::Int=600, scale::Int=3)
    PlotlyKaleido.savefig(p, path; width=width, height=height, scale=scale)
    return path
end

end  # module
