module ADRIAvizPlotlyExt

using ADRIA
using ADRIAviz
using ADRIAviz: AnnotatedOutcomes, _get_scenario_groups, _scenario_clusters,
    _scenario_types, _scenario_rcps, COLORS, labels
using PlotlyBase

include("theme.jl")
include("helpers.jl")
include("viz/scenarios.jl")
include("viz/taxa_dynamics.jl")
include("viz/sensitivity.jl")
include("viz/clustering.jl")
include("viz/data_envelopment.jl")
include("viz/rule_extraction.jl")

# ── !-variant stubs ───────────────────────────────────────────────────────────
function ADRIA.viz.scenarios!(f, args...; kwargs...)
    return error(
        "The PlotlyBase backend does not support mutating (!) operations. " *
        "Use `scenarios(...)` to obtain a new Plot."
    )
end
function ADRIA.viz.clustered_scenarios!(f, args...; kwargs...)
    return error(
        "The PlotlyBase backend does not support mutating (!) operations. " *
        "Use `clustered_scenarios(...)` to obtain a new Plot."
    )
end
function ADRIA.viz.taxonomy!(f, args...; kwargs...)
    return error(
        "The PlotlyBase backend does not support mutating (!) operations. " *
        "Use `taxonomy(...)` to obtain a new Plot."
    )
end
function ADRIA.viz.pawn!(f, args...; kwargs...)
    return error(
        "The PlotlyBase backend does not support mutating (!) operations. " *
        "Use `pawn(...)` to obtain a new Plot."
    )
end
function ADRIA.viz.tsa!(f, args...; kwargs...)
    return error(
        "The PlotlyBase backend does not support mutating (!) operations. " *
        "Use `tsa(...)` to obtain a new Plot."
    )
end
function ADRIA.viz.rsa!(f, args...; kwargs...)
    return error(
        "The PlotlyBase backend does not support mutating (!) operations. " *
        "Use `rsa(...)` to obtain a new Plot."
    )
end
function ADRIA.viz.convergence!(f, args...; kwargs...)
    return error(
        "The PlotlyBase backend does not support mutating (!) operations. " *
        "Use `convergence(...)` to obtain a new Plot."
    )
end
function ADRIA.viz.outcome_map!(f, args...; kwargs...)
    return error(
        "The PlotlyBase backend does not support mutating (!) operations. " *
        "Use `outcome_map(...)` to obtain a new Plot."
    )
end
function ADRIA.viz.rules_scatter!(f, args...; kwargs...)
    return error(
        "The PlotlyBase backend does not support mutating (!) operations. " *
        "Use `rules_scatter(...)` to obtain a new Plot."
    )
end
function ADRIA.viz.data_envelopment_analysis!(f, args...; kwargs...)
    return error(
        "The PlotlyBase backend does not support mutating (!) operations. " *
        "Use `data_envelopment_analysis(...)` to obtain a new Plot."
    )
end

"""
    ADRIA.viz.show_in_browser(p::PlotlyBase.Plot)

Write `p` to a temporary HTML file and open it in the default system browser.
Useful when running Julia from a plain terminal (not Jupyter / Pluto / VS Code).
"""
function ADRIA.viz.show_in_browser(p::PlotlyBase.Plot)
    tmp = tempname() * ".html"
    open(tmp, "w") do io
        show(io, MIME"text/html"(), p)
    end
    if Sys.iswindows()
        run(`cmd /c start "" "$tmp"`; wait=false)
    elseif Sys.isapple()
        run(`open $tmp`; wait=false)
    else
        run(`xdg-open $tmp`; wait=false)
    end
    return tmp
end

end  # module
