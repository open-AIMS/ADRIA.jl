module ADRIAvizPlotlyExt

using ADRIA
using ADRIAviz
using ADRIAviz: AnnotatedOutcomes, _get_scenario_groups, COLORS, labels
using PlotlyLight

function ADRIA.viz.scenarios(ao::AnnotatedOutcomes; kwargs...)
    return error("PlotlyLight backend: scenarios not yet implemented")
end

function ADRIA.viz.taxonomy(ao::AnnotatedOutcomes; kwargs...)
    return error("PlotlyLight backend: taxonomy not yet implemented")
end

end  # module
