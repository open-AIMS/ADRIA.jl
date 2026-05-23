using YAXArrays
using DataFrames
using ADRIA.analysis: series_confint
using OrderedCollections

# ──────────────────────────────────────────────────────────────────────────────
# scenarios(ao::AnnotatedOutcomes) — primary Plotly dispatch
# ──────────────────────────────────────────────────────────────────────────────

"""
    ADRIA.viz.scenarios(ao::AnnotatedOutcomes; summarize=true, by_RCP=false, kwargs...)

Plot scenario outcomes over time from an `AnnotatedOutcomes` object.
"""
function ADRIA.viz.scenarios(
    ao::AnnotatedOutcomes;
    summarize::Bool=true,
    by_RCP::Bool=false,
    sort_by::Union{Symbol,Nothing}=nothing,
    alpha::Real=0.15,
    xlabel::String="Year",
    ylabel::String="",
    title::String="",
    kwargs...
)::PlotlyBase.Plot
    groups = _get_scenario_groups(ao; by_RCP=by_RCP)
    data = ao.data
    x_vals = collect(data.timesteps)
    tickvals, ticktext = _year_ticks(x_vals)
    colors = _group_colors(groups)
    mat = collect(data)  # (timesteps × scenarios) — avoids repeated YAXArray cat calls

    traces = PlotlyBase.AbstractTrace[]
    for (gname, mask) in groups
        any(mask) || continue
        subset = mat[:, mask]

        # subset is (timesteps × scenarios); series_confint expects (timesteps × scenarios)
        if ndims(subset) == 1
            subset = reshape(subset, :, 1)
        end

        color = get(colors, gname, "#4682b4")
        label = string(gname)

        if summarize
            ci = series_confint(subset)
            lo, mid, hi = ci[:, 1], ci[:, 2], ci[:, 3]
            ts = _confint_traces(x_vals, lo, mid, hi; name=label, color=color, alpha=alpha)
            append!(traces, ts)
        else
            for j in axes(subset, 2)
                push!(
                    traces,
                    PlotlyBase.scatter(;
                        x=x_vals, y=subset[:, j], mode="lines",
                        line_color=_hex_to_rgba(color, 0.5),
                        name=label, showlegend=j == 1, type="scatter"
                    )
                )
            end
        end
    end

    _ylabel =
        isempty(ylabel) ? ADRIAviz.outcome_label(data; metadata_key=:metric_name) : ylabel
    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(;
            title_text=xlabel,
            tickvals=tickvals,
            ticktext=ticktext
        ),
        yaxis=PlotlyBase.attr(; title_text=_ylabel)
    )
    return PlotlyBase.Plot(traces, layout)
end

# ──────────────────────────────────────────────────────────────────────────────
# scenarios(rs, outcomes) and scenarios(scens, outcomes) — convenience dispatches
# ──────────────────────────────────────────────────────────────────────────────

function ADRIA.viz.scenarios(
    rs::ADRIA.ResultSet,
    outcomes::YAXArray;
    kwargs...
)::PlotlyBase.Plot
    return ADRIA.viz.scenarios(rs.inputs, outcomes; kwargs...)
end

function ADRIA.viz.scenarios(
    scenarios::DataFrame,
    outcomes::YAXArray;
    by_RCP::Bool=false,
    kwargs...
)::PlotlyBase.Plot
    _scens = scenarios[collect(outcomes.scenarios), :]
    groups = if by_RCP
        _scenario_rcps(_scens)
    else
        _scenario_types(_scens)
    end
    ao = ADRIA.AnnotatedOutcomes(
        outcomes,
        Dict{Symbol,Any}(
            :scenario_type_groups => groups,
            :scenario_rcp_groups => by_RCP ? groups : nothing
        )
    )
    return ADRIA.viz.scenarios(ao; by_RCP=by_RCP, kwargs...)
end
