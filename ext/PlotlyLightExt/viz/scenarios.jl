using Statistics

"""
    ADRIA.viz.scenarios(rs::ResultSet, outcomes::YAXArray)

Indicative scenario outcomes grouped by guidance type.
"""
function ADRIA.viz.scenarios(rs::ResultSet, outcomes::YAXArray)
    # Metadata for scenario types
    guided_id = sort(unique(rs.inputs.guided))
    guided_idx = Int.(guided_id .+ 2)
    guidance_type = _get_guided_labels()[guided_idx]
    guided_colors = _guided_colors()

    # Build plot for each scenario type
    res = plot()

    ts = outcomes.timesteps
    for (i, (g_id, g_str)) in enumerate(zip(guided_id, guidance_type))
        # Select subset of scenarios for each scenario type
        scen_sel = rs.inputs.guided .== g_id
        subset = outcomes[:, scen_sel]

        q50 = median(subset; dims=:scenarios)[:]
        ci_bounds = Matrix(
            mapreduce(
                collect,
                hcat,
                quantile.(eachrow(subset), [_calc_confint(0.95)]))'
        )
        ci_bound_opts = Dict(:width => 0, :color => guided_colors[i])

        res.scatter(;
            x=ts,
            y=ci_bounds[:, 1],
            mode="lines",
            name=g_str,
            line=ci_bound_opts,
            showlegend=false,
            opacity=0.05
        ).scatter(;
            x=ts,
            y=ci_bounds[:, 2],
            mode="lines",
            name=g_str,
            line=ci_bound_opts,
            showlegend=false,
            opacity=0.05,
            fill="tonexty"
        ).scatter(;
            x=ts,
            y=q50,
            name=g_str,
            line=Dict(:color => guided_colors[i])
        )
    end

    res.layout.yaxis.title.text = outcome_label(outcomes)
    res.layout.xaxis.title.text = "Years"

    return res
end
