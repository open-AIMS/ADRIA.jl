function outcome_title(outcomes::YAXArray)::String
    return get(outcomes.properties, :metric_name, "")
end

function set_plot_opts!(
    outcomes::AbstractArray,
    opts::OPT_TYPE,
    opts_key::Symbol;
    metadata_key::Symbol=:metric_feature,
    label_case=titlecase
)
    if !haskey(opts, opts_key) && (outcomes isa YAXArray)
        opts[opts_key] = outcome_label(
            outcomes; metadata_key=metadata_key, label_case=label_case
        )
    end
end

function outcome_label(
    outcomes::YAXArray; metadata_key::Symbol=:metric_feature, label_case=titlecase
)::String
    outcome_metadata = outcomes.properties

    return if all(haskey.([outcome_metadata], [metadata_key, :metric_unit]))
        _metric_feature = label_case(outcome_metadata[metadata_key])
        _metric_label = if !isempty(outcome_metadata[:metric_unit])
            "[$(outcome_metadata[:metric_unit])]"
        else
            ""
        end
        "$(_metric_feature) $(_metric_label)"
    else
        ""
    end
end
