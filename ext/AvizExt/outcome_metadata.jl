function outcome_title(outcomes::YAXArray)::String
    return get(outcomes.properties, :metric_name, "")
end

function outcome_label(outcomes::YAXArray)::String
    outcome_metadata = outcomes.properties

    return if haskey(outcome_metadata, :metric_feature) &&
        haskey(outcome_metadata, :metric_unit)
        "$(outcome_metadata[:metric_feature]) [$(outcome_metadata[:metric_unit])]"
    else
        ""
    end
end
