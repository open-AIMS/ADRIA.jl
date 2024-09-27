function outcome_title(outcomes::YAXArray)::String
    return get(outcomes.properties, :metric_name, "")
end

function outcome_label(outcomes::YAXArray; label_case=titlecase)::String
    outcome_metadata = outcomes.properties

    return if all(haskey.([outcome_metadata], [:metric_feature, :metric_unit]))
        _metric_feature = label_case(outcome_metadata[:metric_feature])
        "$(_metric_feature) [$(outcome_metadata[:metric_unit])]"
    else
        ""
    end
end
