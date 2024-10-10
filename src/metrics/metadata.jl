const METRIC_METADATA = (:metric_name, :metric_feature, :is_relative, :metric_unit)
const AXIS_METADATA = (:axes_names, :axes_units)

"""
    metadata(outcomes::YAXArray)::Dict{Symbol,Any}

Helper function to extract metadata from YAXArrays.
"""
function metadata(outcomes::YAXArray)::Dict{Symbol,Any}
    return outcomes.properties
end

"""
    fill_metadata!(outcomes::YAXArray{T,N,A}, metric::Metric)::YAXArray{T,N,A} where {T,N,A}
    fill_metadata!(outcomes::YAXArray{T,N,A}, metadata::Dict{Symbol,Any})::YAXArray{T,N,A} where {T,N,A}

Fill outcomes YAXArray metadata (`properties` attribute).

# Arguments
- `outcomes` : YAXArray datacube of metric outcomes.
- `metric` : ADRIA.metrics.Metric object.
- `metadata` : Dict to be used to fill outcomes metrics metadata.
"""
function fill_metadata!(
    outcomes::YAXArray{T,N,A}, metric::Metric
)::YAXArray{T,N,A} where {T,N,A}
    fill_axes_metadata!(@views(outcomes))

    outcomes.properties[:metric_name] = to_string(metric; is_titlecase=true)
    outcomes.properties[:metric_feature] = metric.feature
    outcomes.properties[:is_relative] = metric.is_relative
    outcomes.properties[:metric_unit] = metric.unit

    return outcomes
end
function fill_metadata!(
    outcomes::YAXArray{T,N,A}, metadata::Dict{Symbol,Any}
)::YAXArray{T,N,A} where {T,N,A}
    # If `outcomes.properties`` is not a Dict{Symbol,Any} we need to build a new cube
    # because `properties`` type can't be updated
    if !(typeof(outcomes.properties) == Dict{Symbol,Any})
        _axis_names = axes_names(outcomes)
        _axis_labels = axis_labels.([outcomes], _axis_names)
        outcomes = DataCube(outcomes.data; NamedTuple{_axis_names}(_axis_labels)...)
    end

    fill_axes_metadata!(@views(outcomes))

    for metric_metadata in METRIC_METADATA
        if haskey(metadata, metric_metadata)
            outcomes.properties[metric_metadata] = metadata[metric_metadata]
        else
            @warn "Metric metadata \"$metric_metadata\" not found and won't be filled."
        end
    end

    return outcomes
end

"""
    fill_axes_metadata!(outcomes::YAXArray)::Nothing

Fill outcomes axes metadata.
"""
function fill_axes_metadata!(outcomes::YAXArray)::Nothing
    _axes_names::Tuple = axes_names(outcomes)
    outcomes.properties[:axes_names] = collect(
        parentmodule(metrics).human_readable_name.(
            _axes_names
        )
    )
    outcomes.properties[:axes_units] = collect(axes_units(_axes_names))
    return nothing
end

"""
    axes_units(axes_names::Union{Vector{Symbol},Tuple})::Tuple

Get units for each metric axis.
"""
function axes_units(axes_names::Union{Vector{Symbol},Tuple})::Tuple
    return values((timesteps="years", species="", locations="", scenarios="")[axes_names])
end
