"""
    to_string(m::Metric)::String

Get name of metric as a string.
"""
function to_string(m::Metric; is_titlecase=false)::String
    return to_string(m.func; is_titlecase=is_titlecase)
end
function to_string(f::Function; is_titlecase=false)::String
    metric_string::String = join(split(String(Symbol(f))[2:end], "_"), " ")
    return is_titlecase ? titlecase(metric_string) : metric_string
end

"""
    to_symbol(m::Metric)::String

Get name of metric as a symbol.
"""
function to_symbol(m::Metric)::Symbol
    return Symbol(replace(to_string(m), ' ' => '_'))
end

"""
    metric_label(m::Metric)::String
    metric_label(f::Function, unit::String)

Return name of metric in the format: "Title Case [Unit]", suitable for use as a label.

# Example
```julia
m_label = metric_label(scenario_total_cover)
# "Scenario Total Cover [m²]"
```
"""
function metric_label(m::Metric)::String
    return metric_label(m.func, m.unit)
end
function metric_label(f::Function, unit::String)::String
    n = titlecase(to_string(f))
    if length(unit) > 0
        n *= " [$unit]"
    end

    return n
end

"""
    dims(m::Metric)::Tuple

Get dimension names for a given outcome/metric.
"""
function dims(m::Metric)::Tuple
    return m.dims
end

"""
    ndims(m::Metric)::Int64

Infer the number of dimensions for a given outcome/metric.
"""
function Base.ndims(m::Metric)::Int64
    return length(dims(m))
end

"""
    call_metric(metric::Union{Function,Metric}, data::YAXArray, args...; kwargs...)

Convenience method that slices the data in the specified manner.

# Arguments
- `metric` : Function, the metric function to apply to "raw" data.
- `data` : YAXArray, data to pass into `metric`
- `args` : Additional positional arguments to pass into `metric`
- `kwargs` : Additional keyword arguments to pass into `slice_results`
    - `dims` : dummy keyword argument, not used but defined to allow use with other methods
"""
function call_metric(metric::Union{Function,Metric}, data::YAXArray, args...; kwargs...)
    dims = get(kwargs, :dims, nothing)
    if isnothing(dims)
        return metric(slice_results(data; kwargs...), args...)
    else
        return metric(slice_results(data; kwargs...), args...; dims=dims)
    end
end

"""
    slice_results(data::YAXArray; timesteps=(:), species=(:), locations=(:), scenarios=(:))

Slice data as indicated. Dimensions not found in target data are ignored.
"""
function slice_results(
    data::YAXArray; timesteps=(:), groups=(:), sizes=(:), locations=(:), scenarios=(:)
)::YAXArray
    f_dims = (
        timesteps=timesteps, groups=groups, sizes=sizes, locations=locations,
        scenarios=scenarios
    )

    s_names = keys(f_dims)
    d_names = axes_names(data)
    common_dims::Vector{Symbol} = intersect(s_names, d_names)

    selected_slice = (; zip(common_dims, [getfield(f_dims, k) for k in common_dims])...)
    return data[selected_slice...]
end
