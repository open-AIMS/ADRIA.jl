using YAXArrays

"""
    DataCube(data::AbstractArray; kwargs...)::YAXArray

Constructor for YAXArray. When used with `axes_names`, the axes labels will be UnitRanges
from 1 up to that axis length.

# Arguments
- `data` : Array of data to be used when building the YAXArray
- `axes_names` : Tuple of axes names
- `properties` : NamedTuple of metadata to be added to the YAXArray
"""
function DataCube(
    data::AbstractArray; properties::Dict{Symbol,Any}=Dict{Symbol,Any}(), kwargs...
)::YAXArray
    return YAXArray(Tuple(Dim{name}(val) for (name, val) in kwargs), data, properties)
end
function DataCube(
    data::AbstractArray, axes_names::Tuple; properties::Dict{Symbol,Any}=Dict{Symbol,Any}()
)::YAXArray
    return DataCube(
        data; properties=properties, NamedTuple{axes_names}(1:len for len in size(data))...
    )
end
function DataCube(
    data::AbstractArray, axes_names::Tuple, properties::Dict{Symbol,Any}
)::YAXArray
    return DataCube(
        data; properties=properties, NamedTuple{axes_names}(1:len for len in size(data))...
    )
end

"""
    ZeroDataCube(; T::DataType=Float64, kwargs...)::YAXArray
    ZeroDataCube(axes_names::Tuple, axes_sizes::Tuple; T::DataType=Float64)::YAXArray

Constructor for YAXArray with all entries equal zero. When `axes_name` and `axes_sizes`
are passed, all axes labels will be ranges.

# Arguments
- `axes_names` : Tuple of axes names
- `axes_sizes` : Tuple of axes sizes
- `properties` : NamedTuple of metadata to be added to the YAXArray
"""
function ZeroDataCube(;
    T::Type{D}=Float64, properties::Dict{Symbol,Any}=Dict{Symbol,Any}(), kwargs...
)::YAXArray where {D}
    return DataCube(
        zeros(T, [length(val) for (name, val) in kwargs]...); properties=properties,
        kwargs...
    )
end
function ZeroDataCube(
    axes_names::Tuple,
    axes_sizes::Tuple;
    properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    T::Type{D}=Float64
)::YAXArray where {D}
    return ZeroDataCube(;
        T=T, properties=properties, NamedTuple{axes_names}(1:size for size in axes_sizes)...
    )
end
function ZeroDataCube(
    axes_names::Tuple,
    axes_sizes::Tuple,
    properties::Dict{Symbol,Any};
    T::Type{D}=Float64
)::YAXArray where {D}
    return ZeroDataCube(;
        T=T, properties=properties, NamedTuple{axes_names}(1:size for size in axes_sizes)...
    )
end

"""
    axes_names(cube::YAXArray)::Tuple

Tuple of YAXArray axes names.
"""
function axes_names(cube::YAXArray)::Tuple
    return name.(cube.axes)
end

"""
    axis_labels(cube::YAXArray, axis_name::Symbol)::Vector{Any}

Vector of YAXArray axis labels.
"""
function axis_labels(cube::YAXArray, axis_name::Symbol)::Vector{Any}
    idx = axis_index(cube, axis_name)
    return cube.axes[idx].val.data
end

"""
    axis_index(cube::YAXArray, axis_name::Symbol)::Int64

YAXArray axis index.
"""
function axis_index(cube::YAXArray, axis_name::Symbol)::Int64
    if count(axes_names(cube) .== axis_name) > 1
        @warn "There are two or more axis with the same name. Returning the first."
    end
    return findfirst(axes_names(cube) .== axis_name)
end

"""
    sort_axis(cube::YAXArray, axis_name::Symbol)

Sorts axis labels of a YAXArray datacube.
"""
function sort_axis(cube::YAXArray, axis_name::Symbol)::YAXArray
    axis_labels = collect(lookup(cube, axis_name))
    labels_sort_idx::Vector{Int64} = sortperm(axis_labels)
    ordered_labels = axis_labels[labels_sort_idx]

    # Get selector to order data
    _axes_names = axes_names(cube)
    n_axes = length(_axes_names)
    selector::Vector{Union{Colon,Vector{Int64}}} = fill(:, n_axes)
    axis_idx::Int64 = findfirst(x -> x == axis_name, _axes_names)
    selector[axis_idx] = labels_sort_idx

    return cube[selector...]
end

"""
    copy_datacube(cube::YAXArray)::YAXArray

Copy a YAXArray data cube.

Alias for `copy(cube)`, kept to maintain backwards compatibility.
"""
function copy_datacube(cube::YAXArray)::YAXArray
    return copy(cube)
end

"""
    copy(cube::YAXArray)::YAXArray

Copy a YAXArray data cube.
"""
function Base.copy(cube::YAXArray)::YAXArray
    new_axlist = Tuple(ax for ax in deepcopy(cube.axes))
    return YAXArray(new_axlist, copy(cube.data))
end
