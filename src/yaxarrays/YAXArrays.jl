function axes_names(cube::YAXArray)
    return name.(cube.axes)
end

"""
    sort(cube::YAXArray, axis_name::Symbol)

Sorts axis labels of a given YAXArray datacube.
"""
function sort(cube::YAXArray, axis_name::Symbol)::YAXArray
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
