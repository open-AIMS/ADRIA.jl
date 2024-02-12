import Base.copy

"""
    copy(cube::YAXArray)::YAXArray

Use julia's Base.copy function to create a copy of an YAXArray
"""
function copy(cube::YAXArray)::YAXArray
    new_data = copy(cube.data)
    new_axlist = Tuple(ax for ax in deepcopy(cube.axes))
    return YAXArray(new_axlist, new_data)
end
