"""
    human_readable_name(names::Array{String}, title_case::Bool)

Make presentable parameter labels.
Returns a copy of original array so input is not modified.


Parameters
----------
names : parameter names to convert
title_case : boolean, whether to convert to Title Case or not.


Returns
-------
converted_names : array[str], of cleaned parameter names
"""
function human_readable_name(names::Array{String}, title_case::Bool=false)::Array{String}
    converted_names = names[:]
    converted_names = map((x) -> replace(x, "_"=>" "), converted_names)
    
    if title_case
        map!(titlecase, converted_names, converted_names)
    end

    return converted_names
end