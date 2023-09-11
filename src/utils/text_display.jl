"""
    human_readable_name(names::Array{String}, title_case::Bool)

Make presentable parameter labels.
Returns a copy of original array so input is not modified.


# Arguments
- names : parameter names to convert
- title_case : boolean, whether to convert to Title Case or not.


# Returns
- converted_names : array[str], of cleaned parameter names
"""
function human_readable_name(names::Vector{String}; title_case::Bool=false)::Vector{String}
    converted_names::Vector{String} = replace.(names[:], "_" => " ")

    if title_case
        map!(titlecase, converted_names, converted_names)
    end

    return converted_names
end
function human_readable_name(name::String; title_case::Bool=false)::String
    if title_case
        return titlecase(replace(name, "_" => " "))
    end

    return replace(name, "_" => " ")
end