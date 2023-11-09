"""
    human_readable_name(names::Array{String}, title_case::Bool)

Make presentable parameter labels.
Returns a copy of original array so input is not modified.


# Arguments
- `names` : parameter names to convert
- `title_case` : boolean, whether to convert to Title Case or not.


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

"""
    to_scientific(x::Float64; kwargs...)::String
    to_scientific(x::Int64; kwargs...)::String

Converts any Real number to scientific notation.

# Arguments
- `x` : Number to be converted
- `digits` : Optional argument. If present, represent the number of decimal digits. Defaults
to 2.
"""
function to_scientific(x::Float64; digits=2)::String
    x == 0.0 && return 0.0
    a, b = get_scientific_factors(x::Float64; digits=digits)
    return "$a.E$b"
end
function to_scientific(x::Int64; digits=2)::String
    return to_scientific(Float64(x); digits=digits)
end

"""
    get_scientific_factors(x::Float64; digits=2)::Tuple
    get_scientific_factors(x::Int64; digits=2)::Tuple

Return a and b to write any Real number x as x = a * 10^b, where 1.0 <= a < 10.0 and b is
an integer. The idea is that `x = a * 10^b` => `log(x) = log(a) + b*log(10)`. Which means
that `log(a)` is the remainder of `log(x)/log(10)` hence `a = exp((log(x)/log(10)))`.
Because `log(x) % log(10)` returns an approximate value, if we just use this to find `a` we
will end up writing `40` as `3.9E1`. To avoid that, we take this approximated value of `a`,
compute `b` and use `b` to find the actual value of `a`.

# Arguments
- `x` : Number to be converted
- `digits` : Optional argument. If present, represent the number of decimal digits. Defaults
to 2.
"""
function get_scientific_factors(x::Float64; digits=2)::Tuple
    x == 0.0 && return 0.0, 0

    if x < 1
        a::Float64, b::Int64 = _scientific_factors(1 / x)
        a = (1 / a) * 10
        b = -(b + 1)
    else
        a::Float64, b::Int64 = _scientific_factors(x)
    end

    a = round(a; digits=digits)

    return a, b
end
function get_scientific_factors(x::Int64; digits=2)::Tuple
    return to_scientific(Float64(x); digits=digits)
end

function _scientific_factors(x::Float64)
    # Handle negative numbers
    signal = x < 0.0 ? -1.0 : 1.0
    x = abs(x)

    # Handle edgecases x = 1 and x = -1
    x == 1.0 && return signal, 0

    # Use an approximate value of `a` to find `b` from log(a) ≈ log(x) % log(10)
    a_tmp::Float64 = exp(log(x) % log(10))
    b::Int64 = round(log(x / a_tmp) / log(10))

    # Use `b` to find the actual value of `a`
    a::Float64 = (x / (10^b)) * signal

    return a, b
end
