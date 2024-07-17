function unit_power(unit::Symbol)::Int64
    return (m=0, cm=-2)[unit]
end

function linear_scale(from_unit::Symbol, to_unit::Symbol)::Float64
    resulting_power = unit_power(from_unit) - unit_power(to_unit)
    return round(10.0^resulting_power, digits=abs(resulting_power))
end
function linear_scale(number::Real, from_unit::Symbol, to_unit::Symbol)::Float64
    return number * linear_scale(from_unit, to_unit)
end

function quadratic_scale(from_unit::Symbol, to_unit::Symbol)::Float64
    resulting_power = unit_power(from_unit) * 2 - unit_power(to_unit) * 2
    return round(10.0^resulting_power, digits=abs(resulting_power))
end
function quadratic_scale(number::Real, from_unit::Symbol, to_unit::Symbol)::Float64
    return number * quadratic_scale(from_unit, to_unit)
end
