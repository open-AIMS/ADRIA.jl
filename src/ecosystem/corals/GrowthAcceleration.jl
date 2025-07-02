function _growth_acceleration_values(param_type::Symbol)
    if param_type == :steepness
        return growth_acceleration_steepness()
    elseif param_type == :height
        return growth_acceleration_height()
    elseif param_type == :midpoint
        return growth_acceleration_midpoint()
    end
    return error("Unknown `param_type` value")
end

function growth_acceleration_steepness()::Vector{Float64}
    return [
        -17.277489999368434,
        -16.76767223836814,
        -16.451332661396947,
        -15.656316048369732,
        -16.61001063114947,
        -17.56374341581085,
        -18.107227101710777,
        -18.271770553407237,
        -18.644593579949774,
        -17.87971779159888,
        -17.604091070767343,
        -18.930681306439357
    ]
end

function growth_acceleration_height()::Vector{Float64}
    return [
        1.9803441375451356,
        1.7852001282440557,
        1.7183361608960794,
        1.487143863701701,
        1.9813044333254952,
        1.72109896769677,
        1.4214857827821794,
        0.11786067208144713,
        1.8658714263220604,
        1.252973414720653,
        0.8206141621641347,
        1.1345347298346355
    ]
end

function growth_acceleration_midpoint()::Vector{Float64}
    return [
        0.18493420529513024,
        0.11703402591226619,
        0.230613793050364,
        0.0912121377221106,
        0.19989373135201927,
        0.17084115336762462,
        0.2668143900612246,
        0.21219714530869305,
        0.09524267706828961,
        0.14958572591149097,
        0.2984197990350005,
        0.10533468675338784
    ]
end

# Probably need a ClusterGrowthAcceleration as well
function _growth_acceleration_struct(field_defs::OrderedDict)::Nothing
    s = IOBuffer()
    write(s, "Base.@kwdef struct GrowthAcceleration{P} <: EcoModel\n")

    for (f, v) in field_defs
        write(s, "$(f)::P = $(v)\n")
    end

    write(s, "end")
    eval(Meta.parse(String(take!(s))))

    return nothing
end

function create_growth_acceleration_struct(
    bounds_var::Float64=0.1
)::Nothing
    struct_fields = OrderedDict{String,Param}()

    growth_accel_params = [:steepness, :midpoint, :height]
    n_cb_calib_groups = length(growth_acceleration_midpoint())
    cb_calib_groups = 1:n_cb_calib_groups

    factor_val::Float64 = 0.0
    factor_name::String = ""
    for param in growth_accel_params
        for group in cb_calib_groups
            factor_val = _growth_acceleration_values(param)[group]
            factor_name = "growth_acceleration_" * String(param) * "_$(group)"

            lower_bound = factor_val - bounds_var * abs(factor_val)
            upper_bound = factor_val + bounds_var * abs(factor_val)

            struct_fields[factor_name] = Factor(
                factor_val;
                ptype="continuous",
                dist=Uniform,
                dist_params=(lower_bound, upper_bound),
                name=human_readable_name(factor_name; title_case=true),
                description=""
            )
        end
    end

    _growth_acceleration_struct(struct_fields)

    return nothing
end

create_growth_acceleration_struct()
