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
        -18.774799898803458,
        -15.643415963266387,
        -19.008589392588153,
        -18.564180311352413,
        -15.134761869820837,
        -17.47808759125323,
        -18.030748393417703,
        -17.853513329490728,
        -18.538008312231337,
        -15.603118544548202,
        -16.78149971688862,
        -17.485955617875728
    ]
end

function growth_acceleration_height()::Vector{Float64}
    return [
        1.8381981005996322,
        0.8467771352746775,
        0.9764968248915961,
        1.5431188865441479,
        0.5031583441336789,
        1.4008895565174615,
        0.11478237792331683,
        1.5263384987806172,
        1.6189122845714996,
        1.1345938896278691,
        0.12096555701867479,
        0.652335904450534
    ]
end

function growth_acceleration_midpoint()::Vector{Float64}
    return [
        0.2853001842125711
        0.081995019135015
        0.28353778696967435
        0.13625131504441762
        0.14524935506495604
        0.028992241541660448
        0.21723100692054645
        0.19765631776872394
        0.15995152211549393
        0.03738192795775167
        0.19700447286934764
        0.2009025085911272
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
            factor_name = "growth_acceleration_cb_group_$(group)_" * String(param)

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
