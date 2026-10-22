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
        -19.68522515680423,
        -17.591136914659057,
        -16.32686454604817,
        -15.256465203867311,
        -16.66598105955471,
        -17.21205286917644,
        -16.084803868416564,
        -18.384910375942535,
        -17.374491746067957,
        -17.325375928477968,
        -17.359188938307504,
        -17.706302907144558
    ]
end

function growth_acceleration_height()::Vector{Float64}
    return [                    # ! Delte this
        1.5952988993156734,     #  1.7069464690106824
        1.7832100877157462,     #  0.13327757311686647
        1.1220310688452086,     #  0.06938495146223728
        1.4475376403793847,     #  1.9091802669096691
        1.7531081632567782,     #  1.1612734347453693
        0.6128194769859905,     #  0.4870627529189771
        1.693770466602163,      #  0.5622354107101416
        1.7126974382099205,     #  0.739032075821341
        1.9716717190691135,     #  0.03667967196234813
        1.8879796302427316,     #  0.038252983917947986
        1.6094528995868778,     #  0.2117497257815088
        0.3826377911360482      #  0.796979838518419
    ]
end

function growth_acceleration_midpoint()::Vector{Float64}
    return [
        0.24980484137954678,
        0.19230264057930244,
        0.1460734416151903,
        0.2108835165769789,
        0.28737704128742997,
        0.034577445361145834,
        0.11812341093936915,
        0.12478029024723014,
        0.2485926828557181,
        0.09162534286434126,
        0.2741192511849479,
        0.2385261343465052
    ]
end

# Probably need a ClusterGrowthAcceleration as well
function _growth_acceleration_struct(field_defs::OrderedDict)::Nothing
    s = IOBuffer()
    write(s, "Base.@kwdef struct GrowthAcceleration <: EcoModel\n")

    for (f, v) in field_defs
        write(s, "$(f)::Param = $(v)\n")
    end

    write(s, "end")
    eval(Meta.parse(String(take!(s))))

    return nothing
end

function _build_growth_acceleration_fields(
    bounds_var::Float64=0.1;
    overrides::Dict{String,Float64}=Dict{String,Float64}()
)::OrderedDict{String,Param}
    struct_fields = OrderedDict{String,Param}()

    growth_accel_params = [:steepness, :midpoint, :height]
    n_cb_calib_groups = length(growth_acceleration_midpoint())
    cb_calib_groups = 1:n_cb_calib_groups

    factor_val::Float64 = 0.0
    factor_name::String = ""
    for param in growth_accel_params
        for group in cb_calib_groups
            factor_name = "growth_acceleration_cb_group_$(group)_" * String(param)
            factor_val = get(
                overrides, factor_name, _growth_acceleration_values(param)[group]
            )

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

    return struct_fields
end

function create_growth_acceleration_struct(
    bounds_var::Float64=0.1;
    overrides::Dict{String,Float64}=Dict{String,Float64}()
)::Nothing
    _growth_acceleration_struct(
        _build_growth_acceleration_fields(bounds_var; overrides=overrides)
    )
    return nothing
end

function create_growth_acceleration_instance(
    bounds_var::Float64=0.1;
    overrides::Dict{String,Float64}=Dict{String,Float64}()
)::GrowthAcceleration
    fields = _build_growth_acceleration_fields(bounds_var; overrides=overrides)
    return GrowthAcceleration(; (Symbol(k) => v for (k, v) in fields)...)
end

create_growth_acceleration_struct()
