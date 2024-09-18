using ADRIA: YAXArray
using ADRIA: metrics
using ADRIA: axes_names

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

function test_covers()::Vector{YAXArray{Float64,4}}
    n_timesteps::Int64 = length(TEST_RS.env_layer_md.timeframe)
    n_groups::Int64 = length(ADRIA.coral_spec().taxa_names)
    n_sizes::Int64 = length(ADRIA.coral_spec().params.name) / n_groups
    n_group_sizes::Int64 = n_groups * n_sizes
    n_locations::Int64 = length(TEST_RS.coral_dhw_tol_log.locations)
    n_scenarios::Int64 = size(TEST_SCENS, 1)

    cover_params = (
        n_timesteps=n_timesteps,
        n_group_sizes=n_group_sizes,
        n_locations=n_locations,
        n_scenarios=n_scenarios
    )

    coral_cover::YAXArray{Float64,4} = Factories.coral_cover(; cover_params...)
    zero_coral_cover::YAXArray{Float64,4} = Factories.zero_coral_cover(; cover_params...)
    full_coral_cover::YAXArray{Float64,4} = Factories.full_coral_cover(; cover_params...)
    return [coral_cover, zero_coral_cover, full_coral_cover]
end

function k_area()::Vector{Float64}
    n_locations::Int64 = length(TEST_RS.coral_dhw_tol_log.locations)
    return round.(rand(n_locations) .* 100; digits=2)
end

function _test_metrics(
    metric::metrics.Metric, params::Tuple; relative::Bool=false
)::Nothing
    metric_result = metric(params...)
    @test all(0.0 .<= metric_result.data)
    relative && @test all(metric_result.data .<= 1.0)

    _test_properties(metric, metric_result)
    return nothing
end

function _test_properties(metric::metrics.Metric, metric_result::YAXArray)
    prop_labels::NTuple{4,Symbol} = (:metric_name, :metric_unit, :axes_names, :axes_units)
    all_properties_present = all(haskey.([metric_result.properties], prop_labels))
    @test all_properties_present
    if all_properties_present
        result_prop_axes_names::Vector{String} = metric_result.properties[:axes_names]
        result_prop_axes_units::Vector{String} = metric_result.properties[:axes_units]
        result_axes_names = string.(axes_names(metric_result))

        @test metric_result.properties[:metric_name] == metrics.metric_label(metric)

        # All keys in result properties match YAXArray
        @test all(result_prop_axes_names .∈ [result_axes_names]) ||
            "Result prop axes_names: $result_prop_axes_names ;" *
              "Result axes names: $result_axes_names"

        # No key in result properties is  YAXArray
        @test !any(result_prop_axes_names .∉ [result_axes_names]) ||
            "Result prop axes_names: $result_prop_axes_names ;" *
              "Result axes names: $result_axes_names"

        # Number of axes_names match number of axes_units
        @test length(result_prop_axes_names) == length(result_prop_axes_units)
    end
end

function test_relative_metric(metric::metrics.Metric, params::Tuple)::Nothing
    return _test_metrics(metric, params; relative=true)
end

function test_absolute_metric(metric::metrics.Metric, params::Tuple)::Nothing
    return _test_metrics(metric, params; relative=false)
end
