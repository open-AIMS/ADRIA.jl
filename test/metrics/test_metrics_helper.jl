using ADRIA: YAXArray
using ADRIA: metrics
using ADRIA: axes_names

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

function mock_covers()::Vector{YAXArray{Float64,4}}
    coral_spec = ADRIA.default_coral_spec()
    n_timesteps::Int64 = length(TEST_RS.env_layer_md.timeframe)
    n_groups::Int64 = length(coral_spec.taxa_names)
    n_sizes::Int64 = length(coral_spec.params.name) / n_groups
    n_group_sizes::Int64 = n_groups * n_sizes
    n_locations::Int64 = length(TEST_RS.coral_dhw_tol_log.locations)
    n_scenarios::Int64 = size(TEST_SCENS, 1)

    cover_params = (
        n_timesteps=n_timesteps,
        n_group_sizes=n_group_sizes,
        n_locations=n_locations,
        n_scenarios=n_scenarios
    )

    coral_cover::YAXArray{Float64,4} = MockData.coral_cover(; cover_params...)
    zero_coral_cover::YAXArray{Float64,4} = MockData.zero_coral_cover(; cover_params...)
    full_coral_cover::YAXArray{Float64,4} = MockData.full_coral_cover(; cover_params...)
    return [coral_cover, zero_coral_cover, full_coral_cover]
end

function k_area()::Vector{Float64}
    n_locations::Int64 = length(TEST_RS.coral_dhw_tol_log.locations)
    return round.(rand(n_locations) .* 100; digits=2)
end

function test_metric(metric::metrics.Metric, params::Tuple)::Nothing
    metric_result = metric(params...)
    @test all(0.0 .<= metric_result.data)
    metric.is_relative && @test all(metric_result.data .<= 1.0)

    _test_properties(metric, metric_result)
    return nothing
end

function _test_properties(metric::metrics.Metric, outcomes::YAXArray)
    prop_labels::NTuple{4,Symbol} = (:metric_name, :metric_unit, :axes_names, :axes_units)
    outcome_metadata = outcomes.properties
    all_properties_present = all(haskey.([outcome_metadata], prop_labels))
    @test all_properties_present
    if all_properties_present
        result_prop_axes_names::Vector{String} = outcome_metadata[:axes_names]
        result_prop_axes_units::Vector{String} = outcome_metadata[:axes_units]
        result_axes_names = string.(axes_names(outcomes))

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
