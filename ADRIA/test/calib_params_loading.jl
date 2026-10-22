using ADRIA: YAXArrays, DimensionalData
using YAXArrays: open_dataset, savedataset, Dataset
using DimensionalData: Dim

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
    const TEST_DATA_DIR = joinpath(ADRIA_DIR, "test", "data")
end

const MOCK_CALIB_PARAMS_PATH = joinpath(TEST_DATA_DIR, "mock_calibrated_params.nc")

# Pre-schema-v3 calibration output, written before `dist_std` was calibrated
const MOCK_CALIB_PARAMS_NO_DIST_STD_PATH = joinpath(
    TEST_DATA_DIR, "mock_calibrated_params_no_dist_std.nc"
)

@testset "Calibrated parameter loading" begin
    mock_nc_ds = open_dataset(MOCK_CALIB_PARAMS_PATH)

    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45"; calib_params_fn=MOCK_CALIB_PARAMS_PATH)
    ms = ADRIA.model_spec(dom)
    ms_dict = Dict(string(r.fieldname) => r.val for r in eachrow(ms))

    fg_names = string.(ADRIA.functional_group_names())
    bg_ids = collect(
        DimensionalData.lookup(mock_nc_ds["linear_extension_scale"], :cb_calib_group)
    )

    @testset "Coral per-species params" begin
        le = Array(mock_nc_ds["linear_extension"])
        mb = Array(mock_nc_ds["mb_rate"])
        dm = Array(mock_nc_ds["dist_mean"])
        ds = Array(mock_nc_ds["dist_std"])
        @test ms_dict["$(fg_names[1])_1_1_linear_extension"] ≈ le[1, 1]
        @test ms_dict["$(fg_names[3])_3_4_mb_rate"] ≈ mb[3, 4]
        @test ms_dict["$(fg_names[2])_2_5_dist_mean"] ≈ dm[2, 5]
        @test ms_dict["$(fg_names[2])_2_5_dist_std"] ≈ ds[2, 5]
        @test ms_dict["$(fg_names[5])_5_7_dist_std"] ≈ ds[5, 7]
    end

    # Legacy per-(group, size class) dist_std files must still load size-resolved: a loader
    # that collapsed the size dimension would still pass the single-element checks above
    @testset "2-D dist_std stays size-class resolved" begin
        ds = Array(mock_nc_ds["dist_std"])
        loaded = [ms_dict["$(fg_names[3])_3_$(sc)_dist_std"] for sc = 1:size(ds, 2)]
        @test loaded ≈ ds[3, :]
    end

    @testset "Scale factors" begin
        le_scale = Array(mock_nc_ds["linear_extension_scale"])
        mb_scale = Array(mock_nc_ds["mb_rate_scale"])
        @test ms_dict["linear_extension_scale_cb_group_$(bg_ids[3])_$(fg_names[1])"] ≈
            le_scale[1, 3]
        @test ms_dict["mb_rate_scale_cb_group_$(bg_ids[7])_$(fg_names[4])"] ≈ mb_scale[4, 7]
    end

    @testset "Growth acceleration" begin
        ga = Array(mock_nc_ds["growth_acceleration"])
        @test ms_dict["growth_acceleration_cb_group_$(bg_ids[1])_steepness"] ≈ ga[1, 1]
        @test ms_dict["growth_acceleration_cb_group_$(bg_ids[6])_height"] ≈ ga[6, 2]
        @test ms_dict["growth_acceleration_cb_group_$(bg_ids[12])_midpoint"] ≈ ga[12, 3]
    end
end

@testset "Calibration file without dist_std" begin
    dom = ADRIA.load_domain(
        TEST_DOMAIN_PATH, "45"; calib_params_fn=MOCK_CALIB_PARAMS_NO_DIST_STD_PATH
    )
    ms = ADRIA.model_spec(dom)
    ms_dict = Dict(string(r.fieldname) => r.val for r in eachrow(ms))

    legacy_ds = open_dataset(MOCK_CALIB_PARAMS_NO_DIST_STD_PATH)
    fg_names = string.(ADRIA.functional_group_names())

    # Calibrated parameters that are present still load
    @test ms_dict["$(fg_names[2])_2_5_dist_mean"] ≈ Array(legacy_ds["dist_mean"])[2, 5]

    # `dist_std` falls back to the ADRIA default rather than erroring
    default_dist_std = ADRIA.dist_std()
    @test ms_dict["$(fg_names[1])_1_1_dist_std"] ≈ default_dist_std[1]
end

@testset "Per-functional-group dist_mean / dist_std" begin
    fg_names = string.(ADRIA.functional_group_names())
    n_sizes = size(ADRIA.bin_widths(), 2)
    per_group = [1.5, 2.5, 3.5, 4.5, 5.5]
    per_group_mean = [4.0, 5.0, 6.0, 7.0, 8.0]

    ds = Dataset(;
        dist_std=ADRIA.YAXArrays.YAXArray(
            (Dim{:functional_group}(fg_names),), per_group
        ),
        dist_mean=ADRIA.YAXArrays.YAXArray(
            (Dim{:functional_group}(fg_names),), per_group_mean
        ),
        # `_coral_calib_overrides` unconditionally reads these two
        linear_extension_scale=ADRIA.YAXArrays.YAXArray(
            (Dim{:functional_group}(fg_names), Dim{:cb_calib_group}(1:1)),
            ones(length(fg_names), 1)
        ),
        mb_rate_scale=ADRIA.YAXArrays.YAXArray(
            (Dim{:functional_group}(fg_names), Dim{:cb_calib_group}(1:1)),
            zeros(length(fg_names), 1)
        )
    )

    overrides = ADRIA._coral_calib_overrides(ds)

    for (param, expected) in (("dist_std", per_group), ("dist_mean", per_group_mean))
        @testset "$param" begin
            @testset "One value broadcast across all size classes" begin
                for (fg_idx, fg) in enumerate(fg_names)
                    loaded = [
                        overrides["$(fg)_$(fg_idx)_$(sc)_$(param)"] for sc = 1:n_sizes
                    ]
                    @test all(loaded .≈ expected[fg_idx])
                end
            end

            @testset "Every size class is written" begin
                @test count(endswith("_$(param)"), keys(overrides)) ==
                    length(fg_names) * n_sizes
            end

            # Groups must stay distinct - a bug that broadcast a single scalar everywhere,
            # or transposed the repeat, would collapse them
            @testset "Groups remain distinct" begin
                first_of_each = [
                    overrides["$(fg)_$(fg_idx)_1_$(param)"]
                    for (fg_idx, fg) in enumerate(fg_names)
                ]
                @test first_of_each ≈ expected
            end
        end
    end

    # dist_mean and dist_std must not be crossed over
    @testset "mean and std stay separate" begin
        @test overrides["$(fg_names[1])_1_1_dist_mean"] ≈ per_group_mean[1]
        @test overrides["$(fg_names[1])_1_1_dist_std"] ≈ per_group[1]
    end
end

@testset "Depth attenuation" begin
    defaults = ADRIA.DepthAttenuation()

    @testset "Falls back to defaults when absent" begin
        # MOCK_CALIB_PARAMS_PATH predates the depth_attenuation variable
        dom = ADRIA.load_domain(
            TEST_DOMAIN_PATH, "45"; calib_params_fn=MOCK_CALIB_PARAMS_PATH
        )
        ms_dict = Dict(
            string(r.fieldname) => r.val for r in eachrow(ADRIA.model_spec(dom))
        )
        @test ms_dict["eff_dhw_base"] ≈ defaults.eff_dhw_base.val
        @test ms_dict["eff_dhw_mix"] ≈ defaults.eff_dhw_mix.val
    end

    @testset "Loads by label, not position" begin
        # Deliberately reversed relative to fieldnames(DepthAttenuation): a positional
        # loader would silently swap the two values.
        ds = Dataset(;
            depth_attenuation=ADRIA.YAXArrays.YAXArray(
                (Dim{:depth_atten_param}(["eff_dhw_mix", "eff_dhw_base"]),),
                [7.5, 0.123]
            )
        )
        overrides = ADRIA._depth_attenuation_calib_overrides(ds)
        @test overrides["eff_dhw_base"] ≈ 0.123
        @test overrides["eff_dhw_mix"] ≈ 7.5

        instance = ADRIA.create_depth_attenuation_instance(; overrides=overrides)
        @test instance.eff_dhw_base.val ≈ 0.123
        @test instance.eff_dhw_mix.val ≈ 7.5
        # Bounds come from the spec, never from the calibration file
        @test instance.eff_dhw_base.dist_params == defaults.eff_dhw_base.dist_params
    end

    @testset "mixing_scale = 0 is the well-mixed limit, not a NaN" begin
        # 0/0 in the general expression; must degrade to "no attenuation" instead
        @test ADRIA.effective_dhw_at_depth(0.0, 20.0, 0.04, 0.0) == 0.0
        @test ADRIA.effective_dhw_at_depth(9.0, 20.0, 0.04, 0.0) ≈ 9.0
        @test !isnan(ADRIA.effective_dhw_at_depth(0.0, 20.0, 0.04, 0.0))

        # κ_base = 0 likewise removes the depth refuge
        @test ADRIA.effective_dhw_at_depth(9.0, 20.0, 0.0, 12.0) ≈ 9.0
    end
end

@testset "ReefMod domains accept calibrated parameters" begin
    fg_names = string.(ADRIA.functional_group_names())
    dm = Array(open_dataset(MOCK_CALIB_PARAMS_PATH)["dist_mean"])
    ds = Array(open_dataset(MOCK_CALIB_PARAMS_PATH)["dist_std"])

    @testset "ReefModDomain" begin
        dom = ADRIA.load_domain(
            ReefModDomain,
            joinpath(TEST_DATA_DIR, "Reefmod_test_domain"),
            "45";
            calib_params_fn=MOCK_CALIB_PARAMS_PATH
        )
        ms_dict = Dict(
            string(r.fieldname) => r.val for r in eachrow(ADRIA.model_spec(dom))
        )
        @test ms_dict["$(fg_names[2])_2_5_dist_mean"] ≈ dm[2, 5]
        @test ms_dict["$(fg_names[2])_2_5_dist_std"] ≈ ds[2, 5]
    end

    @testset "RMEDomain" begin
        dom = ADRIA.load_domain(
            RMEDomain,
            joinpath(TEST_DATA_DIR, "RME_test_domain"),
            "45";
            calib_params_fn=MOCK_CALIB_PARAMS_PATH
        )
        ms_dict = Dict(
            string(r.fieldname) => r.val for r in eachrow(ADRIA.model_spec(dom))
        )
        @test ms_dict["$(fg_names[2])_2_5_dist_mean"] ≈ dm[2, 5]
        @test ms_dict["$(fg_names[2])_2_5_dist_std"] ≈ ds[2, 5]
    end
end
