using ADRIA: YAXArrays, DimensionalData
using YAXArrays: open_dataset, savedataset, Dataset
using DimensionalData: Dim

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
    const TEST_DATA_DIR = joinpath(ADRIA_DIR, "test", "data")
end

const MOCK_CALIB_PARAMS_PATH = joinpath(TEST_DATA_DIR, "mock_calibrated_params.nc")

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
        @test ms_dict["$(fg_names[1])_1_1_linear_extension"] ≈ le[1, 1]
        @test ms_dict["$(fg_names[3])_3_4_mb_rate"] ≈ mb[3, 4]
        @test ms_dict["$(fg_names[2])_2_5_dist_mean"] ≈ dm[2, 5]
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
