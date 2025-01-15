using ADRIA: metrics
using ADRIA: DataFrame, DataFrameRow, names
using ADRIA: YAXArray, DataCube, ZeroDataCube

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

@testset "metrics.jl" begin
    n_scenarios::Int64 = size(TEST_SCENS, 1)
    n_groups::Int64 = length(ADRIA.default_coral_spec().taxa_names)
    test_scens_datacube::YAXArray{Float64,2} = DataCube(
        Matrix(TEST_SCENS); scenarios=1:n_scenarios, factors=names(TEST_SCENS)
    )

    _test_covers::Vector{YAXArray{Float64,4}} = mock_covers()
    _k_area::Vector{Float64} = k_area()

    @testset "relative_cover" begin
        test_metric(metrics.relative_cover, (TEST_RS,))
        for cover in _test_covers
            test_metric(metrics.relative_cover, (cover,))
        end
    end

    @testset "total_absolute_cover" begin
        test_metric(metrics.total_absolute_cover, (TEST_RS,))
        for cover in _test_covers
            test_metric(
                metrics.total_absolute_cover, (metrics.relative_cover(cover), _k_area)
            )
        end
    end

    @testset "relative_taxa_cover" begin
        test_metric(metrics.relative_taxa_cover, (TEST_RS,))

        for cover in _test_covers
            test_metric(
                metrics.relative_taxa_cover, (cover[:, :, :, 1], _k_area, n_groups)
            )
        end
    end

    @testset "relative_loc_taxa_cover" begin
        for cover in _test_covers
            test_metric(
                metrics.relative_loc_taxa_cover, (cover[:, :, :, 1], _k_area, n_groups)
            )
        end
    end

    @testset "relative_juveniles" begin
        scen_idx = 1
        coral_spec::DataFrame = ADRIA.to_coral_spec(ADRIA.default_coral_params(), TEST_SCENS[scen_idx, :])
        test_metric(metrics.relative_juveniles, (TEST_RS,))
        for cover in _test_covers
            test_metric(
                metrics.relative_juveniles, (cover[:, :, :, scen_idx], coral_spec)
            )
        end
    end

    @testset "absolute_juveniles" begin
        scen_idx = 1
        coral_spec::DataFrame = ADRIA.to_coral_spec(ADRIA.default_coral_params(), TEST_SCENS[scen_idx, :])
        test_metric(metrics.absolute_juveniles, (TEST_RS,))
        for cover in _test_covers
            test_metric(
                metrics.absolute_juveniles, (cover[:, :, :, scen_idx], coral_spec, _k_area)
            )
        end
    end

    @testset "juvenile_indicator" begin
        scen_idx = 1
        coral_spec::DataFrame = ADRIA.to_coral_spec(ADRIA.default_coral_params(), TEST_SCENS[scen_idx, :])
        test_metric(metrics.juvenile_indicator, (TEST_RS,))
        for cover in _test_covers
            test_metric(
                metrics.juvenile_indicator, (cover[:, :, :, scen_idx], coral_spec, _k_area)
            )
        end
    end

    @testset "coral_evenness" begin
        test_metric(metrics.coral_evenness, (TEST_RS,))
        for cover in _test_covers
            test_metric(
                metrics.coral_evenness,
                (
                    metrics.relative_loc_taxa_cover(cover[:, :, :, 1], _k_area, n_groups),
                )
            )
        end
    end

    @testset "absolute_shelter_volume" begin
        test_metric(metrics.absolute_shelter_volume, (TEST_RS,))
        for cover in _test_covers
            test_metric(
                metrics.absolute_shelter_volume, (cover, _k_area, TEST_SCENS)
            )
            test_metric(
                metrics.absolute_shelter_volume,
                (cover[:, :, :, 1], _k_area, TEST_SCENS[1, :])
            )
            test_metric(
                metrics.absolute_shelter_volume, (cover, _k_area, test_scens_datacube)
            )
            test_metric(
                metrics.absolute_shelter_volume,
                (
                    cover[:, :, :, 1], _k_area, test_scens_datacube[1, :]
                )
            )
        end
    end

    @testset "relative_shelter_volume" begin
        # For the two cases below, call three versions of relative_shelter_volumes
        scens_row::DataFrameRow = TEST_SCENS[2, :]
        scens_df::DataFrame = DataFrame(TEST_SCENS[2, :])
        scens_cube::YAXArray = DataCube(
            Matrix(scens_df); scenarios=axes(scens_df, 1), factors=names(scens_df)
        )

        @testset "one scenario case" begin
            for scens in [scens_row, scens_df, scens_cube]
                rsv = ADRIA.metrics.relative_shelter_volume(
                    _test_covers[1][:, :, :, 1], _k_area, scens
                )

                @test all(0.0 .<= rsv .<= 1.0)
                # Warn if all values are very tiny
                # catch Issue #91 : https://github.com/open-AIMS/ADRIA.jl/issues/91)
                @test any(rsv .>= 0.05)
            end
        end

        @testset "multi-scenario case" begin
            for scens in [scens_row, scens_df, scens_cube]
                rsv = ADRIA.metrics.relative_shelter_volume(
                    _test_covers[1], _k_area, TEST_SCENS
                )

                @test all(0.0 .<= rsv .<= 1.0) ||
                    "Min CC: $(minimum(sum(coral_cover, dims=:species)));
                    Max CC: $(maximum(sum(coral_cover, dims=:species))) |
                    $((minimum(rsv), maximum(rsv)))"
                @test any(rsv .>= 0.05)
            end
        end

        # Zero cover case
        zero_rsv = ADRIA.metrics.relative_shelter_volume(
            _test_covers[2], _k_area, DataFrame(TEST_SCENS)
        )
        @test all(zero_rsv .== 0.0)

        # Full cover case
        full_rsv = ADRIA.metrics.relative_shelter_volume(
            _test_covers[3], _k_area, DataFrame(TEST_SCENS)
        )
        @test all(0.0 .<= full_rsv .<= 1.0)

        # TODO Maximum shelter volume case.
        # TODO Check if species 24 is still the one with maximum shelter density
        # TODO Bring that to coral_cover Factories
        # max_coral_cover = DataCube(
        #     rand(n_timesteps, n_group_and_size, n_locations, n_scenarios),
        #     (:timesteps, :species, :locations, :scenarios)
        # )
        # Coral type with maximum shelter density
        # max_coral_cover[species=24, locations=1:3] .= _k_area'
        # max_rsv = ADRIA.metrics.relative_shelter_volume(
        #     max_coral_cover, _k_area, DataFrame(TEST_SCENS)
        # )
        # @test all(max_rsv .== 1.0) ||
        #     "Scenario with complete coral cover does not achieve max RSV |
        #     $(maximum(max_rsv))"
    end
end
