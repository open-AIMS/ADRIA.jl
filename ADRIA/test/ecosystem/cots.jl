using Test
using ADRIA

@testset "COTS Submodel" begin
    @testset "CotsParams factors" begin
        defaults = ADRIA.CotsParams()
        @test defaults.cots_on.val == 0
        for mech in
            (:cots_adult_senescence, :cots_allee_stretched, :cots_juvenile_predation,
            :cots_per_capita_consumption, :cots_preferred_condition,
            :cots_preferred_prey_starvation)
            @test getfield(defaults, mech).val == 0
        end
        @test ADRIA.COTS_FACTOR_NAMES == Set(string.(fieldnames(ADRIA.CotsParams)))

        inst = ADRIA.create_cots_params_instance(;
            overrides=Dict("cots_on" => 1.0, "cots_allee_threshold" => 3.0)
        )
        @test inst.cots_on.val == 1
        @test inst.cots_allee_threshold.val == 3.0
        # Constant at the value by default; the full bounds are kept for sampling studies
        @test inst.cots_allee_threshold.dist_params == [3.0, 3.0, 3.0]
        @test inst.cots_allee_threshold.default_dist_params ==
            defaults.cots_allee_threshold.default_dist_params
        @test defaults.cots_on.dist_params == [0.0, 0.0]
        # Factors for the data package's COTS inputs are neutral by default
        @test defaults.cots_separated_recruitment.val == 0
        @test defaults.cots_larval_conn_mode.val == 0
        @test defaults.cots_external_outbreak_production.val == 0.0
        @test defaults.cots_initial_state_scale.val == 1.0
    end

    @testset "calibrated COTS values load by label" begin
        YA = ADRIA.YAXArrays
        ds = YA.Dataset(;
            cots=YA.YAXArray(
                (
                    ADRIA.DimensionalData.Dim{:cots_param}([
                        "cots_juvenile_predation_capacity", "cots_on"
                    ]),
                ),
                [25.0, 1.0]
            )
        )
        overrides = ADRIA._cots_calib_overrides(ds)
        @test overrides["cots_juvenile_predation_capacity"] == 25.0
        @test overrides["cots_on"] == 1.0
        @test isempty(
            ADRIA._cots_calib_overrides(
                YA.Dataset(;
                    other=YA.YAXArray(
                        (ADRIA.DimensionalData.Dim{:x}(1:2),), [1.0, 2.0])
                )
            )
        )
    end

    @testset "factors map onto COTSMod parameters" begin
        # Every mapped factor exists, and every target is a COTSMod parameter
        @test all(in(fieldnames(ADRIA.CotsParams)), keys(ADRIA.COTS_PARAM_MAP))
        @test all(in(fieldnames(ADRIA.COTSMod.COTSParams)), values(ADRIA.COTS_PARAM_MAP))
        @test allunique(values(ADRIA.COTS_PARAM_MAP))
        @test all(startswith("cots_"), string.(fieldnames(ADRIA.CotsParams)))
        names = collect(string.(fieldnames(ADRIA.CotsParams)))
        vals = [ADRIA.COTS_FACTOR_SPECS[Symbol(n)].val for n in names]
        ps = ADRIA.DataCube(Float64.(vals); factors=names)
        p0 = ADRIA.cots_runtime_params(ps)
        @test !p0.per_capita_consumption && !p0.juvenile_predation
        @test p0.senescence_age isa Int && p0.senescence_age == 6

        vals[findfirst(==("cots_per_capita_consumption"), names)] = 1.0
        vals[findfirst(==("cots_juvenile_predation_capacity"), names)] = 25.0
        p1 = ADRIA.cots_runtime_params(ADRIA.DataCube(Float64.(vals); factors=names))
        @test p1.per_capita_consumption
        @test p1.juvenile_predation_capacity == 25.0
        @test !p1.separated_recruitment

        vals[findfirst(==("cots_separated_recruitment"), names)] = 1.0
        vals[findfirst(==("cots_larval_fecundity"), names)] = 37.9
        p2 = ADRIA.cots_runtime_params(ADRIA.DataCube(Float64.(vals); factors=names))
        @test p2.separated_recruitment
        @test p2.larval_fecundity == 37.9
    end

    @testset "derived RNG seed ignores COTS factors" begin
        YA = ADRIA.YAXArrays
        base_names = ["guided", "eff_dhw_base", "RCP"]
        base_vals = [1.0, 0.63, 45.0]
        cots_names = sort(collect(ADRIA.COTS_FACTOR_NAMES))
        base = ADRIA.DataCube(base_vals; factors=base_names)
        with_cots = ADRIA.DataCube(vcat(base_vals, fill(7.5, length(cots_names)));
            factors=vcat(base_names, cots_names))
        withenv("ADRIA_RNG_SEED" => nothing) do
            @test rand(ADRIA.set_random_seed(base), 4) ==
                rand(ADRIA.set_random_seed(with_cots), 4)
        end
    end
end
