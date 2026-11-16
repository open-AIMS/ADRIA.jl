
if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

@testset "Log element types are Float64 on load" begin
    @test eltype(TEST_RS.ranks) == Float64
    @test eltype(TEST_RS.iv_LvM_log) == Float64
    @test eltype(TEST_RS.iv_CAq_log) == Float64
    @test eltype(TEST_RS.iv_Shd_log) == Float64
    @test eltype(TEST_RS.coral_dhw_tol_log) == Float64

    # coral_cover_log: enabled via log_cover=true in test/config.toml
    @test !isnothing(TEST_RS.coral_cover_log)
    @test eltype(TEST_RS.coral_cover_log) == Float64

    # Verify UInt16 unscaling: cover values must be in [0, 1]
    cover_data = TEST_RS.coral_cover_log.data
    @test all(0.0 .<= cover_data .<= 1.0)
end

@testset "Combine result sets" begin
    combined_rs = ADRIA.combine_results(TEST_RS, TEST_RS)
    @test combined_rs isa ADRIA.ResultSet
end

@testset "Auxiliary logs" begin
    using ADRIA: Zarr

    @testset "no auxiliary logs without COTS" begin
        @test TEST_RS.aux == NamedTuple()
        @test isnothing(TEST_RS.cots_pop_log)  # legacy accessors give nothing
        @test isnothing(TEST_RS.cots_condition_log)
        @test_throws ErrorException TEST_RS.not_a_log
    end

    @testset "COTS logs are declared as LogSpecs with unchanged store attributes" begin
        specs = ADRIA.cots_log_specs(TEST_DOM)
        @test [s.name for s in specs] == [:cots_pop_log, :cots_condition_log]
        tf = length(ADRIA.timesteps(TEST_DOM))
        n_locs = TEST_DOM.coral_growth.n_locs
        @test specs[1].size == [tf, 3, n_locs]
        @test specs[2].size == [tf, n_locs]

        ids = ["a", "b"]
        pop = ADRIA._aux_log_attrs(specs[1], ids)
        @test pop == Dict{Symbol,Any}(
            :structure => ("timesteps", "age_class", "locations", "scenarios"),
            :age_classes => ["recruits", "juveniles", "adults"],
            :units => "COTS ha^-1",
            :unique_loc_ids => ids
        )
        @test ADRIA._aux_log_attrs(specs[2], ids) == Dict{Symbol,Any}(
            :structure => ("timesteps", "locations", "scenarios"),
            :unique_loc_ids => ids
        )
    end

    @testset "load_results discovers aux logs, old fog/shade arrays stay core" begin
        dir = mktempdir()
        root = Zarr.zgroup(Zarr.DirectoryStore(dir))
        grp = Zarr.zgroup(root, ADRIA.LOG_GRP)
        for n in ("fog", "shade", "rankings", "coral_dhw_log")
            Zarr.zcreate(Float32, 2, 2; name=n, path=joinpath(dir, ADRIA.LOG_GRP))
        end
        Zarr.zcreate(
            Float32, 2, 3, 4; name="zz_custom_log", path=joinpath(dir, ADRIA.LOG_GRP),
            attrs=Dict("structure" => ["timesteps", "foo", "scenarios"])
        )
        Zarr.zcreate(
            Float32, 2, 3; name="aa_other_log", path=joinpath(dir, ADRIA.LOG_GRP),
            attrs=Dict("structure" => ["timesteps", "scenarios"])
        )
        aux = ADRIA._load_aux_logs(Zarr.zopen(joinpath(dir, ADRIA.LOG_GRP)))
        @test keys(aux) == (:aa_other_log, :zz_custom_log)
        @test size(aux.zz_custom_log) == (2, 3, 4)
        @test size(aux.zz_custom_log, :foo) == 3

        none = mktempdir()
        r = Zarr.zgroup(Zarr.DirectoryStore(none))
        Zarr.zgroup(r, ADRIA.LOG_GRP)
        Zarr.zcreate(Float32, 2, 2; name="fog", path=joinpath(none, ADRIA.LOG_GRP))
        @test ADRIA._load_aux_logs(Zarr.zopen(joinpath(none, ADRIA.LOG_GRP))) ==
            NamedTuple()
    end

    @testset "per-scenario byte budget keeps the sum and adds aux logs" begin
        tf = length(ADRIA.timesteps(TEST_DOM))
        nl = TEST_DOM.coral_growth.n_locs
        ng = TEST_DOM.coral_growth.n_groups
        n_iv = length(ADRIA.interventions())
        n_met = length(ADRIA.LOC_METRIC_NAMES)

        base = ADRIA._log_bytes_per_scen(TEST_DOM, ADRIA.LogSpec[])
        @test base.total ==
            4 * tf * (nl * n_met + ng + nl * 2 + nl * n_iv + ng * nl * 2)

        cots = ADRIA.cots_log_specs(TEST_DOM)
        with = ADRIA._log_bytes_per_scen(TEST_DOM, cots)
        @test with.total == base.total + 4 * tf * (3 * nl + nl)
        @test base.largest >= 4 * tf * ng * nl  # includes optional logs when enabled
        @test with.largest >= base.largest

        # A large aux log dominates the largest-array check
        big = [ADRIA.LogSpec(:big, [:timesteps, :x], [tf, 10_000_000])]
        b = ADRIA._log_bytes_per_scen(TEST_DOM, big)
        @test b.largest == 4 * tf * 10_000_000
        @test ADRIA._max_chunk_by_size_for(TEST_DOM, big) ==
            max(1, floor(Int, 1_800_000_000 / b.total))
        @test ADRIA._max_chunk_by_size_for(TEST_DOM, big) <
            ADRIA._max_chunk_by_size_for(TEST_DOM, ADRIA.LogSpec[])
    end
end
