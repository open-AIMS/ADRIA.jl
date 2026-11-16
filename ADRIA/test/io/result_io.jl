
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

@testset "Combine result sets with different scenario columns" begin
    using ADRIA: DataFrames
    # A source without two of the factors (as from a run of another model version)
    dropped = [:cots_on, :cots_juvenile_predation]
    rs_a = let r = TEST_RS
        ADRIA.@set! r.inputs = DataFrames.select(r.inputs, DataFrames.Not(dropped))
        r
    end
    @test !(:cots_on in propertynames(rs_a.inputs))

    for (a, b) in ((TEST_RS, rs_a), (rs_a, TEST_RS))
        combined = ADRIA.combine_results(a, b)
        n_a, n_b = size(a.inputs, 1), size(b.inputs, 1)
        @test size(combined.inputs, 1) == n_a + n_b
        # Union of the columns, in order of first appearance, all Float64
        @test names(combined.inputs) ==
            names(DataFrames.vcat(a.inputs, b.inputs; cols=:union))
        @test all(eltype(c) == Float64 for c in eachcol(combined.inputs))
        # Missing factors are NaN for the source lacking them, values kept elsewhere
        for col in dropped
            lacking = a === rs_a ? (1:n_a) : ((n_a + 1):(n_a + n_b))
            has = a === rs_a ? ((n_a + 1):(n_a + n_b)) : (1:n_a)
            @test all(isnan, combined.inputs[lacking, col])
            @test combined.inputs[has, col] == TEST_RS.inputs[:, col]
        end
        @test ADRIA.scenario_sources(combined) == [fill(1, n_a); fill(2, n_b)]
    end

    # The union specification lists every factor once
    combined = ADRIA.combine_results(rs_a, TEST_RS)
    ms = ADRIA.model_spec(combined)
    @test allunique(zip(ms.component, ms.fieldname))
    @test issubset(TEST_RS.model_spec.fieldname, ms.fieldname)
end

@testset "Provenance is stored, loaded, merged and optional" begin
    prov = ADRIA.provenance(TEST_RS)
    @test length(prov) == 1
    @test prov[1]["adria_version"] == string(pkgversion(ADRIA))
    @test prov[1]["models"] == []  # COTS off, nothing attached
    @test prov[1]["overrides"] == []

    merged = ADRIA.combine_results(TEST_RS, TEST_RS)
    @test ADRIA.provenance(merged) == [prov[1], prov[1]]
    @test length(ADRIA.provenance(merged)) == length(merged.source_sizes)

    # A store written without provenance loads with `nothing`; merging it keeps alignment
    old_dir = mktempdir()
    cp(TEST_RS.env_layer_md.dpkg_path, old_dir; force=true)
    fn = joinpath(old_dir, ADRIA.INPUTS, ".zattrs")
    attrs = ADRIA.JSON.parsefile(fn)
    delete!(attrs, "provenance")
    write(fn, ADRIA.JSON.json(attrs))
    old_rs = ADRIA.load_results(old_dir)
    @test ADRIA.provenance(old_rs) == [nothing]
    mixed = ADRIA.combine_results(old_rs, TEST_RS)
    @test ADRIA.provenance(mixed) == [nothing, prov[1]]
end

@testset "Auxiliary logs" begin
    using ADRIA: Zarr

    @testset "no auxiliary logs without COTS" begin
        @test TEST_RS.aux == NamedTuple()
        @test isnothing(TEST_RS.cots_pop_log)  # legacy accessors give nothing
        @test isnothing(TEST_RS.cots_condition_log)
        @test_throws "no property or auxiliary log" TEST_RS.not_a_log
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

    @testset "aux LogSpec round-trips through a store, merged attrs equal the source's" begin
        spec = ADRIA.LogSpec(
            :synthetic_log, [:timesteps, :stage, :locations], [2, 2, 3];
            labels=Dict(:stage => ["a", "b"]), units="u", description="d"
        )
        ids = ["x", "y", "z"]
        dir = mktempdir()
        Zarr.zgroup(Zarr.DirectoryStore(dir), ADRIA.LOG_GRP)
        log_fn = joinpath(dir, ADRIA.LOG_GRP)
        ADRIA._setup_aux_logs(log_fn, [spec], ids, 4, 2; fill=Float32(NaN))
        grp = Zarr.zopen(log_fn)
        spec2 = ADRIA._aux_logspec_from_store(grp, :synthetic_log)
        @test all(getfield(spec2, f) == getfield(spec, f) for f in fieldnames(typeof(spec)))
        # A merged store builds its aux arrays from the rebuilt spec: same attributes
        dir2 = mktempdir()
        Zarr.zgroup(Zarr.DirectoryStore(dir2), ADRIA.LOG_GRP)
        log_fn2 = joinpath(dir2, ADRIA.LOG_GRP)
        ADRIA._setup_aux_logs(log_fn2, [spec2], ids, 8, 2; fill=Float32(NaN))
        a1 = Zarr.zopen(log_fn)["synthetic_log"].attrs
        a2 = Zarr.zopen(log_fn2)["synthetic_log"].attrs
        @test a1 == a2
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

        Zarr.zcreate(Float32, 2, 2; name="no_structure", path=joinpath(none, ADRIA.LOG_GRP))
        @test_throws "no_structure" ADRIA._load_aux_logs(
            Zarr.zopen(joinpath(none, ADRIA.LOG_GRP))
        )
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

        # A large aux log is in the sum, so it limits the chunk through `total`
        big = [ADRIA.LogSpec(:big, [:timesteps, :x], [tf, 10_000_000])]
        b = ADRIA._log_bytes_per_scen(TEST_DOM, big)
        @test b.largest == 4 * tf * 10_000_000
        @test ADRIA._max_chunk_by_size_for(TEST_DOM, big) ==
            max(1, floor(Int, 1_800_000_000 / b.total))
        @test ADRIA._max_chunk_by_size_for(TEST_DOM, big) <
            ADRIA._max_chunk_by_size_for(TEST_DOM, ADRIA.LogSpec[])
    end

    @testset "chunk cap is unchanged for COTS-off runs, also with optional logs" begin
        old_cap(tf, nl, ng, ni, nm) = max(
            1,
            floor(
                Int,
                1_800_000_000 /
                (4 * tf * (nl * nm + ng + nl * 2 + nl * ni + ng * nl * 2))
            )
        )
        tf, nl, ng, ns, ni, nm = 75, 3806, 5, 6, 3, 6
        total = 4 * tf * (nl * nm + ng + nl * 2 + nl * ni + ng * nl * 2)
        # largest array with both optional logs fully enabled (2-byte elements)
        largest = max(4 * tf * nl * nm, 2 * tf * ng * ns * nl)
        @test ADRIA._chunk_cap_by_bytes(total, largest) == old_cap(tf, nl, ng, ni, nm)

        # Same on the test domain with the optional logs switched on
        dom_cap(aux) = ADRIA._max_chunk_by_size_for(TEST_DOM, aux)
        off = dom_cap(ADRIA.LogSpec[])
        withenv("ADRIA_LOG_COVER" => "true", "ADRIA_LOG_DHW_TOLS" => "true") do
            @test dom_cap(ADRIA.LogSpec[]) == off
        end
        @test off == old_cap(
            length(ADRIA.timesteps(TEST_DOM)), TEST_DOM.coral_growth.n_locs,
            TEST_DOM.coral_growth.n_groups, length(ADRIA.interventions()),
            length(ADRIA.LOC_METRIC_NAMES)
        )
    end

    @testset "per-array 2 GB check" begin
        # Every counted array is part of the sum, so with real arrays `largest` can only
        # bind through the uncounted optional logs; test the rule on synthetic sizes.
        @test ADRIA._chunk_cap_by_bytes(1_000, 1_000) == 1_800_000
        @test ADRIA._chunk_cap_by_bytes(1_000, 90_000_000) == 20  # per-array bound binds
        @test ADRIA._chunk_cap_by_bytes(10_000_000_000, 1_000) == 1  # never below 1
    end
end
