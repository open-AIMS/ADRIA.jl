
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
        combined = with_temp_output_dir(() -> ADRIA.combine_results(a, b))
        n_a, n_b = size(a.inputs, 1), size(b.inputs, 1)
        @test size(combined.inputs, 1) == n_a + n_b

        # Union of the columns, in order of first appearance, all Float64: the columns of
        # the first source, then those only the second has, in its order
        in_dropped = in(string.(dropped))
        expected = if a === rs_a
            [names(rs_a.inputs); filter(in_dropped, names(TEST_RS.inputs))]
        else
            names(TEST_RS.inputs)
        end
        @test names(combined.inputs) == expected
        @test all(eltype(c) == Float64 for c in eachcol(combined.inputs))

        # Missing factors are NaN for the source lacking them, values kept elsewhere
        lacking = a === rs_a ? (1:n_a) : ((n_a + 1):(n_a + n_b))
        has = a === rs_a ? ((n_a + 1):(n_a + n_b)) : (1:n_a)
        for col in dropped
            @test all(isnan, combined.inputs[lacking, col])
            @test combined.inputs[has, col] == TEST_RS.inputs[:, col]
        end

        # One model specification per source, in order, and the source of each scenario
        @test length(ADRIA.model_specs(combined)) == 2
        @test ADRIA.scenario_sources(combined) == [fill(1, n_a); fill(2, n_b)]
    end

    # Two sources with different specifications: one lacks the dropped factors, and the
    # two disagree on the default value of a factor they share
    shared = TEST_RS.model_spec.fieldname[1]
    spec_a = filter(r -> !(r.fieldname in dropped), TEST_RS.model_spec)
    spec_a[spec_a.fieldname .== shared, :val] .= 99.0
    rs_b = let r = rs_a
        ADRIA.@set! r.model_spec = spec_a
        ADRIA.@set! r.model_specs = [spec_a]
        r
    end

    combined = with_temp_output_dir(() -> ADRIA.combine_results(rs_b, TEST_RS))
    ms = ADRIA.model_spec(combined)
    val_of(spec, f) = only(spec[spec.fieldname .== f, :val])

    # The union lists every factor of either source once, in order of first appearance
    @test allunique(zip(ms.component, ms.fieldname))
    @test ms.fieldname ==
        [spec_a.fieldname; filter(in(dropped), TEST_RS.model_spec.fieldname)]

    # The first source wins for a shared factor; a factor only the second has is taken
    # from it
    @test val_of(ms, shared) == 99.0
    @test val_of(ms, :cots_on) == val_of(TEST_RS.model_spec, :cots_on)
    flipped = with_temp_output_dir(() -> ADRIA.combine_results(TEST_RS, rs_b))
    @test val_of(ADRIA.model_spec(flipped), shared) == val_of(TEST_RS.model_spec, shared)

    # Each source keeps its own specification
    specs = ADRIA.model_specs(combined)
    @test specs[1].fieldname == spec_a.fieldname && val_of(specs[1], shared) == 99.0
    @test specs[2].fieldname == TEST_RS.model_spec.fieldname
    @test val_of(specs[2], shared) == val_of(TEST_RS.model_spec, shared)

    # Sets that are not merged have a single source
    @test length(ADRIA.model_specs(TEST_RS)) == 1
    @test ADRIA.scenario_sources(TEST_RS) == fill(1, size(TEST_RS.inputs, 1))
end

@testset "Union spec: bounds span the sources, a varying factor is not constant" begin
    ms = TEST_RS.model_spec
    row_of(spec, f) = only(eachrow(spec[spec.fieldname .== f, :]))
    col_vals(rs, f) = rs.inputs[!, string(f)]

    # A factor that varies in TEST_RS, and one that is constant in it
    in_table(f) = string(f) in names(TEST_RS.inputs)
    spread(f) = in_table(f) && length(unique(col_vals(TEST_RS, f))) > 1
    varying = first(filter(spread, ms.fieldname[.!ms.is_constant]))
    constant = :cots_on
    @test row_of(ms, constant).is_constant

    # A source in which `varying` is held at a value outside its bounds in TEST_RS and
    # `constant` is held at another value
    v = row_of(ms, varying).upper_bound + 1.0
    c = row_of(ms, constant).upper_bound + 1.0
    spec_b = copy(ms)
    for (f, x) in ((varying, v), (constant, c))
        i = findfirst(==(f), spec_b.fieldname)
        spec_b[i, :val] = spec_b[i, :lower_bound] = spec_b[i, :upper_bound] = x
        spec_b[i, :is_constant] = true
    end
    inputs_b = copy(TEST_RS.inputs)
    inputs_b[!, string(varying)] .= v
    inputs_b[!, string(constant)] .= c
    rs_b = let r = TEST_RS
        ADRIA.@set! r.inputs = inputs_b
        ADRIA.@set! r.model_spec = spec_b
        ADRIA.@set! r.model_specs = [spec_b]
        r
    end

    for (x, y) in ((rs_b, TEST_RS), (TEST_RS, rs_b))
        merged = with_temp_output_dir(() -> ADRIA.combine_results(x, y))
        ums = ADRIA.model_spec(merged)
        for (f, held) in ((varying, v), (constant, c))
            r = row_of(ums, f)
            @test !r.is_constant
            @test r.lower_bound == min(row_of(ms, f).lower_bound, held)
            @test r.upper_bound == max(row_of(ms, f).upper_bound, held)
        end

        # A factor constant at the same value in both sources stays constant
        others = filter(f -> f != constant && in_table(f), ms.fieldname[ms.is_constant])
        same = first(others)
        @test row_of(ums, same).is_constant
    end

    # Bounds equal in both sources, but the merged column holds two values: not constant;
    # NaN (a source lacking the factor) is not a value
    DF = ADRIA.DataFrames.DataFrame
    spec(lb, ub) = DF(;
        component=["C"], fieldname=[:f], val=[lb], lower_bound=[lb], upper_bound=[ub],
        is_constant=[lb == ub]
    )
    union_row(specs, col) = only(eachrow(ADRIA._union_model_specs(specs, DF(; f=col))))
    @test !union_row([spec(1.0, 1.0), spec(1.0, 1.0)], [1.0, 2.0]).is_constant
    @test union_row([spec(1.0, 1.0), spec(1.0, 1.0)], [1.0, NaN]).is_constant
    @test !union_row([spec(1.0, 1.0), spec(2.0, 2.0)], [1.0, 1.0]).is_constant

    # Specifications without bounds, in one source or all, merge without error
    no_bounds = DF(; component=["C"], fieldname=[:f], val=[3.0])
    r = union_row([spec(1.0, 2.0), no_bounds], [1.0, 3.0])
    @test (r.val, r.lower_bound, r.upper_bound, r.is_constant) == (1.0, 1.0, 2.0, false)
    r = union_row([no_bounds, spec(1.0, 2.0)], [3.0, 1.0])
    @test r.val == 3.0 && (r.lower_bound, r.upper_bound) == (1.0, 2.0)
    @test ADRIA._union_model_specs([no_bounds, no_bounds], DF(; f=[3.0, 3.0])) ==
        no_bounds
end

@testset "Combining rejects result sets whose sources do not match their scenarios" begin
    n = size(TEST_RS.inputs, 1)
    combine(a, b) = with_temp_output_dir(() -> ADRIA.combine_results(a, b))

    too_many = let r = TEST_RS
        ADRIA.@set! r.source_sizes = [n + 1]
        r
    end
    @test_throws DimensionMismatch combine(TEST_RS, too_many)
    @test_throws "sources hold $(n + 1) scenarios" combine(too_many, TEST_RS)

    # One model specification and one provenance record per source
    split_sizes = let r = TEST_RS
        ADRIA.@set! r.source_sizes = [1, n - 1]
        r
    end
    @test_throws DimensionMismatch combine(TEST_RS, split_sizes)
    @test_throws "2 sources but 1 model specifications" combine(split_sizes, TEST_RS)

    # A provenance record too many, with sizes and specifications consistent
    extra_record = let r = TEST_RS
        ADRIA.@set! r.provenance = Any[nothing, nothing]
        r
    end
    @test_throws DimensionMismatch combine(TEST_RS, extra_record)
    @test_throws "2 provenance records" combine(extra_record, TEST_RS)
end

@testset "Loading rejects inconsistent merged stores" begin
    n = size(TEST_RS.inputs, 1)

    """
        edited_store(edit!)

    Copy of the fixture store with its scenario-table attributes edited by `edit!`.
    """
    function edited_store(edit!)
        dir = mktempdir()
        cp(TEST_RS.env_layer_md.dpkg_path, dir; force=true)
        fn = joinpath(dir, ADRIA.INPUTS, ".zattrs")
        attrs = ADRIA.JSON.parsefile(fn)
        edit!(attrs)
        write(fn, ADRIA.JSON.json(attrs))
        return dir
    end

    # Source sizes must add up to the number of scenarios
    dir = edited_store(a -> a["source_sizes"] = [n + 1])
    @test_throws DimensionMismatch ADRIA.load_results(dir)
    @test_throws "inconsistent" ADRIA.load_results(dir)

    # A multi-source store needs one model specification file per source
    dir = edited_store(a -> (a["source_sizes"]=[1, n - 1]; delete!(a, "provenance")))
    @test_throws ArgumentError ADRIA.load_results(dir)
    @test_throws "model_spec_source_1.csv is missing" ADRIA.load_results(dir)

    # Provenance records must match the sources (absent provenance is fine)
    dir = edited_store(a -> a["source_sizes"] = [1, n - 1])
    @test_throws ArgumentError ADRIA.load_results(dir)
    @test_throws "provenance records" ADRIA.load_results(dir)
end

@testset "Provenance is stored, loaded, merged and optional" begin
    prov = ADRIA.provenance(TEST_RS)
    @test length(prov) == 1
    @test prov[1]["adria_version"] == "v" * string(pkgversion(ADRIA))
    @test prov[1]["models"] == []  # COTS off, nothing attached
    @test prov[1]["overrides"] == []

    # A merged set keeps one record per source, in source order (read back from its store)
    merged = with_temp_output_dir(() -> ADRIA.combine_results(TEST_RS, TEST_RS))
    @test ADRIA.provenance(merged) == [prov[1], prov[1]]
    @test length(ADRIA.provenance(merged)) == length(merged.source_sizes)

    # Records of different sources stay distinct and re-merging flattens them
    other = Dict{String,Any}(prov[1]..., "adria_version" => "v0.0.0")
    rs_other = let r = TEST_RS
        ADRIA.@set! r.provenance = Any[other]
        r
    end
    remerged = with_temp_output_dir(() -> ADRIA.combine_results(rs_other, merged))
    @test ADRIA.provenance(remerged) == [other, prov[1], prov[1]]
    @test length(ADRIA.provenance(remerged)) == length(ADRIA.model_specs(remerged))

    # A store written without provenance loads with `nothing`; merging it keeps alignment
    old_dir = mktempdir()
    cp(TEST_RS.env_layer_md.dpkg_path, old_dir; force=true)
    fn = joinpath(old_dir, ADRIA.INPUTS, ".zattrs")
    attrs = ADRIA.JSON.parsefile(fn)
    delete!(attrs, "provenance")
    write(fn, ADRIA.JSON.json(attrs))
    old_rs = ADRIA.load_results(old_dir)
    @test ADRIA.provenance(old_rs) == [nothing]

    mixed = with_temp_output_dir(() -> ADRIA.combine_results(old_rs, TEST_RS))
    @test ADRIA.provenance(mixed) == [nothing, prov[1]]

    # Merging only sources without provenance stores none and loads as `nothing` per source
    none = with_temp_output_dir(() -> ADRIA.combine_results(old_rs, old_rs))
    @test ADRIA.provenance(none) == [nothing, nothing]
end

@testset "Submodel logs" begin
    using ADRIA: Zarr, LogSpec

    """Copy of the fixture store with its log group changed by `edit!(log_fn, n_scens)`."""
    function store_with_logs(edit!)
        dir = mktempdir()
        cp(TEST_RS.env_layer_md.dpkg_path, dir; force=true)
        edit!(joinpath(dir, ADRIA.LOG_GRP), size(TEST_RS.inputs, 1))
        return ADRIA.load_results(dir)
    end

    rs_ids = TEST_RS.loc_ids
    NO_LOGS = Pair{Symbol,Vector{LogSpec}}[]
    tf = length(ADRIA.timesteps(TEST_DOM))
    xspec = LogSpec(:x, [:timesteps], [tf]; units="u")

    """A store whose scenarios hold `val` in log `x` of each submodel in `groups`."""
    function with_x(groups, val)
        return store_with_logs() do log_fn, n
            logs = ADRIA._setup_submodel_logs(
                log_fn, [g => [xspec] for g in groups], rs_ids, n, n
            )
            for g in groups
                logs[g].x[:, :] = fill(Float32(val), tf, n)
            end
        end
    end

    """`combine_results` of `x` and `y`, written to a temporary output directory."""
    combine(x, y) = with_temp_output_dir(() -> ADRIA.combine_results(x, y))

    # Stores and merges shared by the merge tests below (each takes a copy of the fixture)
    n = size(TEST_RS.inputs, 1)
    src_a1, src_a5, src_b2 = with_x([:sub_a], 1), with_x([:sub_a], 5), with_x([:sub_b], 2)
    merged_ab = combine(src_a1, src_b2)

    @testset "no submodel logs without COTS" begin
        @test TEST_RS.logs == NamedTuple()
        # The deprecated accessors warn and give nothing when absent (their values cannot
        # be read under --depwarn=error, where the warning is an error)
        @test_deprecated TEST_RS.cots_pop_log
        @test_deprecated TEST_RS.cots_condition_log
        if Base.JLOptions().depwarn != 2
            @test isnothing(TEST_RS.cots_pop_log)
            @test isnothing(TEST_RS.cots_condition_log)
        end
        @test_throws ArgumentError TEST_RS.not_a_log
        @test propertynames(TEST_RS) == fieldnames(typeof(TEST_RS))
    end

    @testset "COTS logs are declared as LogSpecs and grouped under `cots`" begin
        specs = ADRIA.cots_log_specs(TEST_DOM)
        @test [s.name for s in specs] == [:cots_pop_log, :cots_condition_log]
        tf = length(ADRIA.timesteps(TEST_DOM))
        n_locs = TEST_DOM.coral_growth.n_locs
        @test specs[1].size == [tf, 3, n_locs]
        @test specs[2].size == [tf, n_locs]

        # The mapping to the names in the result set (`rs.logs.cots.pop`)
        @test ADRIA.COTS_LOG_GROUP == :cots
        gspecs = ADRIA.cots_group_log_specs(TEST_DOM)
        @test [s.name for s in gspecs] == [:pop, :condition]
        for (g, s) in zip(gspecs, specs)
            @test (g.dims, g.size, g.labels, g.units, g.description) ==
                (s.dims, s.size, s.labels, s.units, s.description)
        end

        ids = ["a", "b"]
        pop = ADRIA._log_attrs(specs[1], ids)
        @test pop == Dict{Symbol,Any}(
            :structure => ("timesteps", "age_class", "locations", "scenarios"),
            :age_classes => ["recruits", "juveniles", "adults"],
            :units => "COTS ha^-1",
            :unique_loc_ids => ids
        )
        @test ADRIA._log_attrs(specs[2], ids) == Dict{Symbol,Any}(
            :structure => ("timesteps", "locations", "scenarios"),
            :unique_loc_ids => ids
        )
    end

    @testset "LogSpec round-trips through a store, merged attrs equal the source's" begin
        spec = LogSpec(
            :synthetic_log, [:timesteps, :stage, :locations], [2, 2, 3];
            labels=Dict(:stage => ["a", "b"]), units="u", description="d"
        )
        ids = ["x", "y", "z"]

        dir = mktempdir()
        Zarr.zgroup(Zarr.DirectoryStore(dir), ADRIA.LOG_GRP)
        log_fn = joinpath(dir, ADRIA.LOG_GRP)
        logs = ADRIA._setup_submodel_logs(
            log_fn, [:sub => [spec]], ids, 4, 2; fill=Float32(NaN)
        )
        @test keys(logs) == (:sub,) && keys(logs.sub) == (:synthetic_log,)

        # Layout: a subgroup of the log group per submodel
        grp = Zarr.zopen(log_fn)
        @test collect(keys(grp.groups)) == ["sub"] && isempty(grp.arrays)
        @test collect(keys(grp.groups["sub"].arrays)) == ["synthetic_log"]
        spec2 = ADRIA._logspec_from_store(grp.groups["sub"], :synthetic_log)
        @test all(getfield(spec2, f) == getfield(spec, f) for f in fieldnames(typeof(spec)))

        # A merged store builds its arrays from the rebuilt spec: same attributes
        dir2 = mktempdir()
        Zarr.zgroup(Zarr.DirectoryStore(dir2), ADRIA.LOG_GRP)
        log_fn2 = joinpath(dir2, ADRIA.LOG_GRP)
        ADRIA._setup_submodel_logs(log_fn2, [:sub => [spec2]], ids, 8, 2; fill=Float32(NaN))

        a1 = Zarr.zopen(log_fn).groups["sub"]["synthetic_log"].attrs
        a2 = Zarr.zopen(log_fn2).groups["sub"]["synthetic_log"].attrs
        @test a1 == a2
    end

    @testset "submodel and log names are checked" begin
        dir = mktempdir()
        Zarr.zgroup(Zarr.DirectoryStore(dir), ADRIA.LOG_GRP)
        log_fn = joinpath(dir, ADRIA.LOG_GRP)
        setup(pairs) = ADRIA._setup_submodel_logs(log_fn, pairs, ["a"], 2, 1)
        @test_throws "`core` is reserved" setup([:core => [xspec]])
        @test_throws "`rankings` is reserved" setup([:rankings => [xspec]])
        @test_throws ArgumentError setup([Symbol("a b") => [xspec]])
        @test_throws ArgumentError setup([:a => [xspec], :a => [xspec]])
        @test_throws "same name" setup([:a => [xspec, xspec]])

        # The same log name in different submodels is fine
        logs = setup([:a => [xspec], :b => [xspec]])
        @test keys(logs) == (:a, :b) && keys(logs.a) == keys(logs.b) == (:x,)
    end

    @testset "load_results discovers submodels, core arrays stay flat" begin
        rs = with_x([:sub_b, :sub_a], 3)
        @test keys(rs.logs) == (:sub_a, :sub_b)
        @test keys(rs.logs.sub_a) == (:x,)
        @test size(rs.logs.sub_a.x) == (tf, size(TEST_RS.inputs, 1))
        @test all(==(3.0), rs.logs.sub_b.x.data)
        @test_throws ArgumentError rs.x

        # Core logs are untouched
        @test rs.ranks == TEST_RS.ranks

        # A subgroup array needs a `structure` attribute
        @test_throws ArgumentError store_with_logs() do log_fn, n
            Zarr.zgroup(Zarr.DirectoryStore(log_fn), "bad")
            Zarr.zcreate(Float32, 2, 2; name="y", path=joinpath(log_fn, "bad"))
        end
    end

    @testset "stores of the old flat layout load as logs.cots.pop / .condition" begin
        specs = ADRIA.cots_log_specs(TEST_DOM)

        """Flat arrays under the LogSpec names, as written before logs were grouped."""
        function write_flat(log_fn, n; specs=specs)
            for (k, spec) in enumerate(specs)
                arr = Zarr.zcreate(
                    Float32, spec.size..., n; name=String(spec.name), path=log_fn,
                    attrs=ADRIA._log_attrs(spec, rs_ids)
                )
                vals = fill(Float32(k), spec.size..., n)
                arr[fill(Colon(), ndims(vals))...] = vals
            end
        end

        old = store_with_logs(write_flat)
        @test keys(old.logs) == (:cots,)
        @test keys(old.logs.cots) == (:condition, :pop)
        @test size(old.logs.cots.pop) == (specs[1].size..., n)
        @test size(old.logs.cots.pop, :age_class) == 3
        @test all(==(1.0), old.logs.cots.pop.data)
        @test all(==(2.0), old.logs.cots.condition.data)

        # The deprecated accessors see them
        @test_deprecated old.cots_pop_log
        @test_deprecated old.cots_condition_log
        if Base.JLOptions().depwarn != 2
            @test old.cots_pop_log === old.logs.cots.pop
            @test old.cots_condition_log === old.logs.cots.condition
        end

        # A store holding only one of the two flat COTS arrays loads that one
        pop_only = store_with_logs((fn, n) -> write_flat(fn, n; specs=specs[1:1]))
        @test keys(pop_only.logs) == (:cots,) && keys(pop_only.logs.cots) == (:pop,)
        @test all(==(1.0), pop_only.logs.cots.pop.data)
        @test !haskey(pop_only.logs.cots, :condition)

        # A flat array that is not a known log is not mistaken for one, nor are core logs
        custom = store_with_logs(
            (fn, n) -> Zarr.zcreate(Float32, 2, 2; name="zz_custom_log", path=fn)
        )
        @test keys(custom.logs) == ()

        # A log stored both flat and in its subgroup is rejected, naming both (a bare log
        # group, so no copy of the fixture store is needed)
        dir = mktempdir()
        Zarr.zgroup(Zarr.DirectoryStore(dir), ADRIA.LOG_GRP)
        log_fn = joinpath(dir, ADRIA.LOG_GRP)
        write_flat(log_fn, 1; specs=specs[1:1])
        ADRIA._setup_submodel_logs(
            log_fn, [:cots => ADRIA.cots_group_log_specs(TEST_DOM)[1:1]], rs_ids, 1, 1
        )
        both = () -> ADRIA._load_submodel_logs(Zarr.zopen(log_fn))
        @test_throws ArgumentError both()
        @test_throws "as `cots/pop` and as the flat array `cots_pop_log`" both()

        # Merging an old store with a new-layout one keeps the log values, in the new layout
        new = store_with_logs() do log_fn, n
            logs = ADRIA._setup_submodel_logs(
                log_fn, [:cots => ADRIA.cots_group_log_specs(TEST_DOM)], rs_ids, n, n
            )
            logs.cots.pop[:, :, :, :] = fill(Float32(7), specs[1].size..., n)
        end
        merged = combine(old, new)
        @test keys(merged.logs.cots) == (:condition, :pop)
        @test all(==(1.0), merged.logs.cots.pop[scenarios = 1:n].data)
        @test all(==(7.0), merged.logs.cots.pop[scenarios = (n + 1):(2n)].data)
        @test all(==(2.0), merged.logs.cots.condition[scenarios = 1:n].data)
        @test all(==(0.0), merged.logs.cots.condition[scenarios = (n + 1):(2n)].data)

        merged_pop = Zarr.zopen(joinpath(merged.env_layer_md.dpkg_path, ADRIA.LOG_GRP))
        @test merged_pop.groups["cots"]["pop"].attrs["units"] == "COTS ha^-1"
    end

    @testset "merge: same log name in different submodels stays separate" begin
        @test keys(merged_ab.logs) == (:sub_a, :sub_b)
        @test all(==(1.0), merged_ab.logs.sub_a.x[scenarios = 1:n].data)
        @test all(isnan, merged_ab.logs.sub_a.x[scenarios = (n + 1):(2n)].data)
        @test all(isnan, merged_ab.logs.sub_b.x[scenarios = 1:n].data)
        @test all(==(2.0), merged_ab.logs.sub_b.x[scenarios = (n + 1):(2n)].data)

        # Same submodels in all sources: no NaN
        m2 = combine(src_a1, src_a5)
        @test keys(m2.logs) == (:sub_a,)
        @test m2.logs.sub_a.x.data == hcat(fill(1.0, tf, n), fill(5.0, tf, n))
    end

    @testset "merge: a source lacking a submodel is filled with NaN" begin
        cases = ((src_a1, TEST_RS, 1:n, (n + 1):(2n)), (TEST_RS, src_a1, (n + 1):(2n), 1:n))
        for (x, y, has, lacks) in cases
            m = combine(x, y)
            @test keys(m.logs) == (:sub_a,)
            @test all(==(1.0), m.logs.sub_a.x[scenarios = has].data)
            @test all(isnan, m.logs.sub_a.x[scenarios = lacks].data)
        end
    end

    @testset "merge: a source lacking one log of a submodel is filled with NaN for it" begin
        yspec = LogSpec(:y, [:timesteps], [tf]; units="u")
        xy = store_with_logs() do log_fn, n
            specs = [:sub_a => [xspec, yspec]]
            logs = ADRIA._setup_submodel_logs(log_fn, specs, rs_ids, n, n)
            logs.sub_a.x[:, :] = fill(Float32(3), tf, n)
            logs.sub_a.y[:, :] = fill(Float32(4), tf, n)
        end

        m = combine(src_a1, xy)
        @test keys(m.logs) == (:sub_a,) && keys(m.logs.sub_a) == (:x, :y)
        @test m.logs.sub_a.x.data == hcat(fill(1.0, tf, n), fill(3.0, tf, n))
        @test all(isnan, m.logs.sub_a.y[scenarios = 1:n].data)
        @test all(==(4.0), m.logs.sub_a.y[scenarios = (n + 1):(2n)].data)
    end

    @testset "re-merge: logs, NaN fill and sources carry through" begin
        # Merging the merged set with a third source: three sources, aligned records
        abc = combine(merged_ab, src_a5)
        @test abc.source_sizes == [n, n, n]
        @test length(ADRIA.model_specs(abc)) == 3
        @test length(ADRIA.provenance(abc)) == 3
        @test ADRIA.scenario_sources(abc) == [fill(1, n); fill(2, n); fill(3, n)]

        # Values and the NaN fill of the first merge are kept; the new source lacks sub_b
        @test keys(abc.logs) == (:sub_a, :sub_b)
        x_a, x_b = abc.logs.sub_a.x, abc.logs.sub_b.x
        @test all(==(1.0), x_a[scenarios = 1:n].data)
        @test all(isnan, x_a[scenarios = (n + 1):(2n)].data)
        @test all(==(5.0), x_a[scenarios = (2n + 1):(3n)].data)
        @test all(isnan, x_b[scenarios = 1:n].data)
        @test all(==(2.0), x_b[scenarios = (n + 1):(2n)].data)
        @test all(isnan, x_b[scenarios = (2n + 1):(3n)].data)
    end

    @testset "merge: a log must have the same shape, labels and units in all sources" begin
        """A store with log `x` of `sub_a` declared by `spec` instead, filled with 2."""
        function with_spec(spec)
            return store_with_logs() do log_fn, n
                logs = ADRIA._setup_submodel_logs(log_fn, [:sub_a => [spec]], rs_ids, n, n)
                arr = logs.sub_a.x
                arr[fill(Colon(), ndims(arr))...] = fill(Float32(2), size(arr))
            end
        end

        wider = with_spec(LogSpec(:x, [:timesteps, :stage], [tf, 2]; units="u"))
        @test_throws DimensionMismatch combine(src_a1, wider)
        @test_throws "sub_a.x" combine(wider, src_a1)

        other_units = with_spec(LogSpec(:x, [:timesteps], [tf]; units="v"))
        @test_throws ArgumentError combine(src_a1, other_units)
        years = Dict(:timesteps => string.(1:tf))
        labelled = with_spec(LogSpec(:x, [:timesteps], [tf]; units="u", labels=years))
        @test_throws ArgumentError combine(src_a1, labelled)
    end

    @testset "per-scenario byte budget keeps the sum and adds submodel logs" begin
        tf = length(ADRIA.timesteps(TEST_DOM))
        nl = TEST_DOM.coral_growth.n_locs
        ng = TEST_DOM.coral_growth.n_groups
        n_iv = length(ADRIA.interventions())
        n_met = length(ADRIA.LOC_METRIC_NAMES)

        base = ADRIA._log_bytes_per_scen(TEST_DOM, NO_LOGS)
        @test base.total ==
            4 * tf * (nl * n_met + ng + nl * 2 + nl * n_iv + ng * nl * 2)

        cots = [ADRIA.COTS_LOG_GROUP => ADRIA.cots_group_log_specs(TEST_DOM)]
        with = ADRIA._log_bytes_per_scen(TEST_DOM, cots)
        @test with.total == base.total + 4 * tf * (3 * nl + nl)
        @test base.largest >= 4 * tf * ng * nl  # includes optional logs when enabled
        @test with.largest >= base.largest

        # A large submodel log is in the sum, so it limits the chunk through `total`
        big = [:big => [LogSpec(:big, [:timesteps, :x], [tf, 10_000_000])]]
        b = ADRIA._log_bytes_per_scen(TEST_DOM, big)
        @test b.largest == 4 * tf * 10_000_000
        @test ADRIA._max_chunk_by_size_for(TEST_DOM, big) ==
            max(1, floor(Int, 1_800_000_000 / b.total))
        @test ADRIA._max_chunk_by_size_for(TEST_DOM, big) <
            ADRIA._max_chunk_by_size_for(TEST_DOM, NO_LOGS)
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
        dom_cap(sub) = ADRIA._max_chunk_by_size_for(TEST_DOM, sub)
        off = dom_cap(NO_LOGS)
        withenv("ADRIA_LOG_COVER" => "true", "ADRIA_LOG_DHW_TOLS" => "true") do
            @test dom_cap(NO_LOGS) == off
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
