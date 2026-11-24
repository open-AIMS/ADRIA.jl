if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

const _SPEC_COMPONENTS = collect(ADRIA.spec_components())

@testset "Extracting component parameters" begin
    # Run full example to make sure nothing errors
    ADRIA.setup()  # Load and apply configuration options

    # Use a temporary directory for result location
    ENV["ADRIA_OUTPUT_DIR"] = mktempdir()

    x = ADRIA.component_params(ADRIA_DOM_45, Intervention)
    @test size(x, 1) > 0

    x = ADRIA.component_params(ADRIA_DOM_45, ADRIA.Coral)
    @test size(x, 1) > 0

    x = ADRIA.component_params(ADRIA_DOM_45, ADRIA.CAqCriteriaWeights)
    @test size(x, 1) > 0

    x = ADRIA.component_params(ADRIA_DOM_45, ADRIA.FogCriteriaWeights)
    @test size(x, 1) > 0

    # Every component is present, rows are in component order and the slices tile the spec
    ms = ADRIA.model_spec(ADRIA_DOM_45)
    parts = [ADRIA.component_params(ADRIA_DOM_45, c) for c in _SPEC_COMPONENTS]
    @test all(DataFrames.nrow.(parts) .> 0)
    @test reduce(vcat, parts) == ms
    @test ADRIA.component_params(ADRIA_DOM_45, [ADRIA.Coral, ADRIA.CotsParams]) ==
        vcat(parts[7], parts[11])
    @test allunique(ms.fieldname)

    # The spec on a DataFrame filters in the same way
    @test ADRIA.component_params(ms, ADRIA.Intervention) == parts[2]
end

@testset "Model spec layout" begin
    ms = ADRIA.model_spec(ADRIA_DOM_45)
    @test names(ms) == [
        "component", "fieldname", "val", "ptype", "dist", "dist_params",
        "default_dist_params", "direction", "lower_bound", "upper_bound", "is_constant",
        "name", "description"
    ]
    @test eltype(ms.val) == Float64
    @test eltype(ms.dist_params) <: Tuple
    @test all(ms.lower_bound .<= ms.upper_bound)
    @test ms.is_constant == (ms.lower_bound .== ms.upper_bound)
    @test all(isnothing, ms[ms.component .== "Coral", :direction])
    @test all(
        d -> d in (minimum, maximum),
        ADRIA.component_params(ms, ADRIA.CAqCriteriaWeights).direction
    )

    # `model_spec` hands out a copy: editing it does not change the domain
    ms.val .= 0.0
    @test ADRIA.model_spec(ADRIA_DOM_45).val != ms.val
end

@testset "Table rows follow the component definitions" begin
    # The table is built from the component definitions, so a default instance of every
    # component must agree with its rows (growth caps only lower scale-factor bounds).
    ADRIA.create_coral_struct()  # reset bounds so the result is test-order independent
    ADRIA.create_growth_acceleration_struct()
    dom = ADRIA_DOM_45
    n_dhw = length(ADRIA.component_params(dom, ADRIA.EnvironmentalLayer).fieldname)
    @test n_dhw == 3

    insts = Dict(
        ADRIA.Intervention => ADRIA.Intervention(),
        ADRIA.CAqCriteriaWeights => ADRIA.CAqCriteriaWeights(),
        ADRIA.FogCriteriaWeights => ADRIA.FogCriteriaWeights(),
        ADRIA.LvMCriteriaWeights => ADRIA.LvMCriteriaWeights(),
        ADRIA.DepthThresholds => ADRIA.DepthThresholds(),
        ADRIA.Coral => ADRIA.Coral(),
        ADRIA.GrowthAcceleration => ADRIA.GrowthAcceleration(),
        ADRIA.DepthAttenuation => ADRIA.DepthAttenuation(),
        ADRIA.LowDHWBleaching => ADRIA.LowDHWBleaching(),
        ADRIA.CotsParams => ADRIA.CotsParams()
    )
    for (comp, inst) in insts
        rows = ADRIA.component_params(dom, comp)
        @test Symbol.(rows.fieldname) == collect(fieldnames(typeof(inst)))
        for (i, fn) in enumerate(fieldnames(typeof(inst)))
            p = getfield(inst, fn)
            # the test domain has no MCB data: its iv_MCB_* rows are the defaults too
            @test rows.val[i] == p.val
            @test rows.ptype[i] == p.ptype
            @test rows.dist[i] === p.dist
            @test rows.name[i] == p.name
            @test rows.description[i] == p.description
            if !startswith(string(fn), "linear_extension_scale")
                @test collect(rows.dist_params[i]) == p.dist_params
                @test collect(rows.default_dist_params[i]) == p.default_dist_params
            end
        end
    end
end

@testset "create_coral_struct bounds reach subsequently loaded domains" begin
    try
        ADRIA.create_coral_struct((0.5, 1.5))
        ADRIA.create_growth_acceleration_struct(0.5)
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
        row = ADRIA.component_params(dom, ADRIA.Coral)
        i = findfirst(==(:tabular_Acropora_1_1_linear_extension), row.fieldname)
        v = row.val[i]
        @test row.dist_params[i][1] ≈ v * 0.5
        @test row.dist_params[i][2] ≈ v * 1.5
        ga = ADRIA.component_params(dom, ADRIA.GrowthAcceleration)
        j = findfirst(==(:growth_acceleration_cb_group_1_height), ga.fieldname)
        @test ga.dist_params[j][1] ≈ ga.val[j] * 0.5
        @test ga.dist_params[j][2] ≈ ga.val[j] * 1.5
    finally
        ADRIA.create_coral_struct()
        ADRIA.create_growth_acceleration_struct()
    end
end

@testset "create_coral_struct overrides reach subsequently loaded domains" begin
    try
        ADRIA.create_coral_struct(
            (0.9, 1.1);
            overrides=Dict("tabular_Acropora_1_1_linear_extension" => 0.05)
        )
        ADRIA.create_growth_acceleration_struct(
            0.1; overrides=Dict("growth_acceleration_cb_group_1_height" => 0.5)
        )
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
        row = ADRIA.component_params(dom, ADRIA.Coral)
        i = findfirst(==(:tabular_Acropora_1_1_linear_extension), row.fieldname)
        @test row.val[i] == 0.05
        @test collect(row.dist_params[i]) ≈ [0.045, 0.055, 0.05]
        @test row.val[i] ==
            Base.invokelatest(ADRIA.Coral).tabular_Acropora_1_1_linear_extension.val
        ga = ADRIA.component_params(dom, ADRIA.GrowthAcceleration)
        j = findfirst(==(:growth_acceleration_cb_group_1_height), ga.fieldname)
        @test ga.val[j] == 0.5
        @test collect(ga.dist_params[j]) ≈ [0.45, 0.55]
    finally
        ADRIA.create_coral_struct()
        ADRIA.create_growth_acceleration_struct()
    end
    # reset: defaults again
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
    row = ADRIA.component_params(dom, ADRIA.Coral)
    i = findfirst(==(:tabular_Acropora_1_1_linear_extension), row.fieldname)
    @test row.val[i] != 0.05
end

@testset "Parameter table components" begin
    # Every ADRIA struct that holds `Param` factors is a component of the table, so a new
    # component cannot be left out of `spec_components()` unnoticed.
    param_structs = Set{DataType}()
    for m in (ADRIA, ADRIA.decision), n in names(m; all=true)
        isdefined(m, n) || continue
        T = getfield(m, n)
        if T isa DataType && isconcretetype(T) && T.name.module === m &&
            any(ft -> ft <: ADRIA.Param, fieldtypes(T))
            push!(param_structs, T)
        end
    end
    @test param_structs == Set(ADRIA.spec_components())

    ms = ADRIA.model_spec(ADRIA_DOM_45)
    @test unique(ms.component) == collect(ADRIA._component_name.(ADRIA.spec_components()))
end

@testset "Removed `model` field" begin
    err = try
        ADRIA_DOM_45.model
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("model_spec", err.msg) && occursin("component_params", err.msg)
    @test !hasproperty(ADRIA_DOM_45, :model)
end

@testset "Updating the parameter table" begin
    dom = deepcopy(ADRIA_DOM_45)
    @test !hasproperty(dom, :model)

    # Copies of the store are independent; deepcopy of a Domain does not share the table
    st = copy(dom.spec_store)
    st.tbl.val[1] += 1.0
    @test st.tbl.val[1] != dom.spec_store.tbl.val[1]
    @test dom.spec_store.tbl !== ADRIA_DOM_45.spec_store.tbl

    # update! matches rows by name, so a partial spec only touches its own rows
    before = ADRIA.model_spec(dom)
    part = before[before.fieldname .== :cots_on, :]
    part[!, :val] .= 1.0
    ADRIA.update!(dom, part)
    after = ADRIA.model_spec(dom)
    @test after[after.fieldname .== :cots_on, :val] == [1.0]
    @test after[after.fieldname .!= :cots_on, :] == before[before.fieldname .!= :cots_on, :]
    @test ADRIA.param_table(dom).cots_on == [1.0]

    # Full write-back round trips
    ADRIA.update!(dom, after)
    @test ADRIA.model_spec(dom) == after

    # Bounds are written and the derived columns follow
    ms = ADRIA.model_spec(dom)
    i = findfirst(ms.fieldname .== :guided)
    ms[i, :dist_params] = (0.0, 1.0)
    ADRIA.update!(dom, ms)
    g = ADRIA.model_spec(dom)
    g = g[g.fieldname .== :guided, :]
    @test g.lower_bound == [0.0] && g.upper_bound == [1.0]

    # Unknown factors are an error and leave the table untouched
    snapshot = ADRIA.model_spec(dom)
    bad = snapshot[1:2, :]
    bad[!, :fieldname] = [bad.fieldname[1], :not_a_factor]
    bad[!, :val] .= 99.0
    @test_throws ArgumentError ADRIA.update!(dom, bad)
    @test ADRIA.model_spec(dom) == snapshot
end

@testset "update_params! sets factors by name" begin
    dom = deepcopy(ADRIA_DOM_45)
    ms = ADRIA.model_spec(dom)

    # A row holding only some of the factors, in an order unlike the table's
    names_ = [:iv_CAq_N_TA, :guided, :cots_bh_a, :dhw_scenario]
    row = DataFrames.DataFrame(String.(names_) .=> [500_000.0, 1.0, 2.0, 3.0])
    row[!, :RCP] .= 45
    before = ADRIA.model_spec(dom)
    ADRIA.update_params!(dom, row[1, :])
    after = ADRIA.model_spec(dom)
    for (fn, v) in zip(names_, [500_000.0, 1.0, 2.0, 3.0])
        @test after[after.fieldname .== fn, :val] == [v]
    end
    # every other factor is untouched (the positional `Model` update scrambled them)
    others = after.fieldname .∉ Ref(names_)
    @test after[others, :] == before[others, :]

    # Unknown names are an error, and nothing is written
    row = DataFrames.DataFrame("guided" => [-1.0], "not_a_factor" => [1.0])
    @test_throws ArgumentError ADRIA.update_params!(dom, row[1, :])
    @test ADRIA.model_spec(dom) == after
end

@testset "update_params! is atomic when a value cannot be converted" begin
    dom = deepcopy(ADRIA_DOM_45)
    before = ADRIA.model_spec(dom)
    row = DataFrames.DataFrame("guided" => [1.0], "cots_bh_a" => [missing])
    @test_throws MethodError ADRIA.update_params!(dom, row[1, :])
    @test ADRIA.model_spec(dom) == before
end

@testset "A partial spec is checked against the growth limits" begin
    dom = deepcopy(ADRIA_DOM_45)
    before = ADRIA.model_spec(dom)
    fn = :linear_extension_scale_cb_group_1_tabular_Acropora
    i = findfirst(==(fn), before.fieldname)
    row = before[[i], :]
    row[!, :dist_params] = [(before.dist_params[i][1], 50.0)]
    @test_throws ErrorException ADRIA.update!(dom, row)
    @test ADRIA.model_spec(dom) == before
end

@testset "Only Domains without a `model` field get the removed-field error" begin
    struct _ModelFieldDomain <: ADRIA.Domain
        model::Int
    end
    @test _ModelFieldDomain(3).model == 3
end

@testset "ParamSpec to table rows (`spec_rows`)" begin
    PS = ADRIA.ADRIAinterfaces.ParamSpec
    row(ps; kw...) = only(ADRIA.spec_rows([ps]; kw...))[2]

    # Real, triangular: constant by default (collapsed dist_params), full range kept
    r = row(PS(:x_a, 1.5, (0.5, 2.5); dist=:triang, dist_params=[0.5, 2.5, 1.5], label="A"))
    @test r.ptype == "continuous" && r.dist == ADRIA.TriangularDist && r.name == "A"
    @test r.val == 1.5 && r.dist_params == [1.5, 1.5, 1.5]
    @test r.default_dist_params == [0.5, 2.5, 1.5]

    # ... or sampled when not constant; overrides replace the value only
    r = row(
        PS(
            :x_a, 1.5, (0.5, 2.5); dist=:triang, dist_params=[0.5, 2.5, 1.5], constant=false
        );
        overrides=Dict("x_a" => 2.0))
    @test r.val == 2.0 && r.dist_params == r.default_dist_params == [0.5, 2.5, 1.5]
    @test row(PS(:x_a, 1.5, (0.5, 2.5))).name == "x_a"  # label defaults to the name

    # Uniform (the default dist) and the discrete dists
    r = row(PS(:x_u, 1, (0, 2)))
    @test r.dist == ADRIA.Uniform && r.dist_params == [1.0, 1.0]
    @test r.default_dist_params == [0.0, 2.0]
    r = row(
        PS(:x_i, 6, (3, 10); kind=:integer, dist=:discrete_ordered, dist_params=[3, 10, 1])
    )
    @test r.ptype == "ordered discrete" && r.dist == ADRIA.DiscreteOrderedUniformDist
    @test r.dist_params == [6.0, 6.0, 6.0] && r.default_dist_params == [3.0, 10.0, 1.0]
    r = row(PS(:x_d, 0, (0, 1e9); kind=:integer, dist=:discrete_uniform))
    @test r.ptype == "ordered discrete" && r.dist == ADRIA.DiscreteUniform
    @test r.default_dist_params == [0.0, 1.0e9]

    # Bool and categorical: parameters are the levels
    r = row(PS(:x_b, 0, (0, 1); kind=:bool))
    @test r.ptype == "unordered categorical" && r.dist == ADRIA.CategoricalDistribution
    @test r.dist_params == [0.0, 0.0] && r.default_dist_params == [0.0, 1.0]
    r = row(PS(:x_c, 1, (0, 2); kind=:categorical, levels=[0, 1, 2]))
    @test r.dist_params == [1.0, 1.0, 1.0] && r.default_dist_params == [0.0, 1.0, 2.0]

    # Wrong number of parameters or a kind/dist mismatch is an error
    @test_throws ArgumentError ADRIA.spec_rows([PS(:x_t, 1, (0, 2); dist=:triang)])
    @test_throws ArgumentError ADRIA.spec_rows([PS(:x_t, 1, (0, 2); dist=:categorical)])
    @test_throws ArgumentError ADRIA.spec_rows([PS(:x_t, 1, (0, 2); dist=:nope)])

    # COTS declares its factors as ParamSpecs: the row view and `spec_rows` agree
    fs = ADRIA.COTS_FACTOR_SPECS
    @test length(ADRIA.COTS_PARAM_SPECS) == 43 == length(fs)
    @test first.(ADRIA.spec_rows(ADRIA.COTS_PARAM_SPECS)) ==
        collect(fieldnames(ADRIA.CotsParams)) == collect(keys(fs))
    @test (keys(fs)[1], keys(fs)[2], keys(fs)[24], keys(fs)[41], keys(fs)[43]) ==
        (
        :cots_on, :cots_bh_a, :cots_senescence_age, :cots_larval_conn_seed,
        :cots_initial_state_scale
    )

    # Literal rows from before the ParamSpec declaration: triang, discrete_ordered, bool,
    # categorical and discrete_uniform
    @test fs.cots_bh_a == (val=1.5, ptype="continuous", dist=ADRIA.TriangularDist,
        dist_params=(0.5, 2.5, 1.5), name="COTS Beverton-Holt a",
        description="Beverton-Holt recruitment parameter 'a'.")
    @test (fs.cots_senescence_age.val, fs.cots_senescence_age.ptype,
        fs.cots_senescence_age.dist, fs.cots_senescence_age.dist_params) ==
        (6.0, "ordered discrete", ADRIA.DiscreteOrderedUniformDist, (3.0, 10.0, 1.0))
    @test (fs.cots_on.val, fs.cots_on.ptype, fs.cots_on.dist, fs.cots_on.dist_params,
        fs.cots_on.name) ==
        (
        0.0, "unordered categorical", ADRIA.CategoricalDistribution, (0.0, 1.0), "COTS On"
    )
    @test fs.cots_larval_conn_mode.dist_params == (0.0, 1.0, 2.0)
    @test (fs.cots_larval_conn_seed.ptype, fs.cots_larval_conn_seed.dist,
        fs.cots_larval_conn_seed.dist_params) ==
        ("ordered discrete", ADRIA.DiscreteUniform, (0.0, 1.0e9))

    # Edge cases: equal bounds, non-constant bool/categorical, parameter counts
    r = row(PS(:x_e, 1, (1, 1)))
    @test r.dist_params == [1.0, 1.0] == r.default_dist_params
    r = row(PS(:x_b, 1, (0, 1); kind=:bool, constant=false))
    @test r.dist_params == r.default_dist_params == [0.0, 1.0]
    r = row(PS(:x_c, 1, (0, 2); kind=:categorical, levels=[0, 1, 2], constant=false))
    @test r.dist_params == [0.0, 1.0, 2.0]
    @test_throws ArgumentError ADRIA.spec_rows(
        [PS(:x_o, 1, (0, 5); kind=:integer, dist=:discrete_ordered)]
    )
    @test_throws ArgumentError ADRIA.spec_rows(
        [PS(:x_o, 1, (0, 5); dist=:uniform, dist_params=[0, 2, 5])]
    )
end
