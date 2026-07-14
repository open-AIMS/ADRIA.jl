using ADRIA
using ADRIA.Distributions
using ADRIA.DataFrames

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

@testset "sample" begin
    @info "Sample tests"
    dom = deepcopy(ADRIA_DOM_45)
    num_samples = 32
    scens = ADRIA.sample(dom, num_samples)
    ms = ADRIA.model_spec(dom, scens)
    constant_params = ms.is_constant

    @testset "constant params are constant" begin
        @test all(
            values(scens[1, constant_params]) .== values(scens[end, constant_params])
        ) ||
            "Constant params are not constant!"
    end

    @testset "values are within expected bounds" begin
        lb = values(ms[:, :lower_bound])
        ub = values(ms[:, :upper_bound])

        not_cw_mask =
            ms.component .∉
            [(
                "SeedCriteriaWeights",
                "FogCriteriaWeights",
                "MCCriteriaWeights",
                "decision.DepthThresholds"
            )]

        # Reactive params are zeroed when both strategies are Periodic (not Reactive)
        not_reactive_mask = .!(contains.(string.(ms.fieldname), "reactive"))

        iv_fieldnames = ADRIA.component_params(dom, Intervention).fieldname
        not_iv_mask = .!(
        ms.fieldname .∈ [filter(x -> x ∉ [:guided, :heritability], iv_fieldnames)]
)

        to_test_mask =
            not_cw_mask .& not_reactive_mask .& not_iv_mask
        _lb, _ub = lb[to_test_mask], ub[to_test_mask]

        eco = (ms.component .== "Coral") .& .!(constant_params)

        msg = "Sampled values were not in expected bounds!"
        coral_msg = "Sampled coral values were not in expected bounds!"
        for i = 1:num_samples
            # Filter CriteriaWeights factors
            scen_vals = collect(scens[i, :])
            not_cw_scen_vals = scen_vals[to_test_mask]

            if scens[i, :guided] > 0

                # All params should be within bounds
                cond = _lb .<= not_cw_scen_vals .<= _ub
                @test all(cond) ||
                    "$msg | $(ms[to_test_mask,:][.!(cond), :]) | $(not_cw_scen_vals[.!(cond)])"
                continue
            end

            # When no interventions are used, e.g., for counterfactual or unguided scenarios
            # (guided ∈ [-1, 0]) intervention parameters are set to 0 so only check ecological values
            cond = lb[eco] .<= scen_vals[eco] .<= ub[eco]
            @test all(cond) ||
                "$coral_msg | $(ms[.!(cond), :]) | $(scen_vals[eco][.!(cond)])"

            # Note: Test to ensure all intervention factors are set to 0 is covered by the guided
            # sampling test below
        end
    end

    @testset "guided-dependent params are resolved for plain sample()" begin
        cf_rows = scens.guided .== -1
        unguided_rows = scens.guided .== 0
        non_guided_rows = .!(scens.guided .> 0)

        interv_params = string.(
            ADRIA.component_params(ADRIA.model_spec(dom), Intervention).fieldname
        )
        strategy_params = ["seed_strategy", "fog_strategy", "mc_strategy"]
        non_strategy_interv_params = String[
            ip for ip in interv_params if ip ∉ vcat(["guided"], strategy_params)
        ]
        strategy_interv_params = String[ip for ip in strategy_params if ip ∈ interv_params]

        # CF (guided == -1): all intervention params zeroed, strategies sentineled to -1
        if any(cf_rows)
            @test all(
                all.(==(0), eachrow(scens[cf_rows, non_strategy_interv_params]))
            ) || "CF rows in plain sample() have non-zero intervention params"
            @test all(
                all.(.==(-1), eachrow(scens[cf_rows, strategy_interv_params]))
            ) || "CF rows in plain sample() do not have strategy sentinel -1"
            @test all(scens[cf_rows, :depth_min] .== 0) &&
                  all(scens[cf_rows, :depth_offset] .== 0) ||
                "CF rows in plain sample() have non-zero depth thresholds"
        end

        # Criteria weights only matter when guided > 0
        seed_w = ADRIA.component_params(ms, ADRIA.SeedCriteriaWeights).fieldname
        fog_w = ADRIA.component_params(ms, ADRIA.FogCriteriaWeights).fieldname
        mc_w = ADRIA.component_params(ms, ADRIA.MCCriteriaWeights).fieldname
        criteria_cols = string.(vcat(seed_w, fog_w, mc_w))
        if any(non_guided_rows)
            @test all(
                all.(==(0), eachrow(scens[non_guided_rows, criteria_cols]))
            ) || "Non-guided rows in plain sample() have non-zero criteria weights"
        end

        # plan_horizon only matters when guided != 0
        if any(unguided_rows)
            @test all(scens[unguided_rows, :plan_horizon] .== 0) ||
                "Unguided rows in plain sample() have non-zero plan_horizon"
        end
    end
end

@testset "Targeted sampling" begin
    @testset "Counterfactual sampling" begin
        dom = deepcopy(ADRIA_DOM_45)
        num_samples = 32
        scens = ADRIA.sample_cf(dom, num_samples)

        @test all(scens.guided .== -1) || "Intervention scenarios found"

        # Get Intervention params
        interv_params = string.(
            ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Intervention).fieldname
        )

        # Ensure all interventions are deactivated (ignoring "guided").
        # Strategy params (seed_strategy/fog_strategy/mc_strategy) are set to -1.0
        # as a CF sentinel (see §2.3 change 2 in component_dependency_dag_v3.md);
        # all other intervention params are 0.
        strategy_params = ["seed_strategy", "fog_strategy", "mc_strategy"]
        non_strategy_params = String[
            ip for ip in interv_params
                   if ip != "guided" && ip ∉ strategy_params
        ]
        strategy_interv_params = String[ip for ip in strategy_params if ip ∈ interv_params]

        @test all(all.(==(0), eachrow(scens[:, non_strategy_params]))) ||
            "Non-strategy intervention factors with values != 0 found"
        @test all(all.(.==(-1), eachrow(scens[:, strategy_interv_params]))) ||
            "Strategy params not set to -1 sentinel in CF scenarios"
    end

    @testset "Guided sampling" begin
        @info "Guided sampling"
        dom = deepcopy(ADRIA_DOM_45)
        num_samples = 32
        scens = ADRIA.sample_guided(dom, num_samples)
        ms = ADRIA.model_spec(dom, scens)

        @test all(scens.guided .> 0) || "Non-intervention scenarios found"

        # Get Intervention params
        interv_params = string.(
            ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Intervention).fieldname
        )

        # Ignore guided
        interv_params = String[ip for ip in interv_params if ip != "guided"]

        # Ensure at least one intervention is active
        @test all(any.(>(0), eachrow(scens[:, interv_params]))) ||
            "All intervention factors had values <= 0"

        seed_weights = ADRIA.component_params(ms, ADRIA.SeedCriteriaWeights).fieldname
        fog_weights = ADRIA.component_params(ms, ADRIA.FogCriteriaWeights).fieldname

        @test all(abs.(sum(Matrix(scens[:, seed_weights]); dims=2) .- 1.0) .< 10e-6) ||
            "Some seeding weights are not properly normalized."
        @test all(abs.(sum(Matrix(scens[:, fog_weights]); dims=2) .- 1.0) .< 10e-6) ||
            "Some fogging weights are not properly normalized."
    end

    @testset "Specific Intervention strategy" begin
        dom = deepcopy(ADRIA_DOM_45)

        # :guided is now a strict ordinal {-1, 0, 1} (CF / unguided / guided);
        # `set_factor_bounds!` restricts sampling to plain numeric bounds.
        test_bounds = [(-1.0, -1.0), (0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (-1.0, 1.0)]

        for bounds in test_bounds
            ADRIA.set_factor_bounds!(dom, :guided, bounds)
            scens = ADRIA.sample(dom, 32)

            @test all(bounds[1] .<= scens.guided .<= bounds[2])
        end
    end

    @testset "Specific MCDA method" begin
        dom = deepcopy(ADRIA_DOM_45)

        test_inputs = [
            ("Cocoso",), ("Mairca",),
            ("Cocoso", "Moora", "Piv", "Vikor"),
            ("Cocoso", "Mairca", "Moora", "Piv", "Vikor")
        ]

        for inp in test_inputs
            ADRIA.set_factor_bounds!(dom, :mcda_method, inp)
            # Ensure guided is in the "guided" regime so mcda_method isn't
            # sentineled to 0 by the guided<=0 dependency rule.
            ADRIA.set_factor_bounds!(dom, :guided, (1.0, 1.0))
            scens = ADRIA.sample(dom, 32)

            expected = sort(ADRIA.decision.decision_method_encoding.(collect(inp)))
            @test all(scens.mcda_method .∈ [expected])
        end
    end

    @testset "mcda_method sentineled when not guided" begin
        dom = deepcopy(ADRIA_DOM_45)
        ADRIA.set_factor_bounds!(dom, :guided, (-1.0, 0.0))  # CF + unguided only
        scens = ADRIA.sample(dom, 32)

        @test all(scens.mcda_method .== 0.0)
    end

    @testset "guided/mcda_method split" begin
        dom = deepcopy(ADRIA_DOM_45)
        scens = ADRIA.sample(dom, 32)

        @test :guided in propertynames(scens)
        @test :mcda_method in propertynames(scens)
        @test all(scens.guided .∈ [[-1.0, 0.0, 1.0]])

        # mcda_method must be 0 (sentinel) on every row where guided <= 0,
        # and a valid 1..N index on every row where guided == 1.
        not_guided = scens.guided .<= 0
        is_guided = scens.guided .== 1
        @test all(scens.mcda_method[not_guided] .== 0.0)
        n_methods = length(ADRIA.decision.mcda_methods())
        @test !any(is_guided) ||
            all(1.0 .<= scens.mcda_method[is_guided] .<= n_methods)
    end

    @testset "Unguided sampling" begin
        dom = deepcopy(ADRIA_DOM_45)
        num_samples = 32
        scens = ADRIA.sample_unguided(dom, num_samples)

        @test all(scens.guided .== 0) || "Intervention or counterfactual scenarios found"

        # Get Intervention params
        interv_params = string.(
            ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Intervention).fieldname
        )

        # Ignore guided, planning horizon, and reactive params (which are conditionally zeroed)
        interv_params = String[
            ip for ip in interv_params if
            ip ∉ ["plan_horizon", "guided"] && !contains(ip, "reactive")
        ]

        # Ensure at least one intervention is active
        @test all(any.(>(0), eachrow(scens[:, interv_params]))) ||
            "Some intervention params had values <= 0"
    end

    @testset "Site selection sampling" begin
        dom = deepcopy(ADRIA_DOM_45)
        num_samples = 32
        scens = ADRIA.sample_selection(dom, num_samples)

        @test all(scens.guided .> 0) || "Intervention or counterfactual scenarios found"

        # Get Intervention params
        ms = ADRIA.model_spec(dom)
        target_params = string.(
            ADRIA.component_params(
                ms,
                [
                    ADRIA.EnvironmentalLayer,
                    ADRIA.Intervention,
                    ADRIA.SeedCriteriaWeights,
                    ADRIA.FogCriteriaWeights
                ]
            ).fieldname
        )

        # Ignore guided
        target_params = String[ip for ip in target_params if ip != "guided"]

        # Ensure at least one intervention is active
        @test all(any.(>(0), eachrow(scens[:, target_params]))) ||
            "All target factors had values <= 0"

        # Check that all coral parameters are set to their nominated default values
        coral_params = ADRIA.component_params(ms, ADRIA.Coral).fieldname

        @test all([
            all(scens[:, c] .== ms[ms.fieldname .== c, :val][1]) for c in coral_params
        ]) || "Non-default coral parameter value found"
    end
end

@testset "Sample bounds" begin
    @testset "Get sampling bounds" begin
        dom = deepcopy(ADRIA_DOM_45)
        cb_calib_groups::Vector{Int64} = dom.loc_data.CB_CALIB_GROUPS
        ms = ADRIA._filtered_model_spec(ADRIA.model_spec(dom), cb_calib_groups)

        @testset "Continuous variables" begin
            continuous_factors = ms[(ms.ptype .∉ [ADRIA.DISCRETE_FACTOR_TYPES]), :]
            sample_factors = continuous_factors[rand(1:nrow(continuous_factors), 5), :]

            for factor in eachrow(sample_factors)
                fn = factor.fieldname
                @test ADRIA.get_bounds(dom, fn) == factor.dist_params[1:2] ||
                    "get_bounds mismatch for $fn"
                @test ADRIA.get_attr(dom, fn, :default_dist_params) ==
                      factor.default_dist_params ||
                    "get_attr(:default_dist_params) mismatch for $fn"
            end
        end

        @testset "Discrete variables" begin
            discrete_factors = ms[(ms.ptype .∈ [ADRIA.DISCRETE_FACTOR_TYPES]), :]
            sample_factors = discrete_factors[rand(1:nrow(discrete_factors), 5), :]

            for factor in eachrow(sample_factors)
                fn = factor.fieldname
                @test ADRIA.get_bounds(dom, fn)[1] == factor.dist_params[1] ||
                    "get_bounds lower bound mismatch for $fn"
                @test ADRIA.get_bounds(dom, fn)[2] == factor.dist_params[2] ||
                    "get_bounds upper bound mismatch for $fn"
                @test ADRIA.get_attr(dom, fn, :default_dist_params)[1] ==
                      factor.default_dist_params[1] ||
                    "get_attr(:default_dist_params) lower bound mismatch for $fn"
                @test ADRIA.get_attr(dom, fn, :default_dist_params)[2] ==
                      factor.default_dist_params[2] ||
                    "get_attr(:default_dist_params) upper bound mismatch for $fn"
            end
        end
    end

    @testset "Set new sampling bounds" begin
        set_factor_bounds! = ADRIA.set_factor_bounds!

        dom = deepcopy(ADRIA_DOM_45)
        num_samples = 32
        cb_calib_groups::Vector{Int64} = dom.loc_data.CB_CALIB_GROUPS
        ms = ADRIA._filtered_model_spec(ADRIA.model_spec(dom), cb_calib_groups)

        test_components = ["EnvironmentalLayer", "Coral"]

        function _test_bounds(
            scens::DataFrame, factor_mask::BitVector, bounds_ranges::Vector
        )
            # scens[:,collect(factor_fieldnames)] == scens[:,factor_mask]

            filt_scens = Matrix(scens[scens.guided .> 0, factor_mask])
            min_scens, max_scens = vcat.([extrema(x) for x in eachcol(filt_scens)]...)
            min_bounds, max_bounds = vcat.(extrema.(collect.(bounds_ranges))...)

            err_msg = "Sampled continuous factor is outside of specified new bounds."
            # @test all((max_scens .<= max_bounds) .&& (min_scens .>= min_bounds)) || err_msg
            return all((max_scens .<= max_bounds) .&& (min_scens .>= min_bounds)) || err_msg
        end

        @testset "Uniform distributions" begin
            factor_mask = (ms.component .∈ [test_components]) .&& (ms.dist .== Uniform)
            factors = ms[factor_mask, :]
            factor_fieldnames = (factors.fieldname...,)

            @testset "set_factor_bounds!" begin
                bounds_ranges = [range(b[1], b[2], 5) for b in factors.default_dist_params]
                new_bounds = Tuple.(sort.(rand.(bounds_ranges, 2)))
                dom = set_factor_bounds!(dom; NamedTuple{factor_fieldnames}(new_bounds)...)
                scens = ADRIA.sample_guided(dom, num_samples)

                @test _test_bounds(scens, factor_mask, bounds_ranges)
            end

            @testset "set to default bounds" begin
                new_bounds = ADRIA.get_attr.(
                    [dom], factor_fieldnames, [:default_dist_params]
                )
                dom = ADRIA.set_factor_bounds!(
                    dom; NamedTuple{factor_fieldnames}(new_bounds)...
                )

                factor_params = ms[ms.fieldname .∈ [factor_fieldnames], :]
                @test all(factor_params.dist_params .== factor_params.default_dist_params)
                scens = ADRIA.sample(dom, num_samples)

                @test _test_bounds(scens, factor_mask, new_bounds)
            end
        end

        # @testset "DiscreteUniform distributions" begin
        #     factor_mask = ms.component .∈ [test_components] .&& ms.dist .== DiscreteUniform
        #     factors = ms[factor_mask, :]
        #     factor_fieldnames = (factors.fieldname...,)

        #     @testset "set_factor_bounds!" begin
        #         bounds_ranges = [b[1]:b[2] for b in factors.default_dist_params]
        #         new_bounds = Tuple.(sort.(rand.(bounds_ranges, 2)))
        #         dom = set_factor_bounds!(dom; NamedTuple{factor_fieldnames}(new_bounds)...)
        #         scens = ADRIA.sample_guided(dom, num_samples)

        #         @test _test_bounds(scens, factor_mask, bounds_ranges)
        #     end

        #     @testset "get_default_dist_params" begin
        #         new_bounds =
        #             ADRIA.get_attr.([dom], factor_fieldnames, [:default_dist_params])
        #         dom = set_factor_bounds!(dom; NamedTuple{factor_fieldnames}(new_bounds)...)

        #         factor_params = ms[ms.fieldname .∈ [factor_fieldnames], :]
        #         @test all(factor_params.dist_params .== factor_params.default_dist_params)
        #         scens = ADRIA.sample(dom, num_samples)

        #         @test _test_bounds(scens, factor_mask, new_bounds)
        #     end
        # end

        # @testset "DiscreteOrderedUniformDist distributions" begin
        #     factor_mask =
        #         ms.component .∈ [test_components] .&&
        #         ms.dist .== ADRIA.DiscreteOrderedUniformDist .&&
        #         .!contains.(String.(ms.fieldname), "reactive")
        #     factors = ms[factor_mask, :]
        #     factor_fieldnames = (factors.fieldname...,)
        #     @testset "set_factor_bounds!" begin
        #         bounds_ranges = [b[1]:b[3]:b[2] for b in factors.default_dist_params]
        #         new_bounds = Tuple.(sort.(rand.(bounds_ranges, 2)))

        #         new_dist_params = [
        #             (b[1], b[2], old_dist_params[3]) for
        #             (b, old_dist_params) in zip(new_bounds, factors.default_dist_params)
        #         ]
        #         dom = ADRIA.set_factor_bounds!(
        #             dom; NamedTuple{factor_fieldnames}(new_dist_params)...
        #         )
        #         scens = ADRIA.sample_guided(dom, num_samples)

        #         @test _test_bounds(scens, factor_mask, bounds_ranges)
        #     end

        #     @testset "get_default_dist_params" begin
        #         new_bounds =
        #             ADRIA.get_attr.([dom], factor_fieldnames, [:default_dist_params])
        #         dom = set_factor_bounds!(dom; NamedTuple{factor_fieldnames}(new_bounds)...)

        #         factor_params = ms[ms.fieldname .∈ [factor_fieldnames], :]
        #         @test all(factor_params.dist_params .== factor_params.default_dist_params)
        #         scens = ADRIA.sample(dom, num_samples)

        #         @test _test_bounds(scens, factor_mask, new_bounds)
        #     end
        # end

        @testset "TriangularDist distributions" begin
            factor_mask = ms.component .∈ [test_components] .&& ms.dist .== TriangularDist
            factors = ms[factor_mask, :]
            factor_fieldnames = (factors.fieldname...,)
            @testset "set_factor_bounds!" begin
                bounds_ranges = [range(b[1], b[2], 5) for b in factors.default_dist_params]
                new_bounds = Tuple.(sort.(rand.(bounds_ranges, 2)))
                mode_ranges = [range(nb[1], nb[2], 5) for nb in new_bounds]
                new_modes = (rand.(mode_ranges))
                new_dist_params = [(b[1], b[2], p) for (b, p) in zip(new_bounds, new_modes)]
                dom = set_factor_bounds!(
                    dom; NamedTuple{factor_fieldnames}(new_dist_params)...
                )
                scens = ADRIA.sample_guided(dom, num_samples)

                @test _test_bounds(scens, factor_mask, new_bounds)
            end

            @testset "get_default_dist_params" begin
                new_bounds = ADRIA.get_attr.(
                    [dom], factor_fieldnames, [:default_dist_params]
                )
                dom = set_factor_bounds!(dom; NamedTuple{factor_fieldnames}(new_bounds)...)

                factor_params = ms[ms.fieldname .∈ [factor_fieldnames], :]
                @test all(factor_params.dist_params .== factor_params.default_dist_params)
                scens = ADRIA.sample(dom, num_samples)

                @test _test_bounds(scens, factor_mask, new_bounds)
            end
        end

        @testset "DiscreteTriangularDist distributions" begin
            factor_mask =
                (ms.component .∈ [test_components]) .&&
                (ms.dist .== ADRIA.DiscreteTriangularDist) .&&
                .!contains.(String.(ms.fieldname), "reactive")
            factors = ms[factor_mask, :]

            if nrow(factors) == 0
                @test_skip "No non-reactive DiscreteTriangularDist factors to test"
            else
                factor_fieldnames = (factors.fieldname...,)
                @testset "set_factor_bounds!" begin
                    bounds_ranges = [b[1]:b[2] for b in factors.default_dist_params]
                    new_bounds = Tuple.(sort.(rand.(bounds_ranges, 2)))
                    new_mode_ranges = [nb[1]:nb[2] for nb in new_bounds]
                    new_peaks = (rand.(new_mode_ranges))
                    new_dist_params = [
                        (b[1], b[2], p) for (b, p) in zip(new_bounds, new_peaks)
                    ]
                    dom = set_factor_bounds!(
                        dom; NamedTuple{factor_fieldnames}(new_dist_params)...
                    )
                    scens = ADRIA.sample_guided(dom, num_samples)

                    @test _test_bounds(scens, factor_mask, new_bounds)
                end

                @testset "get_default_dist_params" begin
                    new_bounds = ADRIA.get_attr.(
                        [dom], factor_fieldnames, [:default_dist_params]
                    )
                    dom = set_factor_bounds!(
                        dom; NamedTuple{factor_fieldnames}(new_bounds)...
                    )

                    factor_params = ms[ms.fieldname .∈ [factor_fieldnames], :]
                    @test all(
                        factor_params.dist_params .== factor_params.default_dist_params
                    )
                    scens = ADRIA.sample(dom, num_samples)

                    @test _test_bounds(scens, factor_mask, new_bounds)
                end
            end
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Tests for component_dependency_dag_v3 implementation
# (sampling_dependencies.jl + sampling_transforms.jl)
# Ref: .claude/plans/components/component_dependency_dag_v3.md §4 Steps 5,6,8,9
# ─────────────────────────────────────────────────────────────────────────────

@testset "Dependency DAG — _validate_dependencies" begin
    dom = deepcopy(ADRIA_DOM_45)
    spec = ADRIA.model_spec(dom)
    @test_nowarn ADRIA._validate_dependencies(spec)
end

@testset "Dependency DAG — _resolve_conditional_spec!" begin
    PERIODIC = ADRIA.DECISION_STRATEGY[:periodic]

    @testset "CF regime (guided=-1)" begin
        dom = deepcopy(ADRIA_DOM_45)
        spec = ADRIA.model_spec(dom)
        ADRIA._resolve_conditional_spec!(spec, (guided=-1.0,))

        for col in ADRIA._GROUP_MEMBERS[:intervention_group]
            row = spec[spec.fieldname .== col, :]
            isempty(row) && continue
            @test row.is_constant[1] && isequal(row.val[1], 0.0) ||
                "CF: :$col should be fixed to 0.0, got val=$(row.val[1])"
        end

        # PINNED REGRESSION (§2.3 changes 1 & 2): strategy cols must be -1.0 for CF.
        # Previously had negate-sign inversion AND wrong fix_to (0.0 not -1.0).
        for col in ADRIA._GROUP_MEMBERS[:strategy_group]
            row = spec[spec.fieldname .== col, :]
            isempty(row) && continue
            @test row.is_constant[1] && isequal(row.val[1], -1.0) ||
                "CF: :$col should be -1.0 sentinel (§2.3 change 2). " *
                  "Got val=$(row.val[1]). Check negate-sign or fix_to."
        end

        for col in ADRIA._GROUP_MEMBERS[:criteria_weights]
            row = spec[spec.fieldname .== col, :]
            isempty(row) && continue
            @test row.is_constant[1] && isequal(row.val[1], 0.0) ||
                "CF: criteria weight :$col should be 0.0"
        end

        for col in ADRIA._GROUP_MEMBERS[:depth_thresholds]
            row = spec[spec.fieldname .== col, :]
            isempty(row) && continue
            @test row.is_constant[1] && isequal(row.val[1], 0.0) ||
                "CF: depth threshold :$col should be 0.0"
        end
    end

    @testset "Unguided regime (guided=0)" begin
        dom = deepcopy(ADRIA_DOM_45)
        spec = ADRIA.model_spec(dom)
        ADRIA._resolve_conditional_spec!(spec, (guided=0.0,))

        ph_row = spec[spec.fieldname .== :plan_horizon, :]
        if !isempty(ph_row)
            @test ph_row.is_constant[1] && isequal(ph_row.val[1], 0.0) ||
                "Unguided: :plan_horizon should be 0.0"
        end

        for col in ADRIA._GROUP_MEMBERS[:criteria_weights]
            row = spec[spec.fieldname .== col, :]
            isempty(row) && continue
            @test row.is_constant[1] && isequal(row.val[1], 0.0) ||
                "Unguided: criteria weight :$col should be 0.0"
        end

        # PINNED REGRESSION (§2.3 change 1): strategy cols must be PERIODIC for unguided.
        for col in ADRIA._GROUP_MEMBERS[:strategy_group]
            row = spec[spec.fieldname .== col, :]
            isempty(row) && continue
            @test row.is_constant[1] && isequal(row.val[1], Float64(PERIODIC)) ||
                "Unguided: :$col should be PERIODIC ($PERIODIC). " *
                  "Got val=$(row.val[1]). Check negate-sign (§2.3 change 1)."
        end
    end

    @testset "Guided regime (guided=1)" begin
        dom = deepcopy(ADRIA_DOM_45)
        spec = ADRIA.model_spec(dom)
        already_constant = Set(spec[spec.is_constant, :fieldname])
        ADRIA._resolve_conditional_spec!(spec, (guided=1.0,))

        newly_fixed = setdiff(Set(spec[spec.is_constant, :fieldname]), already_constant)
        @test isempty(newly_fixed) ||
            "Guided: unexpected columns newly fixed: $newly_fixed"
    end
end

@testset "Dependency DAG — _transform_group_columns snapshot (E4)" begin
    # Pinned snapshot: strategy cols must NOT appear in prefix-based groups.
    # Their inclusion was the d5871840 regression (overwrites -1.0 CF sentinel).
    dom = deepcopy(ADRIA_DOM_45)
    scens = ADRIA.sample(dom, 8)

    seed_cols = ADRIA._transform_group_columns(scens, :seed_group)
    fog_cols = ADRIA._transform_group_columns(scens, :fog_group)
    mc_cols = ADRIA._transform_group_columns(scens, :mc_group)

    @test :seed_strategy ∉ seed_cols ||
        ":seed_strategy in :seed_group reproduces d5871840 regression"
    @test :fog_strategy ∉ fog_cols ||
        ":fog_strategy in :fog_group reproduces d5871840 regression"
    @test :mc_strategy ∉ mc_cols ||
        ":mc_strategy in :mc_group reproduces d5871840 regression"

    @test !isempty(seed_cols) || ":seed_group resolved empty — prefix filter broken"
    @test !isempty(fog_cols) || ":fog_group resolved empty — prefix filter broken"
    @test !isempty(mc_cols) || ":mc_group resolved empty — prefix filter broken"
end

@testset "Dependency DAG — fog_strategy/mc_strategy gated on intervention activity (Issue #1132)" begin
    # fog_strategy/mc_strategy used to be free to draw reactive even when their
    # own intervention (fogging/N_mc_settlers) was fixed inactive across all
    # rows, keeping reactive_group alive as sampling noise. They should now be
    # fixed to PERIODIC whenever their intervention is inactive.
    PERIODIC = ADRIA.DECISION_STRATEGY[:periodic]

    dom = deepcopy(ADRIA_DOM_45)
    ADRIA.fix_factor!(dom; fogging=0.0, N_mc_settlers=0.0)
    scens = ADRIA.sample(dom, 64)

    @testset "fog_strategy fixed when fogging inactive" begin
        if :fog_strategy in propertynames(scens) && :fogging in propertynames(scens)
            not_cf = scens.fog_strategy .!= -1.0
            @test all(scens[not_cf, :fog_strategy] .== Float64(PERIODIC)) ||
                "fog_strategy should be PERIODIC ($PERIODIC) for all non-CF rows once " *
                  "fogging is fixed inactive"
        end
    end

    @testset "mc_strategy fixed when N_mc_settlers inactive" begin
        if :mc_strategy in propertynames(scens) && :N_mc_settlers in propertynames(scens)
            not_cf = scens.mc_strategy .!= -1.0
            @test all(scens[not_cf, :mc_strategy] .== Float64(PERIODIC)) ||
                "mc_strategy should be PERIODIC ($PERIODIC) for all non-CF rows once " *
                  "N_mc_settlers is fixed inactive"
        end
    end

    @testset "reactive_group dropped when seed is the only reactive-capable lever" begin
        if all(
            c -> c in propertynames(scens),
            [:seed_strategy, :fog_strategy, :mc_strategy, :reactive_response_delay]
        )
            periodic_seed = scens.seed_strategy .== Float64(PERIODIC)
            @test all(scens[periodic_seed, :reactive_response_delay] .== 0.0) ||
                "reactive_group should be 0.0 when seed_strategy is periodic and " *
                  "fog/mc are fixed inactive (no reactive-capable lever remains)"
        end
    end
end

@testset "Dependency DAG — seed_wave_stress ordering (§2.12 pt 1)" begin
    # seed_wave_stress is zeroed AFTER mcda_normalize. For wave_scenario==0
    # guided+seeded rows, remaining 7 seed weights must sum to <1 (not renormalised).
    dom = deepcopy(ADRIA_DOM_45)
    num_samples = 128
    scens = ADRIA.sample_guided(dom, num_samples)
    ms = ADRIA.model_spec(dom)
    seed_weights = ADRIA.component_params(ms, ADRIA.SeedCriteriaWeights).fieldname

    not_seeded = ADRIA.no_seeding(scens)
    guided_seeded_mask = .!not_seeded .& (scens.guided .> 0)
    wave_zero_mask = (scens.wave_scenario .== 0.0)
    target_mask = guided_seeded_mask .& wave_zero_mask

    if any(target_mask)
        @test all(scens[target_mask, :seed_wave_stress] .== 0.0) ||
            "seed_wave_stress should be 0 for wave_scenario==0 guided+seeded rows"

        other_seed_weights = filter(w -> w != :seed_wave_stress, seed_weights)
        row_sums = vec(sum(Matrix(scens[target_mask, other_seed_weights]); dims=2))
        @test all(row_sums .< 1.0) ||
            "Remaining seed weights sum to 1 — ordering wrong: mcda_normalize ran " *
              "AFTER zeroing seed_wave_stress (§2.12 pt 1)"
    else
        @info "No wave_scenario==0 guided+seeded rows in $num_samples samples; " *
            "seed_wave_stress ordering test inconclusive"
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Tests for paired counterfactual sampling (sampling_design_plan_v4.md §5 Step 3)
# ─────────────────────────────────────────────────────────────────────────────

@testset "derive_cf" begin
    dom = deepcopy(ADRIA_DOM_45)
    num_samples = 32
    spec = ADRIA.model_spec(dom)
    intervention_samples = ADRIA.sample_guided(dom, num_samples)
    cf_samples = ADRIA.derive_cf(intervention_samples, spec)

    zeroed_cols = [
        :N_seed_TA, :N_seed_CA, :N_seed_CNA, :N_seed_SM, :N_seed_LM,
        :fogging, :SRM, :N_mc_settlers
    ]
    zeroed_cols = filter(c -> c in propertynames(cf_samples), zeroed_cols)
    strategy_cols = filter(
        c -> c in propertynames(cf_samples), [:seed_strategy, :fog_strategy, :mc_strategy]
    )

    @testset "CF zeroing" begin
        @test all(all(cf_samples[:, c] .== 0.0) for c in zeroed_cols) ||
            "Expected all of $zeroed_cols to be 0.0 in CF rows"
    end

    @testset "strategy sentinel" begin
        @test all(all(cf_samples[:, c] .== -1.0) for c in strategy_cols) ||
            "Expected all of $strategy_cols to be -1.0 (CF sentinel) in CF rows"
    end

    @testset "env params unchanged" begin
        # Columns derive_cf may fix (per _PARAM_DEPENDENCIES/_GROUP_MEMBERS) are not
        # "non-intervention" even though some (criteria weights, depth thresholds)
        # fall outside the Intervention component.
        cf_affected_cols = Set(
            vcat(
                [:guided],
                reduce(
                    vcat,
                    [
                        get(ADRIA._GROUP_MEMBERS, dep.child, [dep.child])
                        for dep in ADRIA._PARAM_DEPENDENCIES
                    ]
                )
            )
        )
        env_cols = [
            c for c in propertynames(intervention_samples) if c ∉ cf_affected_cols
        ]
        @test intervention_samples[:, env_cols] == cf_samples[:, env_cols] ||
            "Non-intervention parameter columns differ between paired rows"
    end

    @testset "parent params zeroed" begin
        nonzero_mask = vec(
            any(Matrix(intervention_samples[:, zeroed_cols]) .> 0; dims=2)
        )
        @test any(nonzero_mask) ||
            "No sampled guided rows had non-zero intervention parent values to test against"
        @test all(all(cf_samples[nonzero_mask, c] .== 0.0) for c in zeroed_cols) ||
            "Expected $zeroed_cols to be 0.0 in CF rows paired to non-zero intervention rows"
    end
end

@testset "sample_paired" begin
    dom = deepcopy(ADRIA_DOM_45)
    n = 16
    scens = ADRIA.sample_paired(dom, n)

    @testset "shape" begin
        @test nrow(scens) == 2n || "Expected 2n rows, got $(nrow(scens))"
        @test count(==("intervention"), scens.run_type) == n
        @test count(==("counterfactual"), scens.run_type) == n
    end

    @testset "pairing" begin
        cf_affected_cols = Set(
            vcat(
                [:guided],
                reduce(
                    vcat,
                    [
                        get(ADRIA._GROUP_MEMBERS, dep.child, [dep.child])
                        for dep in ADRIA._PARAM_DEPENDENCIES
                    ]
                )
            )
        )
        env_cols = [
            c for c in propertynames(scens) if c ∉ cf_affected_cols && c != :run_type
        ]

        interv_rows = scens[scens.run_type .== "intervention", env_cols]
        cf_rows = scens[scens.run_type .== "counterfactual", env_cols]
        @test interv_rows == cf_rows ||
            "Non-intervention parameters do not match between paired rows"
    end

    @testset "run_type type" begin
        @test eltype(scens.run_type) == String
    end
end

@testset "sample_paired_stratified" begin
    dom = deepcopy(ADRIA_DOM_45)
    n_per_stratum = 8
    # Test_domain only ships RCP45 DHW/wave data; repeat the same RCP as two
    # independently-scrambled strata to validate shape/labelling without
    # requiring additional RCP datasets.
    ssp_strata = ["45", "45"]
    scens = ADRIA.sample_paired_stratified(dom, n_per_stratum, ssp_strata)

    @testset "shape" begin
        @test nrow(scens) == 2 * n_per_stratum * length(ssp_strata) ||
            "Expected $(2 * n_per_stratum * length(ssp_strata)) rows, got $(nrow(scens))"
        @test count(==("45"), scens.ssp) == 2 * n_per_stratum * length(ssp_strata)
    end
end

@testset "derive_unguided" begin
    dom = deepcopy(ADRIA_DOM_45)
    num_samples = 32
    spec = ADRIA.model_spec(dom)
    intervention_samples = ADRIA.sample_guided(dom, num_samples)
    ug_samples = ADRIA.derive_unguided(intervention_samples, spec)

    # Unlike CF, unguided rows keep intervention deployment amounts and depth
    # thresholds unchanged (:intervention_group / :depth_thresholds stay active
    # at guided=0.0) — only criteria weights / plan_horizon are zeroed.
    zeroed_cols = filter(
        c -> c in propertynames(ug_samples),
        [
            :seed_heat_stress, :seed_wave_stress, :seed_in_connectivity,
            :fog_heat_stress, :mc_heat_stress, :plan_horizon
        ]
    )
    unchanged_cols = filter(
        c -> c in propertynames(ug_samples),
        [:N_seed_TA, :N_seed_CA, :fogging, :SRM, :N_mc_settlers, :depth_min]
    )
    strategy_cols = filter(
        c -> c in propertynames(ug_samples), [:seed_strategy, :fog_strategy, :mc_strategy]
    )

    @testset "criteria weights/plan_horizon zeroed" begin
        @test all(all(ug_samples[:, c] .== 0.0) for c in zeroed_cols) ||
            "Expected all of $zeroed_cols to be 0.0 in unguided rows"
    end

    @testset "strategy fixed to periodic" begin
        periodic = ADRIA.decision.DECISION_STRATEGY[:periodic]
        @test all(all(ug_samples[:, c] .== periodic) for c in strategy_cols) ||
            "Expected all of $strategy_cols to be $periodic (periodic) in unguided rows"
    end

    @testset "intervention deployment amounts unchanged" begin
        @test intervention_samples[:, unchanged_cols] == ug_samples[:, unchanged_cols] ||
            "Expected intervention deployment amounts to be unchanged in unguided rows"
    end

    @testset "guided column set to 0.0" begin
        @test all(ug_samples.guided .== 0.0)
    end
end

@testset "sample_matched" begin
    dom = deepcopy(ADRIA_DOM_45)
    n = 16

    @testset "default regimes (3-way)" begin
        scens = ADRIA.sample_matched(dom, n)

        @test nrow(scens) == 3n || "Expected 3n rows, got $(nrow(scens))"
        @test count(==("intervention"), scens.run_type) == n
        @test count(==("counterfactual"), scens.run_type) == n
        @test count(==("unguided"), scens.run_type) == n
        @test eltype(scens.run_type) == String

        cf_affected_cols = Set(
            vcat(
                [:guided],
                reduce(
                    vcat,
                    [
                        get(ADRIA._GROUP_MEMBERS, dep.child, [dep.child])
                        for dep in ADRIA._PARAM_DEPENDENCIES
                    ]
                )
            )
        )
        env_cols = [
            c for c in propertynames(scens) if c ∉ cf_affected_cols && c != :run_type
        ]

        interv_rows = scens[scens.run_type .== "intervention", env_cols]
        cf_rows = scens[scens.run_type .== "counterfactual", env_cols]
        ug_rows = scens[scens.run_type .== "unguided", env_cols]
        @test interv_rows == cf_rows ||
            "Non-intervention parameters do not match between guided/CF rows"
        @test interv_rows == ug_rows ||
            "Non-intervention parameters do not match between guided/unguided rows"
    end

    @testset "regimes subsetting" begin
        ug_only = ADRIA.sample_matched(dom, n; regimes=(:unguided,))
        @test nrow(ug_only) == 2n
        @test count(==("counterfactual"), ug_only.run_type) == 0
        @test count(==("unguided"), ug_only.run_type) == n

        cf_only = ADRIA.sample_matched(dom, n; regimes=(:counterfactual,))
        @test nrow(cf_only) == 2n
        @test count(==("unguided"), cf_only.run_type) == 0
        @test count(==("counterfactual"), cf_only.run_type) == n
    end

    @testset "invalid regimes" begin
        @test_throws ArgumentError ADRIA.sample_matched(dom, n; regimes=())
        @test_throws ArgumentError ADRIA.sample_matched(dom, n; regimes=(:bogus,))
    end
end

@testset "sample_matched_stratified" begin
    dom = deepcopy(ADRIA_DOM_45)
    n_per_stratum = 8
    ssp_strata = ["45", "45"]
    scens = ADRIA.sample_matched_stratified(dom, n_per_stratum, ssp_strata)

    @testset "shape" begin
        @test nrow(scens) == 3 * n_per_stratum * length(ssp_strata) ||
            "Expected $(3 * n_per_stratum * length(ssp_strata)) rows, got $(nrow(scens))"
        @test count(==("45"), scens.ssp) == 3 * n_per_stratum * length(ssp_strata)
    end
end
