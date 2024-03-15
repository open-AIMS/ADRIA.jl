using ADRIA
using ADRIA.Distributions

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

@testset "sample" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
    num_samples = 32
    scens = ADRIA.sample(dom, num_samples)
    ms = ADRIA.model_spec(dom)
    constant_params = ms.is_constant

    @testset "constant params are constant" begin
        @test all(values(scens[1, constant_params]) .== values(scens[end, constant_params])) ||
              "Constant params are not constant!"
    end

    @testset "values are within expected bounds" begin
        #interv = (ms.component .== "Intervention") .& .!(constant_params)
        lb = values(ms[:, :lower_bound])
        ub = values(ms[:, :upper_bound])

        not_cw_mask = ms.component .∉ [("SeedCriteriaWeights", "FogCriteriaWeights")]
        not_cw_lb, not_cw_ub = lb[not_cw_mask], ub[not_cw_mask]

        eco = (ms.component .== "Coral") .& .!(constant_params)

        msg = "Sampled values were not in expected bounds!"
        coral_msg = "Sampled coral values were not in expected bounds!"
        for i in 1:num_samples
            # Filter CriteriaWeights factors
            scen_vals = values(scens[i, :])
            not_cw_scen_vals = scen_vals[not_cw_mask]

            if scens[i, :guided] > 0
                cond = not_cw_lb .<= not_cw_scen_vals .<= not_cw_ub
                @test all(cond) || "$msg | $(ms[.!(cond), :]) | $(not_cw_scen_vals[.!(cond)])"
                continue
            end
            # When no interventions are used, e.g., for counterfactual or unguided scenarios
            # (guided ∈ [-1, 0]) intervention parameters are set to 0 so only check ecological values
            cond = lb[eco] .<= scen_vals[eco] .<= ub[eco]
            @test all(cond) || "$coral_msg | $(ms[.!(cond), :]) | $(scen_vals[eco][.!(cond)])"

            # Note: Test to ensure all intervention factors are set to 0 is covered by the guided
            # sampling test below
        end
    end
end

@testset "Targeted sampling" begin
    @testset "Counterfactual sampling" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
        num_samples = 32
        scens = ADRIA.sample_cf(dom, num_samples)

        @test all(scens.guided .== -1) || "Intervention scenarios found"

        # Get Intervention params
        interv_params =
            string.(
                ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Intervention).fieldname
            )

        # Ensure all interventions are deactivated (ignoring the "guided" factor)
        interv_params = String[ip for ip in interv_params if ip != "guided"]
        @test all(all.(==(0), eachcol(scens[:, interv_params]))) ||
              "Intervention factors with values > 0 found"
    end

    @testset "Guided sampling" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
        num_samples = 32
        scens = ADRIA.sample_guided(dom, num_samples)
        ms = ADRIA.model_spec(dom)

        @test all(scens.guided .> 0) || "Non-intervention scenarios found"

        # Get Intervention params
        interv_params =
            string.(
                ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Intervention).fieldname
            )

        # Ignore guided
        interv_params = String[ip for ip in interv_params if ip != "guided"]

        # Ensure at least one intervention is active
        @test all(any.(>(0), eachcol(scens[:, interv_params]))) ||
              "All intervention factors had values <= 0"

        seed_weights = ADRIA.component_params(ms, ADRIA.SeedCriteriaWeights).fieldname
        fog_weights = ADRIA.component_params(ms, ADRIA.FogCriteriaWeights).fieldname

        @test all(abs.(sum(Matrix(scens[:, seed_weights]); dims=2) .- 1.0) .< 10e-6) ||
              "Some seeding weights are not properly normalized."
        @test all(abs.(sum(Matrix(scens[:, fog_weights]); dims=2) .- 1.0) .< 10e-6) ||
              "Some fogging weights are not properly normalized."
    end

    @testset "Unguided sampling" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
        num_samples = 32
        scens = ADRIA.sample_unguided(dom, num_samples)

        @test all(scens.guided .== 0) || "Intervention or counterfactual scenarios found"

        # Get Intervention params
        interv_params =
            string.(
                ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Intervention).fieldname
            )

        # Ignore guided and planning horizon
        interv_params = String[
            ip for ip in interv_params if ip != "plan_horizon" && ip != "guided"
        ]

        # Ensure at least one intervention is active
        @test all(any.(>(0), eachcol(scens[:, interv_params]))) ||
              "All intervention factors had values <= 0"
    end

    @testset "Site selection sampling" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
        num_samples = 32
        scens = ADRIA.sample_selection(dom, num_samples)

        @test all(scens.guided .> 0) || "Intervention or counterfactual scenarios found"

        # Get Intervention params
        ms = ADRIA.model_spec(dom)
        target_params =
            string.(
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
        @test all(any.(>(0), eachcol(scens[:, target_params]))) ||
              "All target factors had values <= 0"

        # Check that all coral parameters are set to their nominated default values
        coral_params = ADRIA.component_params(ms, ADRIA.Coral).fieldname

        @test all([
            all(scens[:, c] .== ms[ms.fieldname.==c, :val][1]) for c in coral_params
        ]) || "Non-default coral parameter value found"
    end
end

@testset "Sample bounds" begin
    @testset "Get sampling bounds" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
        ms = ADRIA.model_spec(dom)

        @testset "Continuous variables" begin
            continuous_factors = ms[(ms.ptype.∉[ADRIA.DISCRETE_FACTOR_TYPES]), :]

            for factor in eachrow(continuous_factors)
                fn = factor.fieldname
                @test ADRIA.get_bounds(dom, fn) == factor.dist_params[1:2]
                @test ADRIA.get_attr(dom, fn, :default_dist_params) ==
                      factor.default_dist_params
            end
        end

        @testset "Discrete variables" begin
            discrete_factors = ms[(ms.ptype.∈[ADRIA.DISCRETE_FACTOR_TYPES]), :]
            for factor in eachrow(discrete_factors)
                fn = factor.fieldname
                @test ADRIA.get_bounds(dom, fn)[1] == factor.dist_params[1]
                @test ADRIA.get_bounds(dom, fn)[2] == factor.dist_params[2]
                @test ADRIA.get_attr(dom, fn, :default_dist_params)[1] ==
                      factor.default_dist_params[1]
                @test ADRIA.get_attr(dom, fn, :default_dist_params)[2] ==
                      factor.default_dist_params[2]
            end
        end
    end

    @testset "Set new sampling bounds" begin
        set_factor_bounds = ADRIA.set_factor_bounds

        dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
        num_samples = 32
        ms = ADRIA.model_spec(dom)

        test_components = ["EnvironmentalLayer", "Intervention", "Coral"]

        function _test_bounds(scens::DataFrame, factor_mask::BitVector, bounds_ranges::Vector)
            filt_scens = Matrix(scens[scens.guided.>0, factor_mask])
            min_scens, max_scens = vcat.([extrema(x) for x in eachcol(filt_scens)]...)
            min_bounds, max_bounds = vcat.(extrema.(collect.(bounds_ranges))...)

            err_msg = "Sampled continuous factor is outside of specified new bounds."
            @test all((max_scens .<= max_bounds) .&& (min_scens .>= min_bounds)) || err_msg
        end

        @testset "Uniform distributions" begin
            factor_mask = (ms.component .∈ [test_components]) .&& (ms.dist .== Uniform)
            factors = ms[factor_mask, :]
            factor_fieldnames = (factors.fieldname...,)

            @testset "set_factor_bounds" begin
                bounds_ranges = [range(b[1], b[2], 5) for b in factors.default_dist_params]
                new_bounds = Tuple.(sort.(rand.(bounds_ranges, 2)))
                dom = set_factor_bounds(dom; NamedTuple{factor_fieldnames}(new_bounds)...)
                scens = ADRIA.sample_guided(dom, num_samples)

                _test_bounds(scens, factor_mask, bounds_ranges)
            end

            @testset "set to default bounds" begin
                new_bounds = ADRIA.get_attr.([dom], factor_fieldnames, [:default_dist_params])
                dom = set_factor_bounds(dom; NamedTuple{factor_fieldnames}(new_bounds)...)

                factor_params = dom.model[ms.fieldname.∈[factor_fieldnames]][1]
                @test all(factor_params.dist_params .== factor_params.default_dist_params)

                scens = ADRIA.sample(dom, num_samples)
                _test_bounds(scens, factor_mask, new_bounds)
            end
        end

        @testset "DiscreteUniform distributions" begin
            factor_mask = ms.component .∈ [test_components] .&& ms.dist .== DiscreteUniform
            factors = ms[factor_mask, :]
            factor_fieldnames = (factors.fieldname...,)

            @testset "set_factor_bounds" begin
                bounds_ranges = [b[1]:b[2] for b in factors.default_dist_params]
                new_bounds = Tuple.(sort.(rand.(bounds_ranges, 2)))
                dom = set_factor_bounds(dom; NamedTuple{factor_fieldnames}(new_bounds)...)
                scens = ADRIA.sample_guided(dom, num_samples)

                _test_bounds(scens, factor_mask, bounds_ranges)
            end

            @testset "get_default_dist_params" begin
                new_bounds = ADRIA.get_attr.([dom], factor_fieldnames, [:default_dist_params])
                dom = set_factor_bounds(dom; NamedTuple{factor_fieldnames}(new_bounds)...)

                factor_params = dom.model[ms.fieldname.∈[factor_fieldnames]][1]
                @test all(factor_params.dist_params .== factor_params.default_dist_params)

                scens = ADRIA.sample(dom, num_samples)

                _test_bounds(scens, factor_mask, new_bounds)
            end
        end

        @testset "DiscreteOrderedUniformDist distributions" begin
            factor_mask = ms.component .∈ [test_components] .&& ms.dist .== ADRIA.DiscreteOrderedUniformDist
            factors = ms[factor_mask, :]
            factor_fieldnames = (factors.fieldname...,)
            @testset "set_factor_bounds" begin
                bounds_ranges = [b[1]:b[2] for b in factors.default_dist_params]
                new_bounds = Tuple.(sort.(rand.(bounds_ranges, 2)))
                new_steps = [round((nb[2] - nb[1]) / 10) for nb in new_bounds]
                new_dist_params = [(b[1], b[2], s) for (b, s) in zip(new_bounds, new_steps)]
                dom = set_factor_bounds(dom; NamedTuple{factor_fieldnames}(new_dist_params)...)
                scens = ADRIA.sample_guided(dom, num_samples)

                _test_bounds(scens, factor_mask, bounds_ranges)
            end

            @testset "get_default_dist_params" begin
                new_bounds = ADRIA.get_attr.([dom], factor_fieldnames, [:default_dist_params])
                dom = set_factor_bounds(dom; NamedTuple{factor_fieldnames}(new_bounds)...)

                factor_params = dom.model[ms.fieldname.∈[factor_fieldnames]][1]
                @test all(factor_params.dist_params .== factor_params.default_dist_params)

                scens = ADRIA.sample(dom, num_samples)
                _test_bounds(scens, factor_mask, new_bounds)
            end
        end

        @testset "TriangularDist distributions" begin
            factor_mask = ms.component .∈ [test_components] .&& ms.dist .== TriangularDist
            factors = ms[factor_mask, :]
            factor_fieldnames = (factors.fieldname...,)
            @testset "set_factor_bounds" begin
                bounds_ranges = [range(b[1], b[2], 5) for b in factors.default_dist_params]
                new_bounds = Tuple.(sort.(rand.(bounds_ranges, 2)))
                mode_ranges = [range(nb[1], nb[2], 5) for nb in new_bounds]
                new_modes = (rand.(mode_ranges))
                new_dist_params = [(b[1], b[2], p) for (b, p) in zip(new_bounds, new_modes)]
                dom = set_factor_bounds(dom; NamedTuple{factor_fieldnames}(new_dist_params)...)
                scens = ADRIA.sample_guided(dom, num_samples)

                _test_bounds(scens, factor_mask, new_bounds)
            end

            @testset "get_default_dist_params" begin
                new_bounds = ADRIA.get_attr.([dom], factor_fieldnames, [:default_dist_params])
                dom = set_factor_bounds(dom; NamedTuple{factor_fieldnames}(new_bounds)...)

                factor_params = dom.model[ms.fieldname.∈[factor_fieldnames]][1]
                @test all(factor_params.dist_params .== factor_params.default_dist_params)

                scens = ADRIA.sample(dom, num_samples)
                _test_bounds(scens, factor_mask, new_bounds)
            end
        end

        @testset "DiscreteTriangularDist distributions" begin
            factor_mask = ms.component .∈ [test_components] .&& ms.dist .== ADRIA.DiscreteTriangularDist
            factors = ms[factor_mask, :]
            factor_fieldnames = (factors.fieldname...,)
            @testset "set_factor_bounds" begin
                bounds_ranges = [b[1]:b[2] for b in factors.default_dist_params]
                new_bounds = Tuple.(sort.(rand.(bounds_ranges, 2)))
                new_mode_ranges = [nb[1]:nb[2] for nb in new_bounds]
                new_peaks = (rand.(new_mode_ranges))
                new_dist_params = [(b[1], b[2], p) for (b, p) in zip(new_bounds, new_peaks)]
                dom = set_factor_bounds(dom; NamedTuple{factor_fieldnames}(new_dist_params)...)
                scens = ADRIA.sample_guided(dom, num_samples)

                _test_bounds(scens, factor_mask, new_bounds)
            end

            @testset "get_default_dist_params" begin
                new_bounds = ADRIA.get_attr.([dom], factor_fieldnames, [:default_dist_params])
                dom = set_factor_bounds(dom; NamedTuple{factor_fieldnames}(new_bounds)...)

                factor_params = dom.model[ms.fieldname.∈[factor_fieldnames]][1]
                @test all(factor_params.dist_params .== factor_params.default_dist_params)

                scens = ADRIA.sample(dom, num_samples)
                _test_bounds(scens, factor_mask, new_bounds)
            end
        end
    end
end
