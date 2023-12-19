using ADRIA


if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "examples", "Test_domain")
end

@testset "sample" begin
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH)
    num_samples = 32
    scens = ADRIA.sample(dom, num_samples)

    ms = ADRIA.model_spec(dom)
    constant_params = ms.is_constant .== true
    @test all(values(scens[1, constant_params]) .== values(scens[end, constant_params])) || "Constant params are not constant!"

    eco = (ms.component .== "Coral") .& .!(constant_params)
    interv = (ms.component .== "Intervention") .& .!(constant_params)
    min_x = values(ms[:, :lower_bound])
    max_x = values(ms[:, :upper_bound])

    msg = "Sampled values were not in expected bounds!"
    coral_msg = "Sampled coral values were not in expected bounds!"
    for i in 1:num_samples
        x = values(scens[i, :])

        if scens[i, :guided] > 0
            cond = min_x .<= x .<= max_x
            @test all(cond) || "$msg | $(ms[.!(cond), :]) | $(x[.!(cond)])"
            continue
        end

        # When no interventions are used, e.g., for counterfactual or unguided scenarios
        # (guided ∈ [-1, 0]) intervention parameters are set to 0 so only check ecological values
        cond = min_x[eco] .<= x[eco] .<= max_x[eco]
        @test all(cond) || "$coral_msg | $(ms[.!(cond), :]) | $(x[eco][.!(cond)])"

        # Note: Test to ensure all intervention factors are set to 0 is covered by the guided
        # sampling test below
    end
end

@testset "Flooring trick" begin
    # Check that continous sampled values are correctly mapped
    # to corresponding discrete values.

    @test ADRIA.map_to_discrete(50.7, 51) == 50

    @test ADRIA.map_to_discrete(10.9, 11) == 10

    # Test value capped to upper bound
    @test ADRIA.map_to_discrete(4, 4) == 3

    x = rand(1:0.01:20, 25)
    expect = min.(floor.(Int64, x), ceil.(Int64, x) .- 1)
    calc = ADRIA.map_to_discrete.(x, Int64.(ceil.(x)))
    msg = """
    Flooring trick failed to produce discrete values
    Expected: $(expect)
    Received: $(calc)
    """
    @test all(expect .== calc) || msg
end

@testset "Targeted sampling" begin
    @testset "Counterfactual sampling" begin
        dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH)
        num_samples = 32
        scens = ADRIA.sample_cf(dom, num_samples)

        @test all(scens.guided .== -1) || "Intervention scenarios found"

        # Get Intervention params
        interv_params = string.(ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Intervention).fieldname)

        # Ensure all interventions are deactivated (ignoring the "guided" factor)
        interv_params = String[ip for ip in interv_params if ip != "guided"]
        @test all(all.(==(0), eachcol(scens[:, interv_params]))) || "Intervention factors with values > 0 found"
    end

    @testset "Guided sampling" begin
        dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH)
        num_samples = 32
        scens = ADRIA.sample_guided(dom, num_samples)

        @test all(scens.guided .> 0) || "Non-intervention scenarios found"

        # Get Intervention params
        interv_params = string.(ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Intervention).fieldname)

        # Ignore guided
        interv_params = String[ip for ip in interv_params if ip != "guided"]

        # Ensure at least one intervention is active
        @test all(any.(>(0), eachcol(scens[:, interv_params]))) || "All intervention factors had values <= 0"

        crit = ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.CriteriaWeights)
        seed_weights = ADRIA.criteria_params(crit, (:seed, :weight)).fieldname
        fog_weights = ADRIA.criteria_params(crit, (:fog, :weight)).fieldname

        @test all(abs.(sum(Matrix(scens[:, seed_weights]); dims=2) .- 1.0) .< 10e-6) ||
            "Some seeding weights are not properly normalized."
        @test all(abs.(sum(Matrix(scens[:, fog_weights]); dims=2) .- 1.0) .< 10e-6) ||
            "Some fogging weights are not properly normalized."

    end

    @testset "Unguided sampling" begin
        dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH)
        num_samples = 32
        scens = ADRIA.sample_unguided(dom, num_samples)

        @test all(scens.guided .== 0) || "Intervention or counterfactual scenarios found"

        # Get Intervention params
        interv_params = string.(ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Intervention).fieldname)

        # Ignore guided and planning horizon
        interv_params = String[ip for ip in interv_params if ip != "plan_horizon" && ip != "guided"]

        # Ensure at least one intervention is active
        @test all(any.(>(0), eachcol(scens[:, interv_params]))) || "All intervention factors had values <= 0"
    end

    @testset "Site selection sampling" begin
        dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH)
        num_samples = 32
        scens = ADRIA.sample_site_selection(dom, num_samples)

        @test all(scens.guided .> 0) || "Intervention or counterfactual scenarios found"

        # Get Intervention params
        ms = ADRIA.model_spec(dom)
        target_params =
            string.(
                ADRIA.component_params(
                    ms,
                    [ADRIA.EnvironmentalLayer, ADRIA.Intervention, ADRIA.CriteriaWeights],
                ).fieldname
            )

        # Ignore guided
        target_params = String[ip for ip in target_params if ip != "guided"]

        # Ensure at least one intervention is active
        @test all(any.(>(0), eachcol(scens[:, target_params]))) || "All target factors had values <= 0"

        # Check that all coral parameters are set to their nominated default values
        coral_params = ADRIA.component_params(ms, ADRIA.Coral).fieldname

        @test all([all(scens[:, c] .== ms[ms.fieldname.==c, :val][1]) for c in coral_params]) || "Non-default coral parameter value found"
    end

    @testset "Set new sampling bounds" begin
        dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH)
        num_samples = 32

        # test continuous factor is sampled within specified range
        bnds = rand(2)
        ADRIA.set_factor_bounds!(
            dom, :deployed_coral_risk_tol, (minimum(bnds), maximum(bnds))
        )
        scens = ADRIA.sample_guided(dom, num_samples)
        @test (maximum(scens[:, "deployed_coral_risk_tol"]) <= maximum(bnds)) ||
            "Sampled continuous factor is outside of specified new bounds."
        @test (minimum(scens[:, "deployed_coral_risk_tol"]) >= minimum(bnds)) ||
            "Sampled continuous factor is outside of specified new bounds."

        # test discrete factor is sampled within specified range and is discrete
        bnds = rand(0.0:1000000.0, 2)
        ADRIA.set_factor_bounds!(dom, :N_seed_TA, (minimum(bnds), maximum(bnds)))
        scens = ADRIA.sample_site_selection(dom, num_samples)
        @test (maximum(scens[:, "N_seed_TA"]) <= maximum(bnds)) ||
            "Sampled discrete factor is outside of specified new bounds."
        @test (minimum(scens[:, "N_seed_TA"]) >= minimum(bnds)) ||
            "Sampled discrete factor is outside of specified new bounds."
        @test all(mod.(scens[:, "N_seed_TA"], 1.0) .== 0.0) ||
            "Sampled discrete factors are not all discrete."
    end

end
