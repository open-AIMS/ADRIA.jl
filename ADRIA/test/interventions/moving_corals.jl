using Test

using ADRIA
using ADRIA: At

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

@testset "MC log matches target location sets" begin
    dom = deepcopy(ADRIA_DOM_45)

    n_locs = length(dom.loc_ids)
    half = n_locs ÷ 2
    locs_1 = dom.loc_ids[1:half]
    locs_2 = dom.loc_ids[(half + 1):end]

    weight_1 = TEST_TARGET_WEIGHT_1
    weight_2 = TEST_TARGET_WEIGHT_2

    ADRIA.set_mc_target_locations!(
        dom, [(weight=weight_1, target_locs=locs_1), (weight=weight_2, target_locs=locs_2)]
    )

    # Fix all MC-relevant intervention params so every sampled scenario
    # guarantees active moving corals: non-zero settlers, starts year 1,
    # runs the full simulation, deploys every year, periodic strategy.
    ADRIA.fix_factor!(
        dom;
        N_mc_settlers=500_000.0,
        mc_year_start=1.0,
        mc_years=75.0,
        mc_deployment_freq=1.0,
        mc_strategy=1.0
    )

    num_samples = 4
    scens = ADRIA.sample_guided(dom, num_samples)
    rs = ADRIA.run_scenarios(dom, scens, "45")

    # Total MC settlers requested per scenario
    N_mc = vec(scens.N_mc_settlers)

    mc_log_1 = dropdims(
        sum(rs.mc_log[locations=dom.loc_ids .∈ [locs_1]]; dims=(:coral_id, :locations));
        dims=(:coral_id, :locations)
    )
    mc_log_2 = dropdims(
        sum(rs.mc_log[locations=dom.loc_ids .∈ [locs_2]]; dims=(:coral_id, :locations));
        dims=(:coral_id, :locations)
    )

    # mc_log is persisted as Float32; use Float32 rtol regardless of in-memory eltype
    fp32_rtol = sqrt(eps(Float32))
    for s in 1:num_samples
        log_1 = mc_log_1[scenarios=At(s)]
        log_2 = mc_log_2[scenarios=At(s)]
        budget_1 = weight_1 * N_mc[s]
        budget_2 = weight_2 * N_mc[s]
        @test all(log_1 .<= budget_1 .|| isapprox.(log_1, budget_1; rtol=fp32_rtol))
        @test all(log_2 .<= budget_2 .|| isapprox.(log_2, budget_2; rtol=fp32_rtol))
    end
end
