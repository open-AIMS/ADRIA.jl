using Test

using ADRIA
using ADRIA: At

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

@testset "MC log matches target location sets" begin
    dom = deepcopy(ADRIA_DOM_45)
    rand_mc_locs = rand(dom.loc_ids, 70, 2)

    locs_1 = rand_mc_locs[:, 1]
    locs_2 = rand_mc_locs[:, 2]

    weight_1 = 0.4
    weight_2 = 0.6

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

    excl_locs_1 = setdiff(locs_1, locs_2)
    mc_log_exc1 = if !isempty(excl_locs_1)
        dropdims(
            sum(
                rs.mc_log[locations=dom.loc_ids .∈ [excl_locs_1]];
                dims=(:coral_id, :locations)
            );
            dims=(:coral_id, :locations)
        )
    else
        []
    end

    excl_locs_2 = setdiff(locs_2, locs_1)
    mc_log_exc2 = if !isempty(excl_locs_2)
        dropdims(
            sum(
                rs.mc_log[locations=dom.loc_ids .∈ [excl_locs_2]];
                dims=(:coral_id, :locations)
            );
            dims=(:coral_id, :locations)
        )
    else
        []
    end

    both_locs = intersect(locs_1, locs_2)
    mc_log_both = dropdims(
        sum(rs.mc_log[locations=dom.loc_ids .∈ [both_locs]]; dims=(:coral_id, :locations));
        dims=(:coral_id, :locations)
    )

    for s in 1:num_samples
        mc_log = mc_log_both[scenarios=At(s)]

        @test all(mc_log .<= N_mc[s] .|| mc_log .≈ N_mc[s]) ||
            "Scenario $s: combined MC log exceeds N_mc_settlers"

        # For locations exclusive to set 1 (never targeted by set 2), logged MC
        # must not exceed weight_1 * N_mc — their full allocated share
        if !isempty(excl_locs_1)
            excl_log_1 = mc_log_exc1[scenarios=At(s)]
            budget_1 = weight_1 * N_mc[s]
            @test all(excl_log_1 .<= budget_1 .|| excl_log_1 .≈ budget_1) ||
                "Scenario $s: MC log for locations exclusive to set 1 exceeds weight_1 * N_mc"
        end

        # For locations exclusive to set 2 (never targeted by set 1), logged MC
        # must not exceed weight_2 * N_mc — their full allocated share
        if !isempty(excl_locs_2)
            excl_log_2 = mc_log_exc2[scenarios=At(s)]
            budget_2 = weight_2 * N_mc[s]
            @test all(excl_log_2 .<= budget_2 .|| excl_log_2 .≈ budget_2) ||
                "Scenario $s: MC log for locations exclusive to set 2 exceeds weight_2 * N_mc"
        end
    end
end
