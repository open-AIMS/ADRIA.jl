if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

@testset "Target locations for seeding and fogging" begin
    target_mask = ADRIA_DOM_45.loc_ids .∈ [ADRIA_DOM_45.loc_ids[(end - 9):end]]

    target_locs = ADRIA_DOM_45.loc_ids[target_mask]
    ADRIA.set_seed_target_locations!(ADRIA_DOM_45, target_locs)
    ADRIA.set_fog_target_locations!(ADRIA_DOM_45, target_locs)

    num_samples = 2
    scens = ADRIA.sample(ADRIA_DOM_45, num_samples)

    rs_raw = ADRIA.run_model(ADRIA_DOM_45, scens[1, :])

    no_seed = (dropdims(sum(rs_raw.seed_log; dims=(1, 2)); dims=(1, 2)) .> 0.0)[.!target_mask]
    @test all(no_seed .== 0.0)

    no_fog = (dropdims(sum(rs_raw.fog_log; dims=1); dims=1) .> 0.0)[.!target_mask]
    @test all(no_fog .== 0.0)
end
