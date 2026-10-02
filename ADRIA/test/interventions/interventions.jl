if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

@testset "Target locations for coral aquaculture and fogging" begin
    target_mask = ADRIA_DOM_45.loc_ids .∈ [ADRIA_DOM_45.loc_ids[(end - 9):end]]

    target_locs = ADRIA_DOM_45.loc_ids[target_mask]
    ADRIA.set_CAq_target_locations!(ADRIA_DOM_45, [(weight=1.0, target_locs=target_locs)])
    ADRIA.set_Fog_target_locations!(ADRIA_DOM_45, target_locs)

    num_samples = 2
    scens = ADRIA.sample(ADRIA_DOM_45, num_samples)

    rs_raw = ADRIA.run_model(ADRIA_DOM_45, scens[1, :])

    no_CAq = (dropdims(sum(rs_raw.CAq_log; dims=(1, 2)); dims=(1, 2)) .> 0.0)[.!target_mask]
    @test all(no_CAq .== 0.0)

    no_Fog = (dropdims(sum(rs_raw.Fog_log; dims=1); dims=1) .> 0.0)[.!target_mask]
    @test all(no_Fog .== 0.0)
end
