using Test
using ADRIA

if !@isdefined(TEST_DOMAIN_PATH)
    const TEST_DOMAIN_PATH = joinpath(pkgdir(ADRIA), "test", "data", "Test_domain")
end

@testset "Low-DHW bleaching reduction in scenario runs" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
    scens = ADRIA.sample(dom, 2)[[1, 1], :]  # two copies of the same scenario

    # Scenario 2 removes bleaching mortality wherever effective DHW is below 12 DHW
    scens.low_dhw_bleach_reduction[2] = 1.0
    scens.low_dhw_bleach_threshold[2] = 12.0

    # The scenarios differ in these factors, so pin the seed: otherwise the derived seed
    # would differ too and the runs would not be comparable
    rs = withenv("ADRIA_RNG_SEED" => "20260930") do
        ADRIA.run_scenarios(dom, scens, "45")
    end

    cover = ADRIA.metrics.scenario_total_cover(rs)
    @test all(isfinite, cover.data)
    @test cover[scenarios = 1].data != cover[scenarios = 2].data
    # Less bleaching loss leaves more coral
    @test sum(cover[scenarios = 2].data) > sum(cover[scenarios = 1].data)
end
