"""Adapter test for the opt-in reef-proportional cull allocation (COTS_CULL_ALLOCATION;
crash_mechanism_protocol.md, 2026-10-19). Run:
julia --project=sandbox sandbox/calibration/test_cull_allocation.jl"""

using Test
using ADRIA
using COTSMod
const MVector = COTSMod.MVector

model(j, a) = COTSHuman(MVector{3,Float64}(1.0, j, a), 1.0, COTSParams())

@testset "Reef-proportional cull allocation" begin
    # Reef 1 = sites 1-3, reef 2 = site 4. Areas in ha.
    groups = [[1, 2, 3], [4]]
    areas = [2.0, 1.0, 1.0, 5.0]
    models = [model(0.0, 10.0), model(0.0, 0.0), model(0.0, 20.0), model(0.0, 1.0)]
    # Reef 1 requests 30 COTS, all at the empty site 2 (site-matched would remove nothing).
    requests = [0.0, 30.0, 0.0, 0.0]
    r = ADRIA.reef_proportional_requests(models, requests, 3, groups, areas)
    stock = 10.0 * 2 + 0.0 + 20.0 * 1
    @test r ≈ [30 * 10 / stock, 0.0, 30 * 20 / stock, 0.0]
    @test sum(r .* areas) ≈ 30.0                              # reef total preserved
    @test all(r .<= [m.N[3] for m in models])                # fully removable (30 <= 40)

    # Request above the reef stock clears every site; reef with no COTS stays unmet.
    big = ADRIA.reef_proportional_requests(models, [0.0, 100.0, 0.0, 0.0], 3, groups, areas)
    @test sum(min.(big, [m.N[3] for m in models]) .* areas) ≈ stock
    empty = ADRIA.reef_proportional_requests([model(0.0, 0.0) for _ in 1:4], requests, 3, groups, areas)
    @test all(empty .== 0.0)

    # Through apply_cots_cull!: site mode removes nothing at site 2; reef mode removes 30 COTS.
    req3 = zeros(1, 4, 4)
    req3[1, 2, :] .= requests          # pre-spawning adults
    for (mode, expected) in ((nothing, 0.0), (groups, 30.0))
        ms = [model(0.0, 10.0), model(0.0, 0.0), model(0.0, 20.0), model(0.0, 1.0)]
        log = zeros(1, 4, 4)
        ADRIA.apply_cots_cull!(ms, req3, log, 1, :pre; groups=mode, areas=isnothing(mode) ? nothing : areas)
        @test sum(log[1, 2, :] .* areas) ≈ expected
        @test sum((m.N[3] for m in ms) .* areas) ≈ stock + 5.0 - expected
    end
    # Juveniles use N[2].
    js = [model(4.0, 0.0), model(0.0, 0.0), model(4.0, 0.0), model(0.0, 0.0)]
    rj = ADRIA.reef_proportional_requests(js, [0.0, 6.0, 0.0, 0.0], 2, groups, areas)
    @test sum(rj .* areas) ≈ 6.0 && rj[2] == 0.0
    @test_throws ArgumentError ADRIA.reef_proportional_requests(models, requests, 1, groups, areas)
    @test_throws DimensionMismatch ADRIA.reef_proportional_requests(models, requests[1:3], 3, groups, areas)
end
