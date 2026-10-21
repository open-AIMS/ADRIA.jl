"""Adapter tests for the opt-in projection mechanisms (crash_mechanism_protocol.md, 2026-10-21):
M1 forward effort culling, M2 projection-start reset, M3 coral-coupled outside supply factor and
M5 outside coral larvae in `settler_cover`. M4 (thermal acclimation) is a one-line min() in the
scenario loop and is checked by the projection runs. Run:
julia --project=sandbox sandbox/calibration/test_projection_mechanisms.jl"""

using Test
using ADRIA
using COTSMod
using Random, SparseArrays
const MVector = COTSMod.MVector

model(r, j, a) = COTSHuman(MVector{3,Float64}(r, j, a), 1.0, COTSParams())
stock(ms, areas) = sum(m.N[3] * a for (m, a) in zip(ms, areas))

@testset "Projection mechanisms" begin
    groups = [[1, 2, 3], [4]]
    areas = [2.0, 1.0, 1.0, 5.0]

    @testset "M1 catch rate" begin
        @test ADRIA.cots_cull_cpue(0.0) == 0.0
        # Inverse of D = 7.874 CPUE^0.549: CPUE 1 at D = 7.874; capped at 20.
        @test ADRIA.cots_cull_cpue(7.874) ≈ 1.0
        @test ADRIA.cots_cull_cpue(7.874 * 2.0^0.549) ≈ 2.0
        @test ADRIA.cots_cull_cpue(1e4) == 20.0
        @test ADRIA.cots_cull_cpue(1e4; cpue_max=50.0) == 50.0
    end

    @testset "M1 effort requests" begin
        ms = [model(0.0, 0.0, 10.0), model(0.0, 0.0, 0.0), model(0.0, 0.0, 20.0), model(0.0, 0.0, 1.0)]
        D1 = (10.0 * 2 + 20.0) / 4.0
        r = ADRIA.effort_cull_requests(ms, [5.0, 0.0], groups, areas)
        removed = 5.0 * ADRIA.cots_cull_cpue(D1)
        @test sum(r .* areas) ≈ removed
        @test r[2] == 0.0 && r[4] == 0.0                      # empty site, zero-effort reef
        @test r[1] / r[3] ≈ 10.0 / 20.0                        # proportional to adults
        # Effort larger than needed clears the reef, never more than the stock.
        big = ADRIA.effort_cull_requests(ms, [1e6, 0.0], groups, areas)
        @test sum(big .* areas) ≈ 40.0
        @test all(big .<= [m.N[3] for m in ms] .+ 1e-12)
        # Through apply_cots_effort_cull!: removals are logged and added to earlier removals.
        log = zeros(1, 4, 4)
        log[1, 2, :] .= 0.5
        before = stock(ms, areas)
        ADRIA.apply_cots_effort_cull!(ms, [5.0, 0.0], groups, areas, log, 1, :pre)
        @test stock(ms, areas) ≈ before - removed
        @test sum((log[1, 2, :] .- 0.5) .* areas) ≈ removed
        @test all(log[1, 4, :] .== 0.0)
    end

    @testset "M1 targeted hours" begin
        ms = [model(0.0, 0.0, 10.0), model(0.0, 0.0, 0.0), model(0.0, 0.0, 20.0), model(0.0, 0.0, 1.0)]
        h = ADRIA.targeted_cull_hours(ms, 100.0, groups, areas; threshold_ha=3.0)
        @test h == [100.0, 0.0]                              # reef 2 (1 adult/ha) is below threshold
        h2 = ADRIA.targeted_cull_hours(ms, 100.0, groups, areas; threshold_ha=0.5)
        @test h2 ≈ 100.0 .* [40.0, 5.0] ./ 45.0               # by stock
        @test ADRIA.targeted_cull_hours(ms, 100.0, groups, areas; threshold_ha=50.0) == [0.0, 0.0]
    end

    @testset "M1 effort table" begin
        path = tempname() * ".csv"
        write(path, "reef_id,hours_per_year\nA,10.0\nA,2.5\nC,4.0\n")
        @test ADRIA.load_cots_cull_effort(path, ["A", "B", "C"]) == [12.5, 0.0, 4.0]
        @test_throws ErrorException ADRIA.load_cots_cull_effort(path, ["A", "B"])
    end

    @testset "M2 reset" begin
        p = COTSParams(; adult_senescence=true)
        ms = [COTSHuman(MVector{3,Float64}(2.0, 4.0, a), 1.0, p) for a in (10.0, 0.0, 20.0, 1.0)]
        ages_before = copy(ms[1].adult_ages)
        factors, other = ADRIA.apply_cots_state_reset!(ms, [[1, 2, 3]], [5.0], areas)
        @test factors ≈ [5.0 / 10.0]                         # group density 40 / 4 = 10
        @test (ms[1].N[3] * 2 + ms[3].N[3]) / 4 ≈ 5.0
        @test ms[1].adult_ages ≈ ages_before .* 0.5
        @test sum(ms[1].adult_ages) ≈ ms[1].N[3]
        @test ms[1].N[1] ≈ 1.0 && ms[1].N[2] ≈ 2.0             # juveniles take the same factor
        @test other == 0.5 && ms[4].N[3] ≈ 0.5                 # others: median factor
        ms2 = [model(2.0, 4.0, a) for a in (10.0, 0.0, 20.0, 1.0)]
        ADRIA.apply_cots_state_reset!(ms2, [[1, 2, 3]], [30.0], areas; juveniles=false, others=:none)
        @test ms2[1].N[3] ≈ 30.0 && ms2[1].N[2] == 4.0 && ms2[4].N[3] == 1.0
        # A group with no adults is left unscaled.
        ms3 = [model(0.0, 1.0, 0.0) for _ in 1:4]
        f3, _ = ADRIA.apply_cots_state_reset!(ms3, [[1, 2]], [5.0], areas)
        @test f3 == [1.0] && ms3[1].N[2] == 1.0
    end

    @testset "M2 reset table" begin
        path = tempname() * ".csv"
        write(path, "site_id,group,year,target_adults_ha\ns1,G1,2026,5.0\ns2,G1,2026,5.0\ns4,G2,2026,2.5\n")
        year, gs, ts, labels = ADRIA.load_cots_state_reset(path, ["s1", "s2", "s3", "s4"])
        @test year == 2026 && gs == [[1, 2], [4]] && ts == [5.0, 2.5] && labels == ["G1", "G2"]
        bad = tempname() * ".csv"
        write(bad, "site_id,group,year,target_adults_ha\ns1,G1,2026,5.0\ns2,G1,2026,6.0\n")
        @test_throws ErrorException ADRIA.load_cots_state_reset(bad, ["s1", "s2"])
    end

    @testset "M3 coral factor" begin
        @test ADRIA.cots_external_coral_factor(0.1, 0.2) ≈ 0.5
        @test ADRIA.cots_external_coral_factor(0.3, 0.2) == 1.0
        @test ADRIA.cots_external_coral_factor(0.0, 0.2) == 0.0
        @test ADRIA.cots_external_coral_factor(0.1, 0.0) == 1.0
    end

    @testset "M5 outside coral larvae" begin
        fec = [10.0 0.0 5.0; 0.0 2.0 0.0]
        conn = sparse([0.5 0.2 0.0; 0.0 0.5 0.1; 0.1 0.0 0.5])
        space = [100.0, 100.0, 100.0]
        α, β, basal = [2.5, 2.5], [5000.0, 5000.0], [0.1, 0.1]
        base = ADRIA.settler_cover(fec, conn, space, α, β, basal, zeros(2, 3); rng=MersenneTwister(1))
        same = ADRIA.settler_cover(fec, conn, space, α, β, basal, zeros(2, 3); rng=MersenneTwister(1),
            external_larvae=nothing)
        @test same == base
        buf = zeros(2, 3)
        ADRIA.settler_cover(fec, conn, space, α, β, basal, buf; external_larvae=fill(3.0, 2, 3))
        @test buf ≈ fec * conn .+ 3.0                         # larvae added before settlement
        empty = ADRIA.settler_cover(zeros(2, 3), conn, space, α, β, basal, zeros(2, 3))
        @test all(empty .== 0.0)
    end
end
