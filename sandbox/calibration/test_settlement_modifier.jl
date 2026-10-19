"""Adapter test for the opt-in reef-habitat settlement modifier (COTS_SETTLEMENT_MODIFIER_CSV;
crash_mechanism_protocol.md, 2026-10-19, wave exposure). Run:
julia --project=sandbox sandbox/calibration/test_settlement_modifier.jl"""

using Test
using ADRIA
using COTSMod
const MVector = COTSMod.MVector

function write_table(path, rows)
    open(path, "w") do io
        println(io, "site_id,multiplier")
        foreach(row -> println(io, join(row, ",")), rows)
    end
    return path
end

@testset "Settlement-modifier adapter" begin
    dir = mktempdir()
    sites = ["A", "B"]
    m = ADRIA.load_cots_settlement_modifier(write_table(joinpath(dir, "s.csv"), [("B", 0.25), ("A", 1.0)]), sites)
    @test m == [1.0, 0.25]

    params = COTSParams(separated_recruitment=true)
    models = [COTSHuman(MVector{3,Float64}(2.0, 4.0, 10.0), 1.0, params) for _ in 1:2]
    for model in models
        model.last_settlement_gate = 0.8
        model.last_background_immigration = 2.0
        model.last_local_fecundity = 7.0
    end
    ADRIA.apply_cots_settlement_modifier!(models, m)
    @test [model.last_settlement_gate for model in models] == [0.8, 0.2]
    @test [model.N[1] for model in models] == [2.0, 0.5]                 # settled background immigration
    @test [model.last_background_immigration for model in models] == [2.0, 0.5]
    @test [model.N[2] for model in models] == [4.0, 4.0]                 # juveniles and adults untouched
    @test [model.N[3] for model in models] == [10.0, 10.0]
    @test [model.last_local_fecundity for model in models] == [7.0, 7.0] # production untouched

    # External supply after the modifier settles at the reduced gate (exactly once).
    ADRIA.COTSMod.apply_external_larval_supply!(models, [10.0, 10.0], [10.0, 10.0])
    @test [model.last_external_immigration for model in models] ≈ [8.0, 2.0]

    @test_throws DimensionMismatch ADRIA.apply_cots_settlement_modifier!(models, [1.0])
    legacy = [COTSHuman(MVector{3,Float64}(2.0, 4.0, 10.0), 1.0, COTSParams(separated_recruitment=false))]
    @test_throws ErrorException ADRIA.apply_cots_settlement_modifier!(legacy, [0.5])
    @test_throws ErrorException ADRIA.load_cots_settlement_modifier(
        write_table(joinpath(dir, "bad_site.csv"), [("A", 1.0), ("Z", 1.0)]), sites)
    @test_throws ErrorException ADRIA.load_cots_settlement_modifier(
        write_table(joinpath(dir, "missing_site.csv"), [("A", 1.0)]), sites)
    @test_throws ErrorException ADRIA.load_cots_settlement_modifier(
        write_table(joinpath(dir, "dup.csv"), [("A", 1.0), ("A", 0.5), ("B", 1.0)]), sites)
    @test_throws ErrorException ADRIA.load_cots_settlement_modifier(
        write_table(joinpath(dir, "bad_value.csv"), [("A", 1.5), ("B", 1.0)]), sites)
    @test_throws ErrorException ADRIA.load_cots_settlement_modifier(joinpath(dir, "missing.csv"), sites)
end
