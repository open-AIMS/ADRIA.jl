"""Adapter test for the opt-in observed-culling switch (COTS_CULL_REMOVALS_CSV;
crash_mechanism_protocol.md, 2026-10-14). Run: julia --project=sandbox sandbox/calibration/test_cull_removals.jl"""

using Test
using ADRIA
using COTSMod
const MVector = COTSMod.MVector

function write_table(path, rows)
    open(path, "w") do io
        println(io, "site_id,report_year,phase,juveniles_ha,adults_ha")
        for row in rows
            println(io, join(row, ","))
        end
    end
    return path
end

toy_models() = [COTSHuman(MVector{3,Float64}(5.0, 4.0, 10.0), 1.0, COTSParams()),
                COTSHuman(MVector{3,Float64}(5.0, 1.0, 2.0), 1.0, COTSParams())]

@testset "Observed culling adapter" begin
    dir = mktempdir()
    sites, years = ["A", "B"], [2012, 2013, 2014]
    path = write_table(joinpath(dir, "cull.csv"), [
        ("A", 2013, "pre", 1.0, 3.0), ("A", 2013, "post", 0.5, 2.0),
        ("A", 2013, "pre", 0.5, 1.0),                 # duplicate rows add
        ("B", 2014, "post", 9.0, 9.0),                # capped at what is present
        ("B", 2030, "pre", 1.0, 1.0)])                # outside the simulated years: ignored
    requests = ADRIA.load_cots_cull_removals(path, sites, years)
    @test size(requests) == (3, 4, 2)
    @test requests[2, :, 1] == [1.5, 4.0, 0.5, 2.0]
    @test requests[3, :, 2] == [0.0, 0.0, 9.0, 9.0]
    @test sum(requests) == 1.5 + 4.0 + 0.5 + 2.0 + 18.0

    # Pre then post within one timestep; the log records what was removed.
    models = toy_models()
    removed = zeros(3, 4, 2)
    ADRIA.apply_cots_cull!(models, requests, removed, 2, :pre)
    @test models[1].N ≈ [5.0, 2.5, 6.0]
    ADRIA.apply_cots_cull!(models, requests, removed, 2, :post)
    @test models[1].N ≈ [5.0, 2.0, 4.0]
    @test removed[2, :, 1] == [1.5, 4.0, 0.5, 2.0]
    @test models[2].N == MVector{3,Float64}(5.0, 1.0, 2.0)   # no request at B in 2013
    ADRIA.apply_cots_cull!(models, requests, removed, 3, :post)
    @test removed[3, :, 2] == [0.0, 0.0, 1.0, 2.0]             # capped
    @test models[2].N ≈ [5.0, 0.0, 0.0]                         # recruits never culled

    # A zero table is the identity.
    zero_path = write_table(joinpath(dir, "zero.csv"), [("A", 2013, "pre", 0.0, 0.0)])
    zero_models = toy_models()
    zero_requests = ADRIA.load_cots_cull_removals(zero_path, sites, years)
    ADRIA.apply_cots_cull!(zero_models, zero_requests, zeros(3, 4, 2), 2, :pre)
    @test [m.N for m in zero_models] == [m.N for m in toy_models()]

    # Validation.
    @test_throws ErrorException ADRIA.load_cots_cull_removals(
        write_table(joinpath(dir, "bad_site.csv"), [("Z", 2013, "pre", 1.0, 1.0)]), sites, years)
    @test_throws ErrorException ADRIA.load_cots_cull_removals(
        write_table(joinpath(dir, "bad_phase.csv"), [("A", 2013, "mid", 1.0, 1.0)]), sites, years)
    @test_throws ErrorException ADRIA.load_cots_cull_removals(
        write_table(joinpath(dir, "bad_value.csv"), [("A", 2013, "pre", -1.0, 1.0)]), sites, years)
    @test_throws ErrorException ADRIA.load_cots_cull_removals(joinpath(dir, "missing.csv"), sites, years)
    @test_throws ArgumentError ADRIA.apply_cots_cull!(toy_models(), requests, removed, 2, :mid)
end
