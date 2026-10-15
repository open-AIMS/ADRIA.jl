"""Adapter test for the opt-in larval-survival index (COTS_LARVAL_SURVIVAL_CSV;
crash_mechanism_protocol.md, 2026-10-15). Run: julia --project=sandbox sandbox/calibration/test_larval_survival_index.jl"""

using Test
using ADRIA
using COTSMod
const MVector = COTSMod.MVector

function write_table(path, rows)
    open(path, "w") do io
        println(io, "site_id,year,multiplier")
        foreach(row -> println(io, join(row, ",")), rows)
    end
    return path
end

@testset "Larval-survival index adapter" begin
    dir = mktempdir()
    sites, years = ["A", "B"], [2018, 2019, 2020]
    m = ADRIA.load_cots_larval_survival(write_table(joinpath(dir, "s.csv"),
        [("A", 2019, 1.5), ("B", 2019, 0.5), ("A", 2030, 9.0)]), sites, years)
    @test size(m) == (3, 2)
    @test m[2, :] == [1.5, 0.5]
    @test all(m[[1, 3], :] .== 1.0)                    # absent site-years and out-of-window rows -> 1

    models = [COTSHuman(MVector{3,Float64}(5.0, 4.0, 10.0), 1.0, COTSParams()) for _ in 1:2]
    for (model, f) in zip(models, (8.0, 2.0))
        model.last_local_fecundity = f
    end
    before = [copy(model.N) for model in models]
    ADRIA.apply_cots_larval_survival!(models, m[2, :])
    @test [model.last_local_fecundity for model in models] == [12.0, 1.0]
    @test [model.N for model in models] == before        # only production is scaled
    ADRIA.apply_cots_larval_survival!(models, [1.0, 1.0])
    @test [model.last_local_fecundity for model in models] == [12.0, 1.0]   # identity

    @test_throws DimensionMismatch ADRIA.apply_cots_larval_survival!(models, [1.0])
    @test_throws ErrorException ADRIA.load_cots_larval_survival(
        write_table(joinpath(dir, "bad_site.csv"), [("Z", 2019, 1.0)]), sites, years)
    @test_throws ErrorException ADRIA.load_cots_larval_survival(
        write_table(joinpath(dir, "bad_value.csv"), [("A", 2019, -1.0)]), sites, years)
    @test_throws ErrorException ADRIA.load_cots_larval_survival(joinpath(dir, "missing.csv"), sites, years)
end
