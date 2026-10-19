"""Adapter test for the opt-in year-specific external-supply survival (COTS_EXTERNAL_SURVIVAL_CSV;
crash_mechanism_protocol.md, 2026-10-19, twenty-fourth round). Run:
julia --project=sandbox sandbox/calibration/test_external_survival.jl"""

using Test
using ADRIA

function write_table(path, rows)
    open(path, "w") do io
        println(io, "year,multiplier")
        foreach(row -> println(io, join(row, ",")), rows)
    end
    return path
end

@testset "External-supply survival loader" begin
    dir = mktempdir()
    years = [2018, 2019, 2020]
    m = ADRIA.load_cots_external_survival(write_table(joinpath(dir, "s.csv"),
        [(2019, 1.5), (2020, 0.5), (2031, 9.0)]), years)
    @test m == [1.0, 1.5, 0.5]                       # absent year -> 1; out-of-window ignored
    @test_throws ErrorException ADRIA.load_cots_external_survival(
        write_table(joinpath(dir, "dup.csv"), [(2019, 1.0), (2019, 2.0)]), years)
    @test_throws ErrorException ADRIA.load_cots_external_survival(
        write_table(joinpath(dir, "neg.csv"), [(2019, -1.0)]), years)
    @test_throws ErrorException ADRIA.load_cots_external_survival(joinpath(dir, "missing.csv"), years)
end
