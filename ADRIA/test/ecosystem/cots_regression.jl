# Fixed-seed regression test of the COTS submodel (COTS on and off).
#
# The reference values were printed by the code at `main` ac75f51a, before COTS was moved or
# refactored. A refactor that is meant to leave results unchanged must keep them: sums and
# extrema agree to rtol=1e-12 (exact equality of floating-point sums can differ across CPUs
# and Julia versions; the bit-exact comparison is the out-of-repo comparison script, see
# docs/src/development/testing_cots.md).
#
# To regenerate the reference values after an intended change in results (this needs
# sign-off: see docs/src/development/testing_cots.md), run this file with the environment
# variable COTS_REGRESSION_PRINT=true and paste the printed `EXPECTED` block over the old one.
#
# Scenarios use derived per-scenario seeds unless stated. The results do not depend on the
# number of Julia threads, so the same values must hold for serial and threaded test runs.
using Random
using Test
using ADRIA

if !@isdefined(TEST_DOMAIN_PATH)
    const TEST_DOMAIN_PATH = joinpath(pkgdir(ADRIA), "test", "data", "Test_domain")
end

# Printed by this file (COTS_REGRESSION_PRINT=true) at ac75f51a: (sum, minimum, maximum)
const EXPECTED = Dict{String,NTuple{3,Float64}}(
    "td_cover_off1" => (194.04794697218858, 0.0, 0.6223124861717224),
    "td_cover_off2" => (120.91912182149827, 0.0, 0.6092947721481323),
    "td_cover_on1" => (191.06593707474713, 0.0, 0.6129801273345947),
    "td_cover_on2" => (120.30090800719336, 0.0, 0.6092978715896606),
    "td_pop_on1" => (5.948446251433873, 0.0, 0.10599999874830246),
    "td_pop_on2" => (1.8208273302040847, 0.0, 0.09001278877258301),
    "td_cond_on1" => (450.16309234686196, 0.0, 0.9999964237213135),
    "td_cond_on2" => (294.7367363628, 0.0, 0.9999994039535522),
    "td_cover_off1_seed123" => (194.03546339116838, 0.0, 0.622288703918457),
    "td_cover_on1_seed123" => (190.47171702803462, 0.0, 0.61238694190979),
    "td_pop_on1_seed123" => (6.063600710916484, 0.0, 0.10599999874830246),
    "pkg_pop1" => (925.3882593872204, 0.0, 21.92469024658203),
    "pkg_cover1" => (56.91687519138202, 0.0, 0.44789180159568787),
    "pkg_pop2" => (523.1784016289832, 0.0, 6.610030174255371),
    "pkg_cover2" => (61.98234200383456, 0.0, 0.44789180159568787),
    "pkg_pop3" => (524.0871908137464, 0.0, 6.589095115661621),
    "pkg_cover3" => (62.08844968718688, 0.0, 0.44789180159568787),
    "pkg_pop4" => (522.0162394908672, 0.0, 7.239533424377441),
    "pkg_cover4" => (61.819966158324405, 0.0, 0.44789180159568787),
    "pkg_cover_off" => (133.36539362321275, 0.0, 0.5821441411972046)
)
@testset "COTS regression (fixed seeds)" begin
    # (sum, minimum, maximum) of an array, accumulated in Float64
    summary_of(x) = (sum(Float64, x), Float64(minimum(x)), Float64(maximum(x)))
    printing = get(ENV, "COTS_REGRESSION_PRINT", "false") == "true"
    printed = Pair{String,Any}[]
    function check(name, actual; rtol=1e-12)
        if printing
            push!(printed, name => actual)
        else
            @test all(
                isapprox.(collect(actual), collect(EXPECTED[name]); rtol=rtol, atol=1e-12)
            )
        end
    end
    function run_in_scratch(dom, scens)
        return withenv("ADRIA_OUTPUT_DIR" => mktempdir()) do
            ADRIA.run_scenarios(dom, scens, "45")
        end
    end
    cover(rs, i) = rs.outcomes[:relative_cover][scenarios = i].data
    # Total relative cover over time and locations of scenario `i`
    total_cover(rs, i) = sum(Float64, cover(rs, i))

    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")

    # Rows 1-2 with COTS off and on (row 4 also with juvenile predation); the on/off pairs
    # share their parameters other than COTS factors and therefore their derived seeds.
    Random.seed!(20261022)
    base = ADRIA.sample(dom, 2)
    scens = vcat(base, base)
    scens.cots_on[3:4] .= 1
    scens.cots_juvenile_predation[4] = 1
    rs = run_in_scratch(dom, scens)

    @testset "Test_domain, derived seeds" begin
        @test !isnothing(rs.cots_pop_log)
        @test all(rs.cots_pop_log[scenarios = 1].data .== 0.0)
        @test all(rs.cots_pop_log[scenarios = 2].data .== 0.0)
        check("td_cover_off1", summary_of(cover(rs, 1)))
        check("td_cover_off2", summary_of(cover(rs, 2)))
        check("td_cover_on1", summary_of(cover(rs, 3)))
        check("td_cover_on2", summary_of(cover(rs, 4)))
        check("td_pop_on1", summary_of(rs.cots_pop_log[scenarios = 3].data))
        check("td_pop_on2", summary_of(rs.cots_pop_log[scenarios = 4].data))
        check("td_cond_on1", summary_of(rs.cots_condition_log[scenarios = 3].data))
        check("td_cond_on2", summary_of(rs.cots_condition_log[scenarios = 4].data))
    end

    @testset "COTS lowers coral cover" begin
        # The pairs differ only in COTS (same derived seed), so COTS is the only cause
        for (off, on) in ((1, 3), (2, 4))
            @test total_cover(rs, on) < total_cover(rs, off)
        end
    end

    @testset "Test_domain, fixed ADRIA_RNG_SEED" begin
        rs_seed = withenv("ADRIA_RNG_SEED" => "123") do
            run_in_scratch(dom, scens[[1, 3], :])
        end
        check("td_cover_off1_seed123", summary_of(cover(rs_seed, 1)))
        check("td_cover_on1_seed123", summary_of(cover(rs_seed, 2)))
        check("td_pop_on1_seed123", summary_of(rs_seed.cots_pop_log[scenarios = 2].data))
    end

    @testset "data package inputs: forcing modes, derived and fixed larval seeds" begin
        @isdefined(make_cots_test_package) ||
            include(joinpath(@__DIR__, "..", "io", "cots_test_package.jl"))
        pkg_dom = ADRIA.load_domain(
            make_cots_test_package(joinpath(mktempdir(), "Test_domain")), "45"
        )
        Random.seed!(20261024)
        pkg = ADRIA.sample(pkg_dom, 2)[1:1, :]
        pscens = vcat(pkg, pkg, pkg, pkg, pkg)
        pscens.cots_on[1:4] .= 1
        pscens.cots_separated_recruitment[2:4] .= 1
        pscens.cots_larval_fecundity[2:4] .= 20.0
        pscens.cots_larval_conn_mode[3:4] .= [1, 2]
        pscens.cots_larval_conn_seed[4] = 20260930  # row 5 (COTS off) also has the package
        prs = run_in_scratch(pkg_dom, pscens)
        for i in 1:4
            check("pkg_pop$i", summary_of(prs.cots_pop_log[scenarios = i].data))
            check("pkg_cover$i", summary_of(cover(prs, i)))
        end
        check("pkg_cover_off", summary_of(cover(prs, 5)))
        @test all(prs.cots_pop_log[scenarios = 5].data .== 0.0)
        # COTS lowers cover with the package inputs too (row 5 is row 1 with COTS off)
        @test total_cover(prs, 1) < total_cover(prs, 5)
    end

    if printing
        println("const EXPECTED = Dict{String,NTuple{3,Float64}}(")
        for (k, v) in printed
            println("    ", repr(k), " => ", repr(v), ",")
        end
        println(")")
    end
end
