using Test
using ADRIA

if !@isdefined(TEST_DOMAIN_PATH)
    const TEST_DOMAIN_PATH = joinpath(pkgdir(ADRIA), "test", "data", "Test_domain")
end

@testset "setup() keeps an ADRIA_OUTPUT_DIR set before run_scenarios" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
    scens = ADRIA.sample(dom, 4)

    out_dir = mktempdir()
    work_dir = mktempdir()  # no config.toml here

    # Unset what `setup()` sets, as if it had not run yet in this process, so the
    # `setup()` call inside `run_scenarios` applies the defaults. `withenv` restores them.
    setup_vars = (
        "ADRIA_HAS_SETUP", "ADRIA_THRESHOLD", "ADRIA_DEBUG", "ADRIA_LOG_DHW_TOLS",
        "ADRIA_LOG_COVER"
    )
    rs = cd(work_dir) do
        withenv("ADRIA_OUTPUT_DIR" => out_dir, (v => nothing for v in setup_vars)...) do
            res = ADRIA.run_scenarios(dom, scens, "45"; show_progress=false)
            @test ENV["ADRIA_OUTPUT_DIR"] == out_dir
            res
        end
    end

    @test isdir(joinpath(out_dir, ADRIA.store_name(rs)))
    @test !ispath(joinpath(work_dir, "Outputs"))
end
