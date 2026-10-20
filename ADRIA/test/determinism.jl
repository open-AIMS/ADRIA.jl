using Test
using ADRIA

if !@isdefined(TEST_DOMAIN_PATH)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DATA_DIR = joinpath(ADRIA_DIR, "test", "data")
    const TEST_DOMAIN_PATH = joinpath(TEST_DATA_DIR, "Test_domain")
end

function _run_twice(dom, samples)
    prev_output_dir = get(ENV, "ADRIA_OUTPUT_DIR", nothing)
    try
        # Separate output directories so back-to-back runs cannot collide
        ENV["ADRIA_OUTPUT_DIR"] = mktempdir()
        rs_1 = ADRIA.run_scenarios(dom, samples, "45")
        ENV["ADRIA_OUTPUT_DIR"] = mktempdir()
        rs_2 = ADRIA.run_scenarios(dom, samples, "45")

        return rs_1, rs_2
    finally
        isnothing(prev_output_dir) ? delete!(ENV, "ADRIA_OUTPUT_DIR") :
        (ENV["ADRIA_OUTPUT_DIR"] = prev_output_dir)
    end
end

@testset "Consecutive runs are deterministic" begin
    dom = @isdefined(TEST_DOM) ? TEST_DOM : ADRIA.load_domain(TEST_DOMAIN_PATH, "45")

    # Sample once, so only the model runs are compared
    cases = [
        ("counterfactual", ADRIA.sample_cf(dom, 8)),
        ("unguided", ADRIA.sample_unguided(dom, 8)),
        ("guided", ADRIA.sample_guided(dom, 16))
    ]

    for (name, samples) in cases
        @testset "$(name)" begin
            rs_1, rs_2 = _run_twice(dom, samples)

            cover_1 = Array(ADRIA.metrics.scenario_total_cover(rs_1))
            cover_2 = Array(ADRIA.metrics.scenario_total_cover(rs_2))

            @test isequal(cover_1, cover_2)
            @test isequal(Array(rs_1.ranks), Array(rs_2.ranks))
        end
    end
end

@testset "Scenario state does not leak between runs" begin
    dom = @isdefined(TEST_DOM) ? TEST_DOM : ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
    samples = ADRIA.sample_cf(dom, 4)

    @testset "to_coral_spec leaves the cached coral spec untouched" begin
        before = deepcopy(ADRIA.coral_spec().params)
        ADRIA.to_coral_spec(samples[1, :])
        @test isequal(ADRIA.coral_spec().params, before)
    end

    @testset "Reused sim_cache matches a fresh one" begin
        n_locs = dom.coral_growth.n_locs
        n_groups = dom.coral_growth.n_groups
        n_sizes = dom.coral_growth.n_sizes
        edges = ADRIA.bin_edges(; unit=:m).data
        make_fgs() = [
            ADRIA.FunctionalGroup.(
                eachrow(edges[:, 1:(end - 1)]),
                eachrow(edges[:, 2:end]),
                eachrow(zeros(n_groups, n_sizes))
            )
            for _ = 1:n_locs
        ]

        ref = ADRIA.run_model(dom, samples[2, :], make_fgs()).raw

        # Run a different scenario first so the reused cache holds stale state
        fgs = make_fgs()
        cache = ADRIA.setup_cache(dom)
        ADRIA.run_model(dom, samples[1, :], fgs; sim_cache=cache)
        reused = ADRIA.run_model(dom, samples[2, :], fgs; sim_cache=cache).raw

        @test isequal(reused, ref)
    end
end
