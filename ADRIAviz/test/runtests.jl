using Test
using ADRIAviz

@testset "ADRIAviz" begin
    pkg_ts = @testset "Package loads without backend" begin
        @test isdefined(ADRIAviz, :COLORS)
        @test isdefined(ADRIAviz, :_get_scenario_groups)
        @test isdefined(ADRIAviz, :_scenario_types)
        @test isdefined(ADRIAviz, :outcome_probability)
        @test isdefined(ADRIAviz, :OPT_TYPE)
        @test isdefined(ADRIAviz, :_time_labels)
    end

    makie_ts = @testset "No Makie in loaded modules" begin
        makie_loaded = any(
            contains(k.name, "Makie") for k in keys(Base.loaded_modules)
        )
        @test !makie_loaded
    end

    let tc1 = Test.get_test_counts(pkg_ts), tc2 = Test.get_test_counts(makie_ts)
        if tc1.fails + tc1.errors + tc2.fails + tc2.errors > 0
            @error "ADRIAviz: package integrity checks failed — aborting viz tests"
            exit(1)
        end
    end

    # Spatial utilities tests (lightweight, no backend required)
    include("spatial_utils.jl")

    # Spatial decoration tests: validate_extent, _calc_gridsize (lightweight);
    # _figure_size/_adaptive_gap (backend-gated); integration (domain-gated)
    include("spatial_decorations.jl")

    if get(ENV, "ADRIA_RUN_VIZ_TESTS", "0") == "1"
        include("annotated_outcomes.jl")
        include("taxa_dynamics.jl")
        include("sensitivity.jl")
        include("clustering.jl")
        # spatial.jl deferred -- requires real domain geometry
    end

    # Plotly backend tests.
    # Gate: ADRIA_RUN_PLOTLY_TESTS=1
    # Constraint: must run in a separate process / before any Makie `using`
    #             statement (Plotly and Makie extensions are mutually exclusive).
    if get(ENV, "ADRIA_RUN_PLOTLY_TESTS", "0") == "1"
        include("plotly.jl")
    end

    # Plotly end-to-end integration test.
    # Gate: ADRIA_RUN_PLOTLY_INTEGRATION=1
    # Requires:
    #   - ADRIA_TEST_DOMAIN set to a valid domain directory path
    #   - Heavy packages not in this Project.toml: MLJ, SIRUS, ADRIAanalysis
    #     These must be available in the active Julia project environment.
    # Runs examples/plotly_viz_check.jl and asserts all figures succeed.
    if get(ENV, "ADRIA_RUN_PLOTLY_INTEGRATION", "0") == "1"
        dom_path = get(ENV, "ADRIA_TEST_DOMAIN", "")
        if isempty(dom_path) || !isdir(dom_path)
            @warn "ADRIA_RUN_PLOTLY_INTEGRATION=1 but ADRIA_TEST_DOMAIN is not set " *
                "or invalid; skipping integration test"
        else
            check_script = joinpath(
                @__DIR__, "..", "..", "ADRIA", "docs", "scripts", "plotly_viz_check.jl"
            )
            @testset "Plotly integration (end-to-end)" begin
                include(check_script)
                n_fail = count(r -> !r[2], RESULTS)
                @test n_fail == 0
            end
        end
    end
end
