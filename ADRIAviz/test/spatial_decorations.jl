using Test
using ADRIAviz.viz

# =============================================================================
# Lightweight tests — no backend required
# =============================================================================

@testset "validate_extent" begin
    @testset "Valid extents pass" begin
        @test isnothing(validate_extent((140.0, 155.0, -25.0, -10.0)))
        @test isnothing(validate_extent((-180.0, 180.0, -85.0, 85.0)))
    end

    @testset "NaN values throw" begin
        @test_throws ErrorException validate_extent((NaN, 155.0, -25.0, -10.0))
        @test_throws ErrorException validate_extent((140.0, NaN, -25.0, -10.0))
        @test_throws ErrorException validate_extent((140.0, 155.0, NaN, -10.0))
        @test_throws ErrorException validate_extent((140.0, 155.0, -25.0, NaN))
        @test_throws ErrorException validate_extent((NaN, NaN, NaN, NaN))
    end

    @testset "Degenerate longitude throws" begin
        @test_throws ErrorException validate_extent((155.0, 140.0, -25.0, -10.0))
        @test_throws ErrorException validate_extent((145.0, 145.0, -25.0, -10.0))
    end

    @testset "Degenerate latitude throws" begin
        @test_throws ErrorException validate_extent((140.0, 155.0, -10.0, -25.0))
        @test_throws ErrorException validate_extent((140.0, 155.0, -20.0, -20.0))
    end

    @testset "Wrong tuple length throws" begin
        @test_throws ErrorException validate_extent((140.0, 155.0, -25.0))
        @test_throws ErrorException validate_extent((140.0, 155.0, -25.0, -10.0, 0.0))
    end

    @testset "Near-pole extent warns" begin
        # mean latitude must exceed the 85 deg threshold in validate_extent
        # for the warning to fire; (80, 89) has a mean of 84.5 and never warns.
        @test_logs (:warn,) validate_extent((0.0, 10.0, 82.0, 89.0))
    end
end

@testset "_calc_gridsize" begin
    @test _calc_gridsize(1) == (1, 1)
    @test _calc_gridsize(2) == (1, 2)
    # n=3: ceil(sqrt(3))=2, cols=ceil(3/2)=2 -> (2,2) landscape
    @test _calc_gridsize(3) == (2, 2)
    @test _calc_gridsize(4) == (2, 2)
    # n=5: ceil(sqrt(5))=3, cols=ceil(5/3)=2 -> swap to (2,3) landscape
    n_r, n_c = _calc_gridsize(5)
    @test n_r <= n_c   # landscape: cols >= rows
    @test n_r * n_c >= 5
    # n=6: (2,3)
    @test _calc_gridsize(6) == (2, 3)
    # n>=7: square-ish
    n_r7, n_c7 = _calc_gridsize(7)
    @test n_r7 * n_c7 >= 7
    @test _calc_gridsize(9) == (3, 3)
    @test _calc_gridsize(16) == (4, 4)
    # Error on n < 1
    @test_throws ErrorException _calc_gridsize(0)
    @test_throws ErrorException _calc_gridsize(-1)
end

# =============================================================================
# Backend-dependent tests — require Makie
# =============================================================================

if get(ENV, "ADRIA_RUN_VIZ_TESTS", "0") == "1"
    # _figure_size and _adaptive_gap live in the Makie extension, which requires
    # all three of Makie, GeoMakie, and GraphMakie to be loaded to trigger.
    # Package extensions aren't `using`-able by name; fetch the module handle
    # via Base.get_extension instead.
    using CairoMakie, GeoMakie, GraphMakie
    _makie_ext = Base.get_extension(ADRIAviz, :ADRIAvizMakieExt)
    _figure_size = _makie_ext._figure_size
    _adaptive_gap = _makie_ext._adaptive_gap

    @testset "_figure_size" begin
        w, h = _figure_size(1, 1)
        @test w == 1000 && h == 800   # single panel landscape

        # Multi-panel: 600*ncols + 20*(ncols-1) wide, 500*nrows + 20*(nrows-1) tall
        w2, h2 = _figure_size(1, 2)
        @test w2 == 2 * 600 + 20
        @test h2 == 500

        w3, h3 = _figure_size(2, 3)
        @test w3 == 3 * 600 + 2 * 20
        @test h3 == 2 * 500 + 20

        # All sizes are positive
        for nr = 1:4, nc = 1:4
            w, h = _figure_size(nr, nc)
            @test w > 0 && h > 0
        end
    end

    @testset "_adaptive_gap" begin
        rg1, cg1 = _adaptive_gap(1)
        rg4, cg4 = _adaptive_gap(4)
        rg9, cg9 = _adaptive_gap(9)

        # Gaps are non-negative and shrink or stay flat as panel count grows
        @test rg1 >= rg4 >= rg9
        @test cg1 >= cg4 >= cg9
        # Minimum floor constraints
        @test rg9 >= 10
        @test cg9 >= 8
    end
end

# =============================================================================
# Integration test — requires domain data and Makie backend
# =============================================================================

if get(ENV, "ADRIA_RUN_MAKIE_SPATIAL_INTEGRATION", "0") == "1"
    dom_path = get(ENV, "ADRIA_TEST_DOMAIN", "")
    if isempty(dom_path) || !isdir(dom_path)
        @warn "ADRIA_RUN_MAKIE_SPATIAL_INTEGRATION=1 but ADRIA_TEST_DOMAIN not set or invalid; skipping"
    else
        using CairoMakie, GeoMakie, GraphMakie, ADRIA

        @testset "Spatial map decorations (integration)" begin
            dom = ADRIA.load_domain(dom_path)

            base_data = dom.loc_data.k .* 100.0

            @testset "map with coastlines only" begin
                f = ADRIA.viz.map(
                    dom, base_data; opts=Dict{Symbol,Any}(
                        :show_coastlines => true
                    )
                )
                @test f isa Makie.Figure
            end

            @testset "map with all decorations" begin
                f = ADRIA.viz.map(
                    dom,
                    base_data;
                    opts=Dict{Symbol,Any}(
                        :show_coastlines => true,
                        :show_coastal_places => true,
                        :show_scale_bar => true,
                        :show_north_arrow => true
                    )
                )
                @test f isa Makie.Figure
            end

            @testset "map decorations off by default" begin
                f = ADRIA.viz.map(dom, base_data)
                @test f isa Makie.Figure
            end
        end
    end
end
