using Test
using ADRIA

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

function _with_scale_ub(spec, fn::Symbol, ub::Float64)
    new_spec = copy(spec)
    new_spec[!, :dist_params] = [
        name == fn ? (dp[1], ub) : dp for
        (name, dp) in zip(spec.fieldname, spec.dist_params)
    ]
    return new_spec
end

@testset "Size class growth limits" begin
    dom = deepcopy(ADRIA_DOM_45)
    spec = ADRIA.model_spec(dom)
    terms = ADRIA._growth_limit_terms(spec)
    fg_names = string.(ADRIA.functional_group_names())

    cb, g, _ = first(sort(collect(keys(terms))))
    scale_fn = Symbol("linear_extension_scale_cb_group_$(cb)_$(fg_names[g])")
    scale_row = spec[spec.fieldname .== scale_fn, :]
    lb, ub = scale_row.dist_params[1][1:2]

    @testset "Domain creation leaves the worst case within the class widths" begin
        @test !isempty(terms)
        @test all(ADRIA._growth_excess(t) <= 1e-12 for t in values(terms))
    end

    @testset "Scale cap is the tightest class" begin
        caps = ADRIA._scale_ub_caps(terms)
        expected = minimum(
            t.width / (t.le_ub * t.gmax) for ((c, f, _), t) in terms if (c, f) == (cb, g)
        )
        @test caps[(cb, g)] ≈ expected
    end

    @testset "Bounds that skip a size class are rejected and leave the model untouched" begin
        before = ADRIA.get_bounds(dom, scale_fn)
        @test_throws ErrorException ADRIA.set_factor_bounds!(dom, scale_fn, (lb, 50.0))
        @test ADRIA.get_bounds(dom, scale_fn) == before
    end

    @testset "Bounds within the limit are accepted" begin
        @test ADRIA.set_factor_bounds!(dom, scale_fn, (lb, (lb + ub) / 2)) isa ADRIA.Domain
        @test ADRIA.set_factor_bounds!(dom, scale_fn, (lb, ub)) isa ADRIA.Domain
    end

    @testset "Violations that already existed are tolerated, worse ones are not" begin
        existing = _with_scale_ub(spec, scale_fn, 50.0)
        @test isnothing(ADRIA._validate_growth_limits(existing, existing))
        @test_throws ErrorException ADRIA._validate_growth_limits(
            existing, _with_scale_ub(spec, scale_fn, 100.0)
        )
    end
end
