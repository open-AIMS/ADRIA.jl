using Test
using Random
using ADRIA.decision.JMcDM
using ADRIA.decision:
    CAqPreferences,
    decision_matrix,
    select_locations,
    juvenile_share,
    JUVENILE_SHARE_MIN_COVER

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

@testset "Juvenile share criterion" begin
    n_groups, n_sizes = 5, 7

    # Cover of shape [groups ⋅ sizes ⋅ locations], relative to habitable area.
    # Each location gets total cover `tot` with a juvenile share `jsh`, where juveniles are
    # the size classes ADRIA defines as such (`ADRIA.juvenile_size_classes`).
    is_juv = ADRIA.juvenile_size_classes(n_sizes)
    function make_cover(tot::Vector{Float64}, jsh::Vector{Float64})
        C = zeros(n_groups, n_sizes, length(tot))
        for l in eachindex(tot)
            C[:, is_juv, l] .= tot[l] * jsh[l] / (n_groups * count(is_juv))
            C[:, .!is_juv, l] .= tot[l] * (1 - jsh[l]) / (n_groups * count(.!is_juv))
        end
        return C
    end

    @testset "Rank transform and neutral value" begin
        tot = [0.30, 0.30, 0.30, 0.30, 0.001, 0.0]
        jsh = [0.20, 0.01, 0.05, 0.05, 0.00, 0.0]
        scores = juvenile_share(make_cover(tot, jsh))

        # Eligible locations are ranked by share, scaled to (0, 1), ties share a rank
        @test scores[1:4] ≈ [3.5, 0.5, 2.0, 2.0] ./ 4
        @test all(0.0 .< scores[1:4] .< 1.0)

        # Locations below the minimum cover are neutral, however low their share
        @test tot[5] < JUVENILE_SHARE_MIN_COVER
        @test scores[5:6] == [0.5, 0.5]

        # Fewer than two eligible locations gives no ranking signal
        @test all(juvenile_share(make_cover([0.3, 0.001], [0.1, 0.0])) .== 0.5)

        # Threshold is adjustable
        @test juvenile_share(make_cover(tot, jsh); min_cover=0.0)[5] < scores[2]
    end

    @testset "Prefers adult-dominated locations" begin
        dom = deepcopy(ADRIA_DOM_45)
        n_locs = length(dom.loc_ids)

        # Juvenile share dominates; a weakly weighted second criterion avoids a location
        # being worst on every criterion, which makes COCOSO scores degenerate (NaN/Inf).
        w = ADRIA.component_params(dom.model, ADRIA.CAqCriteriaWeights)
        ADRIA.fix_factor!(dom; Dict(Symbol(f) => 0.0 for f in w.fieldname)...)
        ADRIA.fix_factor!(dom; iv_CAq_juvenile_share=0.9, iv_CAq_depth=0.1)
        sp = CAqPreferences(dom)

        jsh = collect(range(0.5, 0.01; length=n_locs))
        dm = decision_matrix(
            dom.loc_ids,
            sp.names;
            iv_CAq_juvenile_share=juvenile_share(make_cover(fill(0.3, n_locs), jsh)),
            iv_CAq_depth=mod.((1:n_locs) .* 0.618, 1.0)
        )

        # MOORA is excluded: its reference-point ranking here follows depth alone, and it
        # does the same with `iv_CAq_coral_cover` in place of this criterion.
        methods = filter(m -> m != JMcDM.MooraMethod, ADRIA.mcda_methods())
        @testset "$(method)" for method in methods
            selected = select_locations(sp, dm, method, dom.loc_ids, 3)
            @test selected == reverse(dom.loc_ids)[1:3] ||
                "Locations with the lowest juvenile share were not preferred"
        end
    end

    @testset "Zero weight disables the criterion" begin
        dom = deepcopy(ADRIA_DOM_45)
        n_locs = length(dom.loc_ids)
        rng = Random.Xoshiro(42)
        others = (
            iv_CAq_heat_stress=rand(rng, n_locs), iv_CAq_in_connectivity=rand(rng, n_locs),
            iv_CAq_out_connectivity=rand(rng, n_locs), iv_CAq_depth=rand(rng, n_locs),
            iv_CAq_coral_cover=rand(rng, n_locs)
        )
        juv = rand(rng, n_locs)

        # Exact zero, and the ~1e-17 that a fixed 0 becomes after `gamma_to_dirichlet`
        @testset "weight=$(w)" for w in (0.0, 1e-17)
            ADRIA.fix_factor!(dom; iv_CAq_juvenile_share=w)
            sp = CAqPreferences(dom)
            dm_with = decision_matrix(
                dom.loc_ids, sp.names; others..., iv_CAq_juvenile_share=juv
            )
            dm_without = decision_matrix(dom.loc_ids, sp.names; others...)

            @testset "$(method)" for method in ADRIA.mcda_methods()
                @test select_locations(sp, dm_with, method, dom.loc_ids, n_locs) ==
                    select_locations(sp, dm_without, method, dom.loc_ids, n_locs)
            end
        end
    end
end
