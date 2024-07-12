using Test
using ADRIA

include("../test_helpers.jl")

@testset "rankings" begin
    rs = get_example_resultset()

    @testset "Standalone ranking" begin
        dom = get_example_domain()
        scens = ADRIA.sample_selection(dom, 8)

        # Use rank_locations to get ranks
        n_corals = 1_000_000  # number of corals to deploy
        ranks = ADRIA.decision.rank_locations(dom, n_corals, scens)
        scores = ADRIA.decision.selection_score(ranks, :seed)

        @test all(0.0 .<= scores .< 1.0) || "Selection scores are not ∈ [0, 1]"
    end

    @testset "Scenario-based rankings" begin
        rs = get_example_resultset()

        # rs.ranks
        n_seeded_locs = ADRIA.metrics.n_seed_locations(rs)
        @test any(n_seeded_locs .> 0) || "No locations were found to be seeded"

        n_seeded_locs = ADRIA.metrics.n_seed_locations(rs; timesteps=1:10, scenarios=1:2)
        @test any(n_seeded_locs .> 0) || "No locations were found to be seeded in subset"

        n_fogged_locs = ADRIA.metrics.n_fog_locations(rs)
        @test any(n_fogged_locs .> 0) || "No locations were found to be fogged"

        n_fogged_locs = ADRIA.metrics.n_fog_locations(rs; timesteps=1:10, scenarios=1:2)
        @test any(n_fogged_locs .> 0) || "No locations were found to be fogged in subset"

        seed_summary = ADRIA.decision.deployment_summary_stats(rs.ranks, :seed)
        @test any(seed_summary[:, 2] .> 0) || "Mean number of seeded locations were all zero"

        seed_summary = ADRIA.decision.deployment_summary_stats(rs.ranks, :fog)
        @test any(seed_summary[:, 2] .> 0) || "Mean number of fogged locations were all zero"

        seed_rankings = ADRIA.metrics.seed_ranks(rs)
        @test any(seed_rankings .> 0) || "No seed rankings found"

        seed_rankings = ADRIA.metrics.seed_ranks(rs; timesteps=1:10, scenarios=1:2)
        @test any(seed_rankings .> 0) || "No seed rankings found in subset"

        fog_rankings = ADRIA.metrics.fog_ranks(rs)
        @test any(fog_rankings .> 0) || "No fog rankings found"

        fog_rankings = ADRIA.metrics.fog_ranks(rs; timesteps=1:10, scenarios=1:2)
        @test any(fog_rankings .> 0) || "No fog rankings found in subset"

        # Not relevant for now
        # ADRIA.metrics.shade_ranks(rs)
        # ADRIA.metrics.shade_ranks(rs; timesteps=1:10, scenarios=1:2)

        n_locs = 5
        top_locs = ADRIA.metrics.top_n_locs(rs, n_locs; stat=mean)
        @test size(top_locs, 2) == n_locs || "Number of rank-columns do not match number of requested locations"

        top_locs = ADRIA.metrics.top_n_locs(rs, 10; stat=median, metric=ADRIA.metrics.absolute_shelter_volume)

        # ADRIA.metrics.top_n_seeded_sites(rs, 10);

        loc_rankings = ADRIA.decision.selection_ranks(rs.ranks, :seed)
        loc_sel_freq = ADRIA.decision.selection_frequency(rs.ranks, :seed)

        loc_sel_score = ADRIA.decision.selection_score(rs.ranks, :seed)
        @test all(0.0 .< loc_sel_score .< 1.0) || "Selection score values out of bounds"

        loc_sel_score = ADRIA.decision.selection_score(rs.ranks, :seed; keep_time=false)
        @test all(0.0 .< loc_sel_score .< 1.0) || "Selection score values out of bounds"
    end
end
