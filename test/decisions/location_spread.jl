using Test
using Random
using Statistics
using ADRIA.Distributions
using ADRIA.decision.JMcDM
using ADRIA: centroids, farthest_point_ordering
using ADRIA.decision:
    subtypes,
    SeedPreferences,
    cluster_diversity,
    geographic_separation,
    decision_matrix,
    select_locations

"""
Minimum pairwise distance within a set of locations, used as a measure of how spread
out a selection is.
"""
function _min_pairwise_distance(coords::Vector{Tuple{Float64,Float64}})::Float64
    n = length(coords)
    n < 2 && return Inf

    return minimum(
        ADRIA.Distances.haversine(coords[i], coords[j])
        for i in 1:n for j in (i + 1):n
    )
end

@testset "Location selection" begin
    """
    Tests to ensure cluster diversity and geographic separation of intervention sites
    are applied correctly.
    """

    @testset "Farthest point ordering" begin
        # Regular 5x5 grid, small enough that haversine distances are near-planar
        grid = vec([(Float64(x), Float64(y)) for x in 0:4, y in 0:4])
        order = farthest_point_ordering(grid)

        # Ordering is a permutation of all locations
        @test sort(order) == collect(1:25)

        # First selection is the location closest to the centroid (2, 2)
        @test grid[order[1]] == (2.0, 2.0)

        # Next four selections are the four corners, in some order
        corners = Set([(0.0, 0.0), (0.0, 4.0), (4.0, 0.0), (4.0, 4.0)])
        @test Set(grid[order[2:5]]) == corners

        # Each prefix is at least as spread out as the next one
        spreads = [_min_pairwise_distance(grid[order[1:k]]) for k in 2:25]
        @test issorted(spreads; rev=true)

        # Edge cases
        @test farthest_point_ordering(Tuple{Float64,Float64}[]) == Int64[]
        @test farthest_point_ordering([(1.0, 1.0)]) == [1]
    end

    @testset "Geographic separation scores" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        loc_coords = centroids(dom.loc_data)
        n_locs = length(loc_coords)

        separation_scores = geographic_separation(loc_coords)

        # Scores are a permutation of a decay logarithmic in the rank, spanning [0, 1]
        @test sort(separation_scores) ≈
            sort(1.0 .- log.(1:n_locs) ./ log(n_locs))
        @test maximum(separation_scores) ≈ 1.0
        @test minimum(separation_scores) ≈ 0.0

        # The decay steepens towards the top of the ranking, so the criteria keeps its
        # influence around the selection boundary whatever the number of locations
        # deployed to
        ranked_scores = sort(separation_scores; rev=true)
        @test issorted(diff(ranked_scores))

        # Highest scoring location is the one closest to the centroid
        centroid = (mean(first.(loc_coords)), mean(last.(loc_coords)))
        closest_to_centroid = argmin([
            ADRIA.Distances.haversine(coord, centroid) for coord in loc_coords
        ])
        @test argmax(separation_scores) == closest_to_centroid

        # Any prefix of the ranking is more spread out than a random subset of the
        # same size (a prefix holding every location is trivially identical, so is
        # not informative here)
        Random.seed!(101)
        ranked_idx = sortperm(separation_scores; rev=true)
        for k in filter(<(n_locs), [3, 5, 10])
            top_k_spread = _min_pairwise_distance(loc_coords[ranked_idx[1:k]])
            random_spreads = [
                _min_pairwise_distance(loc_coords[randperm(n_locs)[1:k]])
                for _ in 1:100
            ]

            @test top_k_spread > median(random_spreads)
        end
    end

    @testset "Select under-represented cluster" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        loc_data = dom.loc_data
        n_locs = length(dom.loc_ids)
        min_locs = 10  # select at least 10 locations

        # Setup with clearly imbalanced clusters
        c1 = fill(1, ceil(Int64, n_locs * 0.6))
        n_c1 = length(c1)
        c2 = fill(2, ceil(Int64, n_locs * 0.2))
        c3 = fill(3, ceil(Int64, n_locs * 0.1) + 1)
        cluster_ids = vcat(c1, c2, c3)[1:n_locs]  # truncate to n_locs

        diversity_scores = cluster_diversity(cluster_ids)
        separation_scores = geographic_separation(centroids(loc_data))

        # Test that underrepresented clusters get higher scores
        @test diversity_scores[n_locs - 1] > diversity_scores[1] ||
            "Underrepresented cluster was not assigned a higher score"
        @test diversity_scores[n_locs - 2] > diversity_scores[1] ||
            "Underrepresented cluster was not assigned a higher score"
        @test diversity_scores[n_c1 + 2] > diversity_scores[n_c1 - 1] ||
            "Underrepresented cluster was not assigned a higher score"

        method = first(ADRIA.mcda_methods())

        ADRIA.fix_factor!(dom; seed_cluster_diversity=1.0)

        sp = ADRIA.decision.SeedPreferences(dom)
        dm = decision_matrix(
            dom.loc_ids,
            sp.names;
            seed_depth=loc_data.depth_med,
            seed_in_connectivity=zeros(n_locs),
            seed_out_connectivity=zeros(n_locs),
            seed_heat_stress=zeros(n_locs),
            seed_coral_cover=Float64.(rand(1:100, n_locs)),
            seed_cluster_diversity=diversity_scores,
            seed_geographic_separation=separation_scores
        )

        # When only cluster diversity matters, should select from underrepresented clusters
        result = select_locations(
            sp,
            dm,
            method,
            collect(1:n_locs),
            min_locs
        )

        # Should include locations from clusters 2 and 3
        @test "reef_$(lpad(n_c1+2, 2, "0"))" ∈ result ||
            "Under-represented cluster 2 not selected!"
        @test "reef_$(lpad(n_locs-2, 2, "0"))" ∈ result ||
            "Under-represented cluster 3 not selected!"
    end

    @testset "Select geographically spread out locations" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        loc_data = dom.loc_data
        loc_coords = centroids(loc_data)
        n_locs = length(dom.loc_ids)
        min_locs = 5  # select at least 5 locations

        # Make geographic coverage the only criteria that matters
        ADRIA.fix_factor!(dom; seed_geographic_separation=1.0)

        diversity_scores = cluster_diversity(loc_data.cluster_id)
        separation_scores = geographic_separation(loc_coords)

        method = first(ADRIA.mcda_methods())

        sp = SeedPreferences(dom)
        dm = decision_matrix(
            dom.loc_ids,
            sp.names;
            seed_depth=loc_data.depth_med,
            seed_in_connectivity=zeros(n_locs),
            seed_out_connectivity=zeros(n_locs),
            seed_heat_stress=zeros(n_locs),
            seed_coral_cover=Float64.(rand(1:100, n_locs)),
            seed_cluster_diversity=diversity_scores,
            seed_geographic_separation=separation_scores
        )

        result = select_locations(
            sp,
            dm,
            method,
            collect(1:n_locs),
            min_locs
        )

        # Selected locations should be more spread out than a random selection of the
        # same size
        selected_idx = findall(in(result), dom.loc_ids)
        selected_spread = _min_pairwise_distance(loc_coords[selected_idx])

        Random.seed!(101)
        random_spreads = [
            _min_pairwise_distance(loc_coords[randperm(n_locs)[1:length(selected_idx)]])
            for _ in 1:100
        ]

        @test selected_spread > median(random_spreads)
    end
end
