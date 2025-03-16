using Test
using ADRIA.Distributions
using ADRIA.decision.JMcDM
using ADRIA.decision:
    subtypes,
    SeedPreferences,
    cluster_diversity,
    geographic_separation,
    decision_matrix,
    select_locations

@testset "Location selection" begin
    """
    Tests to ensure cluster diversity and geographic separation of intervention sites
    are applied correctly.
    """

    @testset "Select under-represented cluster" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
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
        separation_scores = geographic_separation(loc_data.mean_to_neighbor)

        # Test that underrepresented clusters get higher scores
        @test diversity_scores[n_locs - 1] > diversity_scores[1] ||
            "Underrepresented cluster was not assigned a higher score"
        @test diversity_scores[n_locs - 2] > diversity_scores[1] ||
            "Underrepresented cluster was not assigned a higher score"
        @test diversity_scores[n_c1 + 2] > diversity_scores[n_c1 - 1] ||
            "Underrepresented cluster was not assigned a higher score"

        method = first(ADRIA.mcda_methods())

        ADRIA.fix_factor!(dom; seed_cluster_diversity=1.0)
        # ADRIA.fix_factor!(dom, seed_geographic_separation=1.0)

        sp = ADRIA.decision.SeedPreferences(dom)
        dm = decision_matrix(
            dom.loc_ids,
            sp.names;
            seed_depth=loc_data.depth_med,
            seed_in_connectivity=zeros(n_locs),
            seed_out_connectivity=zeros(n_locs),
            seed_heat_stress=zeros(n_locs),
            seed_wave_stress=zeros(n_locs),
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

    @testset "Select location closest to other locations" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
        loc_data = dom.loc_data
        n_locs = length(dom.loc_ids)
        min_locs = 10  # select at least 10 locations

        # Make selecting a location that is close to other locations more important
        ADRIA.fix_factor!(dom; seed_geographic_separation=1.0)

        diversity_scores = cluster_diversity(loc_data.cluster_id)
        separation_scores = geographic_separation(loc_data.mean_to_neighbor)

        closest_location = argmin(separation_scores)

        method = first(ADRIA.mcda_methods())

        sp = SeedPreferences(dom)
        dm = decision_matrix(
            dom.loc_ids,
            sp.names;
            seed_depth=loc_data.depth_med,
            seed_in_connectivity=zeros(n_locs),
            seed_out_connectivity=zeros(n_locs),
            seed_heat_stress=zeros(n_locs),
            seed_wave_stress=zeros(n_locs),
            seed_coral_cover=Float64.(rand(1:100, n_locs)),
            seed_cluster_diversity=diversity_scores,
            seed_geographic_separation=separation_scores
        )

        # If geographic separation is preferred, then should select locations that are
        # close to other locations
        result = select_locations(
            sp,
            dm,
            method,
            collect(1:n_locs),
            min_locs
        )

        msg = "Location with lowest mean distance to neighbors was not selected!"
        @test "reef_$(lpad(closest_location, 2, "0"))" ∈ result || msg
    end
end
