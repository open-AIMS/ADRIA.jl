using Test

using ADRIA
using ADRIA.Distributions
using ADRIA: distribute_seeded_corals, location_k, seed_corals!

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

@testset "Seeding" begin
    # first test function on example domain
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)

    # extract inputs for function
    total_site_area = site_area(dom)
    k = location_k(dom)
    current_cover = zeros(size(total_site_area))

    # calculate available space
    available_space = vec((total_site_area .* k) .- current_cover)

    # Randomly generate seeded area
    seeded_area = ADRIA.DataCube(
        rand(Uniform(0.0, 500.0), 3); taxa=["N_seed_TA", "N_seed_CA", "N_seed_SM"]
    )

    @testset "Check coral seed distribution ($i)" for i in 1:10
        seed_locs = rand(1:length(total_site_area), 5)

        # evaluate seeding distributions
        seed_dist = distribute_seeded_corals(
            total_site_area[seed_locs], available_space[seed_locs], seeded_area
        )

        # Area to be seeded for each site
        total_area_seed = seed_dist .* total_site_area[seed_locs]'

        # total area of seeded corals
        total_area_coral_out = sum(total_area_seed, dims=2)

        # absolute available area to seed for selected sites
        selected_avail_space = available_space[seed_locs]

        # Index of max proportion for available space
        # The selected location (row number) should match
        abs_seed_area = seed_dist' .* total_site_area[seed_locs]
        max_ind_out = argmax(abs_seed_area)[1]
        max_ind = argmax(selected_avail_space)

        # Index of min proportion for available space
        # As `abs_seed_area` is a matrix, Cartesian indices are returned
        # when finding the mininum/maximum argument (`argmin()`).
        # The selected location (row number) should match.
        min_ind_out = argmin(abs_seed_area)[1]
        min_ind = argmin(selected_avail_space)

        # Distributed areas and seeded areas
        area_TA = total_area_coral_out[taxa=At("N_seed_TA"), locations=1][1]
        seed_TA = seeded_area[taxa=At("N_seed_TA")][1]

        area_CA = total_area_coral_out[taxa=At("N_seed_CA"), locations=1][1]
        seed_CA = seeded_area[taxa=At("N_seed_CA")][1]

        area_SM = total_area_coral_out[taxa=At("N_seed_SM"), locations=1][1]
        seed_SM = seeded_area[taxa=At("N_seed_SM")][1]

        approx_zero(x) = abs(x) + one(1.0) ≈ one(1.0)
        @test approx_zero(seed_TA - area_TA) && approx_zero(seed_CA - area_CA) && approx_zero(seed_CA - area_CA) || "Area of corals seeded not equal to (colony area) * (number or corals)"
        @test all(seed_dist .< 1.0) || "Some proportions of seeded corals greater than 1"
        @test all(seed_dist .>= 0.0) || "Some proportions of seeded corals less than zero"
        @test all(total_area_seed .< selected_avail_space') || "Area seeded greater than available area"
        @test (max_ind_out == max_ind) || "Maximum distributed proportion of seeded coral not seeded in largest available area."
        @test (min_ind_out == min_ind) || "Minimum distributed proportion of seeded coral not seeded in smallest available area."
    end

    @testset "DHW distribution priors" begin
        C_cover_t = rand(36, 10)  # size class, locations
        a_adapt = rand(2.0:6.0, 36)
        total_location_area = fill(5000.0, 10)

        seed_locs = rand(1:10, 5)  # Pick 5 random locations

        leftover_space_m² = fill(500.0, 10)

        Yseed = zeros(2, 3, 10)
        seed_sc = BitVector([i ∈ [2, 8, 15] for i in 1:36])

        # Initial distributions
        d = truncated(Normal(1.0, 0.15), 0.0, 3.0)
        c_dist_t = rand(d, 36, 10)
        orig_dist = copy(c_dist_t)

        dist_std = rand(36)
        seed_corals!(C_cover_t, total_location_area, leftover_space_m², seed_locs, seeded_area, seed_sc,
            a_adapt, @view(Yseed[1, :, :]), dist_std, c_dist_t)

        # Ensure correct priors/weightings for each location
        for loc in seed_locs
            for (i, sc) in enumerate(findall(seed_sc))
                prior1 = Yseed[1, i, loc] ./ C_cover_t[sc, loc]
                expected = [prior1, 1.0 - prior1]
                @test c_dist_t[sc, loc] > orig_dist[sc, loc] || "Expected mean of distribution to shift | SC: $sc ; Location: $loc"
            end
        end
    end
end
