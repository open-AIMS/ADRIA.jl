using Test

using ADRIA
using ADRIA.Distributions
using ADRIA: distribute_seeded_corals, location_k, update_tolerance_distribution!
using ADRIA: vary_locations, vary_n_corals, vary_seed_density, seed_cap_density
using ADRIA: At

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end
@testset "Seeding" begin

    # dummy density, using vary density strategy will recalculate density
    DUMMY_DENSITY = 5.0

    # extract inputs for function
    total_loc_area = loc_area(ADRIA_DOM_45)
    k = location_k(ADRIA_DOM_45)
    current_cover = zeros(size(total_loc_area))

    # calculate available space
    available_space = vec((total_loc_area .* k) .- current_cover)

    # Randomly generate seeded area
    seeded_area = ADRIA.DataCube(
        rand(Uniform(0.0, 500.0), 3); taxa=["N_seed_TA", "N_seed_CA", "N_seed_SM"]
    )
    # Aprroximate number of corals deployed from the seeded_area
    seeded_volume = seeded_area.data ./ (pi * (3.5 / 2)^2)
    @testset "Check coral seed distribution ($i)" for i in 1:10
        seed_locs = rand(1:length(total_loc_area), 5)

        # evaluate seeding distributions
        seed_dist, _ = distribute_seeded_corals(
            :VARY_SEED_DENSITY,
            seed_locs,
            total_loc_area,
            available_space,
            seeded_area,
            seeded_volume,
            DUMMY_DENSITY
        )

        # Area to be seeded for each site
        total_area_seed = seed_dist .* total_loc_area[seed_locs]'

        # total area of seeded corals
        total_area_coral_out = sum(total_area_seed; dims=2)

        # absolute available area to seed for selected sites
        selected_avail_space = available_space[seed_locs]

        # Index of max proportion for available space
        # The selected location (row number) should match
        abs_seed_area = seed_dist' .* total_loc_area[seed_locs]
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
        @test approx_zero(seed_TA - area_TA) && approx_zero(seed_CA - area_CA) &&
              approx_zero(seed_CA - area_CA) ||
            "Area of corals seeded not equal to (colony area) * (number or corals)"
        @test all(seed_dist .< 1.0) || "Some proportions of seeded corals greater than 1"
        @test all(seed_dist .>= 0.0) || "Some proportions of seeded corals less than zero"
        @test all(total_area_seed .< selected_avail_space') ||
            "Area seeded greater than available area"
        @test (max_ind_out == max_ind) ||
            "Maximum distributed proportion of seeded coral not seeded in largest available area."
        @test (min_ind_out == min_ind) ||
            "Minimum distributed proportion of seeded coral not seeded in smallest available area."
    end

    @testset "DHW distribution priors" begin
        C_cover_t = rand(5, 7, 10)  # size class, locations
        a_adapt = rand(2.0:6.0, 5, 7)
        total_location_area = fill(5000.0, 10)

        seed_locs = rand(1:10, 5)  # Pick 5 random locations

        leftover_space_m² = fill(500.0, 10)
        seed_sc = falses(5, 7)
        seed_sc[[2, 3, 5], 1] .= true

        # Initial distributions
        d = truncated(Normal(1.0, 0.15), 0.0, 3.0)
        c_dist_t = rand(d, 5, 7, 10)
        orig_dist = copy(c_dist_t)

        dist_std = rand(5, 7)

        # Absolute number of corals seeded is not required
        proportional_increase, _ = distribute_seeded_corals(
            :VARY_SEED_DENSITY,
            seed_locs,
            total_location_area,
            leftover_space_m²,
            seeded_area,
            seeded_volume,
            DUMMY_DENSITY
        )

        update_tolerance_distribution!(
            proportional_increase,
            C_cover_t,
            c_dist_t,
            dist_std,
            seed_locs,
            seed_sc,
            a_adapt
        )

        # Ensure correct priors/weightings for each location
        for loc in seed_locs
            for (i, sc) in enumerate(findall(seed_sc))
                @test c_dist_t[sc, loc] > orig_dist[sc, loc] ||
                    "Expected mean of distribution to shift | SC: $sc ; Location: $loc"
            end
        end
    end

    @testset "Seeding Strategies" begin
        @testset "Vary Locations" begin
            target_density = 5.0
            available_space = fill(10.0, 5)
            cumulative_space = cumsum(available_space)
            # deploy [49, 98, 147, 196, ...] corals
            base_deployment = 49
            test_in_corals = [
                fill(base_deployment / 3, 3) .* i for i in 1:length(available_space)
            ]
            # The expected deployment density should be [49, 98, 147, ...] ./ [50, 100, ...]
            expected_density = [
                base_deployment * i / sp for (i, sp) in enumerate(cumulative_space)
            ]

            for (exp, inp) in zip(expected_density, test_in_corals)
                new_density, n_corals, n_iv_locs = vary_locations(
                    available_space, target_density, inp, 5
                )
                @test all(n_corals .== inp)
                @test all([new_density, n_iv_locs] .≈ [exp, n_iv_locs])
            end
            # similar test case as previous but with different density
            target_density = 10.0
            available_space = fill(30.0, 4)
            cumulative_space = cumsum(available_space)
            base_deployment = 299
            test_in_corals = [
                fill(base_deployment / 3, 3) .* i for i in 1:length(available_space)
            ]
            expected_density = [
                base_deployment * i / sp for (i, sp) in enumerate(cumulative_space)
            ]

            for (exp, inp) in zip(expected_density, test_in_corals)
                new_density, n_corals, n_iv_locs = vary_locations(
                    available_space, target_density, inp, 5
                )
                @test all(n_corals .== inp)
                @test all([new_density, n_iv_locs] .== [exp, n_iv_locs])
            end

            target_density = 1.0
            available_space = [
                1.0, 149.0
            ]
            n_corals = [80.0, 80.0, 80.0]
            expected_density = 240 / sum(available_space)
            out_density, out_corals, out_iv_locs = vary_locations(
                available_space, target_density, n_corals, 2
            )
            @test expected_density == out_density
            @test out_corals == n_corals
            @test 2 == out_iv_locs
        end

        @testset "Vary # Corals" begin
            target_density = [3.0, 5.0, 7.0, 10.0, 12.0]
            available_space = [
                10.0, 20.0, 30.0, 40.0, 50.0
            ]
            n_corals = [50.0, 50.0, 50.0]
            n_iv_locs = length(available_space)

            # Test density is held constant for a range of target densities
            # Test new number of corals satisfies available space
            for t_den in target_density
                new_den, new_n_corals, n_iv_locs = vary_n_corals(
                    available_space, t_den, n_corals, n_iv_locs
                )

                @test all(new_den == t_den)
                @test sum(new_n_corals) / new_den <= sum(available_space)
            end

            target_density = 5.0
            available_space = rand(Uniform(10.0, 60.0), 20)
            n_iv_locs = 5:5:20

            # Test new number of corals satisfies available space for variety of n_sites and available space
            for n_iv in n_iv_locs
                new_den, new_n_corals, n_iv_locs = vary_n_corals(
                    available_space[1:n_iv], target_density, n_corals, n_iv
                )

                @test all(new_den == target_density)
                @test sum(new_n_corals) / new_den <= sum(available_space)
            end
        end

        @testset "Vary Seed Density" begin
            # Target density is unused
            target_density = 5.0
            available_space = [
                10.0, 20.0, 30.0, 40.0, 50.0
            ]
            n_corals = [50.0, 50.0, 50.0]
            expected_densities = 150.0 ./ [
                10.0, 30.0, 60.0, 100.0, 150.0
            ]

            for (exp, n_iv) in zip(expected_densities, 1:5)
                dens, n_c, n_l = vary_seed_density(
                    available_space, target_density, n_corals, n_iv
                )
                @test dens == exp
                @test all(n_c .== n_corals)
                @test n_iv == n_l
            end
        end

        @testset "Cap Density" begin
            target_density = [3.0, 5.0, 7.0, 10.0, 12.0]
            available_space = [
                10.0, 20.0, 30.0, 40.0, 50.0
            ]
            n_corals = [50.0, 50.0, 50.0]

            # For a range of target densities, check updated density is <= target and new density satisfies avaialble space
            for t_den in target_density
                new_den, new_n_corals, n_iv_locs = ADRIA.seed_cap_density(
                    available_space, t_den, n_corals, 5
                )

                @test (new_den <= t_den)
                @test sum(new_n_corals) / new_den <= sum(available_space)
            end
        end
    end
end
