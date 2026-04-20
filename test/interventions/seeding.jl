using Test

using ADRIA
using ADRIA.Distributions
using ADRIA: distribute_seeded_corals, location_k, update_tolerance_distribution!
using ADRIA: At

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

# @testset "Seeding" begin
#     n_groups = 5
#     n_sizes = 7
#     seeding_devices_per_m2::Float64 = 9.0

#     # extract inputs for function
#     total_loc_area = loc_area(ADRIA_DOM_45)
#     k = location_k(ADRIA_DOM_45)
#     current_cover = zeros(size(total_loc_area))

#     # calculate available space
#     available_space = vec((total_loc_area .* k) .- current_cover)

#     # Area of a single 1-year-old colony (diameter 3.5 cm), same for all groups
#     colony_area_m² = pi * (3.5 / 2)^2

#     # Randomly generate total seeded area per functional group
#     seeded_area = ADRIA.DataCube(
#         rand(Uniform(0.0, 500.0), n_groups);
#         taxa=["N_seed_TA", "N_seed_CA", "N_seed_CNA", "N_seed_SM", "N_seed_LM"]
#     )
#     # Approximate number of corals deployed from the seeded_area
#     seeded_volume = seeded_area.data ./ colony_area_m²
#     colony_areas = fill(colony_area_m², n_groups)
#     @testset "Check coral seed distribution ($i)" for i in 1:10
#         seed_locs = rand(1:length(total_loc_area), 5)

#         # evaluate seeding distributions
#         seed_dist, _ = distribute_seeded_corals(
#             total_loc_area[seed_locs],
#             available_space[seed_locs],
#             seeded_volume,
#             colony_areas,
#             seeding_devices_per_m2
#         )

#         # Area to be seeded for each site
#         total_area_seed = seed_dist .* total_loc_area[seed_locs]'

#         # total area of seeded corals
#         total_area_coral_out = sum(total_area_seed; dims=2)

#         # absolute available area to seed for selected sites
#         selected_avail_space = available_space[seed_locs]

#         # Index of max proportion for available space
#         # The selected location (row number) should match
#         abs_seed_area = seed_dist' .* total_loc_area[seed_locs]
#         max_ind_out = argmax(abs_seed_area)[1]
#         max_ind = argmax(selected_avail_space)

#         # Index of min proportion for available space
#         # As `abs_seed_area` is a matrix, Cartesian indices are returned
#         # when finding the mininum/maximum argument (`argmin()`).
#         # The selected location (row number) should match.
#         min_ind_out = argmin(abs_seed_area)[1]
#         min_ind = argmin(selected_avail_space)

#         # Distributed areas and seeded areas
#         # total_area_coral_out has integer taxa indices (1:n_groups) because
#         # distribute_seeded_corals builds its axis from plain Vector inputs
#         area_TA = total_area_coral_out[taxa=At(1), locations=1][1]
#         seed_TA = seeded_area[taxa=At("N_seed_TA")][1]

#         area_CA = total_area_coral_out[taxa=At(2), locations=1][1]
#         seed_CA = seeded_area[taxa=At("N_seed_CA")][1]

#         area_SM = total_area_coral_out[taxa=At(4), locations=1][1]
#         seed_SM = seeded_area[taxa=At("N_seed_SM")][1]

#         approx_zero(x) = abs(x) + one(1.0) ≈ one(1.0)
#         @test approx_zero(seed_TA - area_TA) && approx_zero(seed_CA - area_CA) &&
#               approx_zero(seed_CA - area_CA) ||
#             "Area of corals seeded not equal to (colony area) * (number or corals)"
#         @test all(seed_dist .< 1.0) || "Some proportions of seeded corals greater than 1"
#         @test all(seed_dist .>= 0.0) || "Some proportions of seeded corals less than zero"
#         @test all(total_area_seed .< selected_avail_space') ||
#             "Area seeded greater than available area"
#         @test (max_ind_out == max_ind) ||
#             "Maximum distributed proportion of seeded coral not seeded in largest available area."
#         @test (min_ind_out == min_ind) ||
#             "Minimum distributed proportion of seeded coral not seeded in smallest available area."
#     end

#     @testset "DHW distribution priors" begin
#         n_locs = 10
#         C_cover_t = rand(n_groups, n_sizes, n_locs)  # size class, locations
#         a_adapt = rand(2.0:6.0, n_groups)
#         total_location_area = fill(5000.0, n_locs)

#         seed_locs = rand(1:n_locs, 5)  # Pick 5 random locations

#         leftover_space_m² = fill(500.0, n_locs)
#         seed_sc = ADRIA.seed_size_groups(n_groups, n_sizes)

#         # Initial distributions
#         d = truncated(Normal(1.0, 0.15), 0.0, 3.0)
#         c_dist_t = rand(d, 5, n_sizes, n_locs)
#         orig_dist = copy(c_dist_t)

#         dist_std = rand(n_groups, n_sizes)

#         # Absolute number of corals seeded is not required
#         proportional_increase, _ = distribute_seeded_corals(
#             total_location_area[seed_locs],
#             leftover_space_m²[seed_locs],
#             seeded_volume,
#             colony_areas,
#             seeding_devices_per_m2
#         )

#         update_tolerance_distribution!(
#             proportional_increase,
#             C_cover_t,
#             c_dist_t,
#             c_dist_t[:, 1, :],
#             dist_std,
#             seed_locs,
#             seed_sc,
#             a_adapt
#         )

#         # Ensure correct priors/weightings for each location
#         for loc in seed_locs
#             for (i, sc) in enumerate(findall(seed_sc))
#                 @test c_dist_t[sc, loc] > orig_dist[sc, loc] ||
#                     "Expected mean of distribution to shift | SC: $sc ; Location: $loc"
#             end
#         end
#     end
# end

@testset "Seed log matches target location sets" begin
    dom = deepcopy(ADRIA_DOM_45)
    rand_seed_locs = rand(dom.loc_ids, 70, 2)

    locs_1 = rand_seed_locs[:, 1]
    locs_2 = rand_seed_locs[:, 2]

    weight_1 = 0.4
    weight_2 = 0.6

    ADRIA.set_seed_target_locations!(
        dom, [(weight=weight_1, target_locs=locs_1), (weight=weight_2, target_locs=locs_2)]
    )

    ADRIA.fix_factor!(
        dom;
        N_seed_TA=500_000.0,
        N_seed_CA=500_000.0,
        N_seed_CNA=500_000.0,
        N_seed_SM=500_000.0,
        N_seed_LM=500_000.0,
        seed_year_start=1.0,
        seed_years=75.0,
        seed_deployment_freq=1.0,
        seed_strategy=1.0
    )

    num_samples = 4
    scens = ADRIA.sample_guided(dom, num_samples)
    rs = ADRIA.run_scenarios(dom, scens, "45")

    # Total seeds requested per scenario (sum across all coral species)
    seed_cols = names(scens, contains.(names(scens), "N_seed"))
    N_seed = vec(sum(Matrix(scens[:, seed_cols]); dims=2))

    excl_locs_1 = setdiff(locs_1, locs_2)
    seed_log_exc1 = if !isempty(excl_locs_1)
        dropdims(
            sum(
                rs.seed_log[locations=dom.loc_ids .∈ [excl_locs_1]];
                dims=(:coral_id, :locations)
            );
            dims=(:coral_id, :locations)
        )
    else
        []
    end

    excl_locs_2 = setdiff(locs_2, locs_1)
    seed_log_exc2 = if !isempty(excl_locs_2)
        dropdims(
            sum(
                rs.seed_log[locations=dom.loc_ids .∈ [excl_locs_2]];
                dims=(:coral_id, :locations)
            );
            dims=(:coral_id, :locations)
        )
    else
        []
    end

    both_locs = intersect(locs_1, locs_2)
    seed_log_both = dropdims(
        sum(
            rs.seed_log[locations=dom.loc_ids .∈ [both_locs]]; dims=(:coral_id, :locations)
        );
        dims=(:coral_id, :locations)
    )

    for s in 1:num_samples
        seed_log = seed_log_both[scenarios=At(s)]
        @test all(seed_log .<= N_seed[s] .|| seed_log .≈ N_seed[s]) ||
            "Scenario $s: combined seed log exceeds N_seed"

        # For locations exclusive to set 1 (never targeted by set 2), logged seeds
        # must not exceed weight_1 * N_seed — their full allocated share
        if !isempty(excl_locs_1)
            excl_log_1 = seed_log_exc1[scenarios=At(s)]
            budget_1 = weight_1 * N_seed[s]
            @test all(excl_log_1 .<= budget_1 .|| excl_log_1 .≈ budget_1) ||
                "Scenario $s: seed log for locations exclusive to set 1 exceeds weight_1 * N_seed"
        end

        # For locations exclusive to set 2 (never targeted by set 1), logged seeds
        # must not exceed weight_2 * N_seed — their full allocated share
        if !isempty(excl_locs_2)
            excl_log_2 = seed_log_exc2[scenarios=At(s)]
            budget_2 = weight_2 * N_seed[s]
            @test all(excl_log_2 .<= budget_2 .|| excl_log_2 .≈ budget_2) ||
                "Scenario $s: seed log for locations exclusive to set 2 exceeds weight_2 * N_seed"
        end
    end
end
