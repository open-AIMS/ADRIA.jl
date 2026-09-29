using Test

using ADRIA
using ADRIA.Distributions
using ADRIA: distribute_CAq_corals, location_k, update_tolerance_distribution!
using ADRIA: At

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

# @testset "Coral aquaculture" begin
#     n_groups = 5
#     n_sizes = 7
#     iv_CAq_devices_per_m2::Float64 = 9.0

#     # extract inputs for function
#     total_loc_area = loc_area(ADRIA_DOM_45)
#     k = location_k(ADRIA_DOM_45)
#     current_cover = zeros(size(total_loc_area))

#     # calculate available space
#     available_space = vec((total_loc_area .* k) .- current_cover)

#     # Area of a single 1-year-old colony (diameter 3.5 cm), same for all groups
#     colony_area_m² = pi * (3.5 / 2)^2

#     # Randomly generate total coral aquaculture area per functional group
#     CAq_area = ADRIA.DataCube(
#         rand(Uniform(0.0, 500.0), n_groups);
#         taxa=["iv_CAq_N_TA", "iv_CAq_N_CA", "iv_CAq_N_CNA", "iv_CAq_N_SM", "iv_CAq_N_LM"]
#     )
#     # Approximate number of corals deployed from the CAq_area
#     CAq_volume = CAq_area.data ./ colony_area_m²
#     colony_areas = fill(colony_area_m², n_groups)
#     @testset "Check coral aquaculture distribution ($i)" for i in 1:10
#         CAq_locs = rand(1:length(total_loc_area), 5)

#         # evaluate coral aquaculture distributions
#         CAq_dist, _ = distribute_CAq_corals(
#             total_loc_area[CAq_locs],
#             available_space[CAq_locs],
#             CAq_volume,
#             colony_areas,
#             iv_CAq_devices_per_m2
#         )

#         # Area to be deployed via coral aquaculture for each site
#         total_area_CAq = CAq_dist .* total_loc_area[CAq_locs]'

#         # total area of coral aquaculture corals
#         total_area_coral_out = sum(total_area_CAq; dims=2)

#         # absolute available area for coral aquaculture at selected sites
#         selected_avail_space = available_space[CAq_locs]

#         # Index of max proportion for available space
#         # The selected location (row number) should match
#         abs_CAq_area = CAq_dist' .* total_loc_area[CAq_locs]
#         max_ind_out = argmax(abs_CAq_area)[1]
#         max_ind = argmax(selected_avail_space)

#         # Index of min proportion for available space
#         # As `abs_CAq_area` is a matrix, Cartesian indices are returned
#         # when finding the mininum/maximum argument (`argmin()`).
#         # The selected location (row number) should match.
#         min_ind_out = argmin(abs_CAq_area)[1]
#         min_ind = argmin(selected_avail_space)

#         # Distributed areas and coral aquaculture areas
#         # total_area_coral_out has integer taxa indices (1:n_groups) because
#         # distribute_CAq_corals builds its axis from plain Vector inputs
#         area_TA = total_area_coral_out[taxa=At(1), locations=1][1]
#         CAq_TA = CAq_area[taxa=At("iv_CAq_N_TA")][1]

#         area_CA = total_area_coral_out[taxa=At(2), locations=1][1]
#         CAq_CA = CAq_area[taxa=At("iv_CAq_N_CA")][1]

#         area_SM = total_area_coral_out[taxa=At(4), locations=1][1]
#         CAq_SM = CAq_area[taxa=At("iv_CAq_N_SM")][1]

#         approx_zero(x) = abs(x) + one(1.0) ≈ one(1.0)
#         @test approx_zero(CAq_TA - area_TA) && approx_zero(CAq_CA - area_CA) &&
#               approx_zero(CAq_CA - area_CA) ||
#             "Area of coral aquaculture corals not equal to (colony area) * (number or corals)"
#         @test all(CAq_dist .< 1.0) || "Some proportions of coral aquaculture corals greater than 1"
#         @test all(CAq_dist .>= 0.0) || "Some proportions of coral aquaculture corals less than zero"
#         @test all(total_area_CAq .< selected_avail_space') ||
#             "Area deployed via coral aquaculture greater than available area"
#         @test (max_ind_out == max_ind) ||
#             "Maximum distributed proportion of coral aquaculture corals not deployed in largest available area."
#         @test (min_ind_out == min_ind) ||
#             "Minimum distributed proportion of coral aquaculture corals not deployed in smallest available area."
#     end

#     @testset "DHW distribution priors" begin
#         n_locs = 10
#         C_cover_t = rand(n_groups, n_sizes, n_locs)  # size class, locations
#         iv_CAq_a_adapt = rand(2.0:6.0, n_groups)
#         total_location_area = fill(5000.0, n_locs)

#         CAq_locs = rand(1:n_locs, 5)  # Pick 5 random locations

#         leftover_space_m² = fill(500.0, n_locs)
#         CAq_sc = ADRIA.CAq_size_groups(n_groups, n_sizes)

#         # Initial distributions
#         d = truncated(Normal(1.0, 0.15), 0.0, 3.0)
#         c_dist_t = rand(d, 5, n_sizes, n_locs)
#         orig_dist = copy(c_dist_t)

#         dist_std = rand(n_groups, n_sizes)

#         # Absolute number of corals deployed via coral aquaculture is not required
#         proportional_increase, _ = distribute_CAq_corals(
#             total_location_area[CAq_locs],
#             leftover_space_m²[CAq_locs],
#             CAq_volume,
#             colony_areas,
#             iv_CAq_devices_per_m2
#         )

#         update_tolerance_distribution!(
#             proportional_increase,
#             C_cover_t,
#             c_dist_t,
#             c_dist_t[:, 1, :],
#             dist_std,
#             CAq_locs,
#             CAq_sc,
#             iv_CAq_a_adapt
#         )

#         # Ensure correct priors/weightings for each location
#         for loc in CAq_locs
#             for (i, sc) in enumerate(findall(CAq_sc))
#                 @test c_dist_t[sc, loc] > orig_dist[sc, loc] ||
#                     "Expected mean of distribution to shift | SC: $sc ; Location: $loc"
#             end
#         end
#     end
# end

@testset "Coral aquaculture log matches target location sets" begin
    dom = deepcopy(ADRIA_DOM_45)

    n_locs = length(dom.loc_ids)
    half = n_locs ÷ 2
    locs_1 = dom.loc_ids[1:half]
    locs_2 = dom.loc_ids[(half + 1):end]

    weight_1 = TEST_TARGET_WEIGHT_1
    weight_2 = TEST_TARGET_WEIGHT_2

    ADRIA.set_CAq_target_locations!(
        dom, [(weight=weight_1, target_locs=locs_1), (weight=weight_2, target_locs=locs_2)]
    )

    ADRIA.fix_factor!(
        dom;
        iv_CAq_N_TA=500_000.0,
        iv_CAq_N_CA=500_000.0,
        iv_CAq_N_CNA=500_000.0,
        iv_CAq_N_SM=500_000.0,
        iv_CAq_N_LM=500_000.0,
        iv_CAq_year_start=1.0,
        iv_CAq_years=75.0,
        iv_CAq_deployment_freq=1.0,
        iv_CAq_strategy=1.0
    )

    num_samples = 4
    scens = ADRIA.sample_guided(dom, num_samples)
    rs = ADRIA.run_scenarios(dom, scens, "45")

    # Total coral aquaculture requested per scenario (sum across all coral species)
    CAq_cols = names(scens, contains.(names(scens), "N_CAq"))
    N_CAq_total = vec(sum(Matrix(scens[:, CAq_cols]); dims=2))

    CAq_log_1 = dropdims(
        sum(
            rs.CAq_log[locations = dom.loc_ids .∈ [locs_1]]; dims=(:coral_id, :locations)
        );
        dims=(:coral_id, :locations)
    )
    CAq_log_2 = dropdims(
        sum(
            rs.CAq_log[locations = dom.loc_ids .∈ [locs_2]]; dims=(:coral_id, :locations)
        );
        dims=(:coral_id, :locations)
    )

    # CAq_log is persisted as Float32; use Float32 rtol regardless of in-memory eltype
    fp32_rtol = sqrt(eps(Float32))
    for s = 1:num_samples
        log_1 = CAq_log_1[scenarios = At(s)]
        log_2 = CAq_log_2[scenarios = At(s)]
        budget_1 = weight_1 * N_CAq_total[s]
        budget_2 = weight_2 * N_CAq_total[s]
        @test all(log_1 .<= budget_1 .|| isapprox.(log_1, budget_1; rtol=fp32_rtol))
        @test all(log_2 .<= budget_2 .|| isapprox.(log_2, budget_2; rtol=fp32_rtol))
    end
end
