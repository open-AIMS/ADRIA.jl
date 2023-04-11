using Test
using Distributions
using ADRIA: mcda_normalize, filter_decision_matrix, create_intervention_matrix, create_criteria_store, create_tolerances_store


@testset "Create decision matrix" begin

    # Dummy data to create decision matrix from
    n_locations = rand(5:30)
    area = rand(Uniform(100.0, 3000.0), 1, n_locations)
    centr_out = rand(1, n_locations)
    centr_in = rand(1, n_locations)

    sum_cover = rand(1, n_locations)
    max_cover = rand(1, n_locations)

    centr_out = sum_cover .* area .* centr_out
    centr_in = sum_cover .* area .* centr_in
    coral_space = (max_cover .- sum_cover) .* area
    coral_space[coral_space.<0] .= 0.0
    location_ids = collect(1:n_locations)
    args = (iv__heat_stress=rand(1, n_locations),
        iv__wave_stress=rand(1, n_locations), iv__zones=rand(Uniform(0.0, 2.0), 1, n_locations),
        iv__coral_space=coral_space, iv__in_connectivity=centr_in, iv__out_connectivity=centr_out,
        iv__priority=zeros(n_locations))
    criteria_store = create_criteria_store(location_ids, args)

    @test size(criteria_store, 1) == n_locations || "Some locations missing in decision matrix storage."
    @test size(criteria_store, 2) == length(args) || "Some criteria missing in decision matrix storage."
    @test all(criteria_store.locations .== location_ids) || "Location names do not match location ids."
    @test all(criteria_store.criteria .== keys(args)) || "Criteria names do not match input criteria names."

    param_set = NamedDimsArray(rand(Float64, 4), factors=["iv__heat_stress__seed_shade", "iv__wave_stress__seed_shade", "iv__coral_space__seed", "iv__in_connectivity__seed"])
    area_to_seed = 962.11  # area of seeded corals in m^2

    # define functions for tolerances
    f_coral_cover(param) = area_to_seed * param

    tolerances = (iv__coral_space=(>, f_coral_cover(rand())),
        iv__heat_stress=(>, 1 - rand()),
        iv__wave_stress=(>, 1 - rand()))

    tolerances = create_tolerances_store(tolerances)
    filtered = findall(tolerances(:iv__coral_space).(criteria_store(:iv__coral_space)) .& tolerances(:iv__heat_stress).(criteria_store(:iv__heat_stress)) .& tolerances(:iv__wave_stress).(criteria_store(:iv__wave_stress)))
    filtered_criteria_store = ADRIA.filter_decision_matrix(criteria_store, tolerances)

    @test all(filtered_criteria_store.locations .== filtered) || "Some locations which should have been filtered out were not filtered."

    SE, wse = create_intervention_matrix(filtered_criteria_store, param_set, "seed")

    @test size(SE, 2) == size(criteria_store, 2) - 3 || "More or less criteria included for intervention than required."
    @test size(SE, 1) == size(filtered_criteria_store, 1) || "Number of locations in filtered set and decision matrix do not match."

end

@testset "MCDA normalisation" begin
    # randomised weights
    w = vec(rand(Uniform(0, 1), 7))

    # randomised decision matrix
    A = zeros(5, 7)
    A[:, 1] = [1.0, 2.0, 3.0, 4.0, 5.0]
    A[:, 2:6] = rand(Uniform(0, 1), (5, 5))
    A[:, 7] = rand(Uniform(100, 1000), (5, 1))

    norm_A = mcda_normalize(A[:, 2:end])
    norm_w = mcda_normalize(w)

    @test all((sqrt.(sum(norm_A .^ 2, dims=1)) .- 1.0) .< 0.0001) || "Decision matrix normalization not giving column sums = 1."
    @test (sum(norm_w) - 1.0) <= 0.001 || "MCDA weights not summing to one."
end
