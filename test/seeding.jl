using ADRIA
using ADRIA: distribute_seeded_corals, location_k
using Test
using Distributions

@testset "Seeding distribution" begin
    # first test function on example domain
    dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)

    # extract inputs for function
    total_location_area = location_area(dom)
    k = location_k(dom)
    current_cover = zeros(size(total_location_area))

    # calculate available space
    available_space = vec((total_location_area .* k) .- current_cover)

    @testset "Test for Seeding Distribtion ($i)" for i in 1:10
        prefseedlocations = rand(1:length(total_location_area), 5)

        # Randomly generate seeded area
        tmp = rand(Uniform(0.0, 500.0), 2)
        seeded_area = (TA=tmp[1], CA=tmp[2])

        # evaluate seeding distributions
        seed_dist = distribute_seeded_corals(total_location_area, prefseedlocations,
            available_space, seeded_area)

        # proportions of coral
        total_area_coral_TA = seeded_area[1]
        total_area_coral_CA = seeded_area[2]

        # Area to be seeded for each location
        area_TA = seed_dist.TA .* total_location_area[prefseedlocations]
        area_CA = seed_dist.CA .* total_location_area[prefseedlocations]

        # total area of seeded corals
        total_area_coral_TA_out = sum(area_TA)
        total_area_coral_CA_out = sum(area_CA)

        # absolute available area to seed for selected locations
        selected_avail_space = available_space[prefseedlocations]

        abs_seed_area = seed_dist.TA .* total_location_area[prefseedlocations]

        # index of max proportion for available space
        max_ind_out = findfirst(abs_seed_area .== maximum(abs_seed_area))
        max_ind = findfirst(selected_avail_space .== maximum(selected_avail_space))

        # index of min proportion for available space
        min_ind_out = findfirst(abs_seed_area .== minimum(abs_seed_area))
        min_ind = findfirst(selected_avail_space .== minimum(selected_avail_space))

        @test ((total_area_coral_TA - total_area_coral_TA_out) < 10^-5) && ((total_area_coral_CA - total_area_coral_CA_out) < 10^-5) || "Area of corals seeded not equal to (colony area) * (number or corals)"
        @test all((seed_dist.TA .< 1) .&& (seed_dist.CA .< 1)) || "Some proportions of seeded corals greater than 1"
        @test all((seed_dist.TA .>= 0) .&& (seed_dist.CA .>= 0)) || "Some proportions of seeded corals less than zero"
        @test all((area_TA .<= selected_avail_space) .&& (area_CA .<= selected_avail_space)) || "Area seeded greater than available area"
        @test (max_ind_out == max_ind) || "Maximum distributed proportion of seeded coral not seeded in largest available area."
        @test (min_ind_out == min_ind) || "Minimum distributed proportion of seeded coral not seeded in smallest available area."
    end
end
