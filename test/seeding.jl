using ADRIA
import ADRIA: distribute_seeded_corals
using Test
using Distributions

@testset "Seeding distribution function" begin
    # first test function on example domain
    dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)

    # extract inputs for function
    total_site_area = site_area(dom)
    k = dom.site_data.k / 100
    current_cover = zeros(size(total_site_area))
    # calculate available space
    available_space = vec(((total_site_area .* k) .- current_cover))

    prefseedsites = rand(1:length(total_site_area), 5)

    # Randomly generate seeded area
    tmp = rand(Uniform(0.0, 500.0), 2)
    seeded_area = (TA=tmp[1], CA=tmp[2])

    # evaluate seeding distributions
    seed_dist = distribute_seeded_corals(total_site_area, prefseedsites, available_space, seeded_area)

    # proportions of coral 
    total_area_coral_TA = seeded_area[1]
    total_area_coral_CA = seeded_area[2]

    # Area to be seeded for each site
    area_TA = seed_dist.TA .* total_site_area[prefseedsites]
    area_CA = seed_dist.CA .* total_site_area[prefseedsites]

    # total area of seeded corals
    total_area_coral_TA_out = sum(area_TA)
    total_area_coral_CA_out = sum(area_CA)

    # absolute available area to seed for selected sites
    selected_avail_space = available_space[prefseedsites]

    # index of max proportion for available space
    max_ind_out = findfirst(seed_dist.TA .== maximum(seed_dist.TA))
    max_ind = findfirst(selected_avail_space .== maximum(selected_avail_space))

    max_ind_out2 = area_TA[max_ind_out]

    # index of min proportion for available space
    min_ind_out = findfirst(seed_dist.TA .== minimum(seed_dist.TA))
    min_ind = findfirst(selected_avail_space .== minimum(selected_avail_space))

    # run tests
    @test ((total_area_coral_TA - total_area_coral_TA_out) < 10^-5) && ((total_area_coral_CA - total_area_coral_CA_out) < 10^-5) || "Area of corals seeded not equal to (colony area) * (number or corals)"
    @test all((seed_dist.TA .< 1) .&& (seed_dist.CA .< 1)) || "Some proportions of seeded corals greater than 1"
    @test all((seed_dist.TA .>= 0) .&& (seed_dist.CA .>= 0)) || "Some proportions of seeded corals less than zero"
    @test all((area_TA .<= selected_avail_space) .&& (area_CA .<= selected_avail_space)) || "Area seeded greater than available area"

    @test max_ind_out == max_ind || "Maximum distributed proportion of seeded coral not seeded in largest available area."
    @test min_ind_out == min_ind || "Minimum distributed proportion of seeded coral not seeded in smallest available area."
end