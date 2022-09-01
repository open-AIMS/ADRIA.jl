using ADRIA
import ADRIA: distribute_seeded_corals
using Test
using Distributions
using Infiltrator

@testset "seeding" begin
    # first test function on example domain
    dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)

    @infiltrate
    # extract inputs for function
    total_site_area = site_area(dom)
    k = dom.site_data.k / 100
    current_cover = zeros(size(total_site_area))
    # calculate available space
    available_space = vec((k .- current_cover))
    prefseedsites = vec(rand(1:length(total_site_area),1,5))
    col_area_seed = vec(rand(Uniform(0.0,0.005),1,2))
    n_to_seed = vec(rand(20000:100000,1,2))

    # evaluate seeding distributions
    seed_dist = distribute_seeded_corals(total_site_area,
        prefseedsites, available_space,
        n_to_seed, col_area_seed)

    # proportions of coral 
    total_area_coral_TA = sum((n_to_seed[1] .* col_area_seed[1]))
    total_area_coral_CA = sum((n_to_seed[2] .* col_area_seed[2]))

    # Area to be seeded for each site
    area_TA = seed_dist[1] .* total_site_area[prefseedsites]
    area_CA = seed_dist[2] .* total_site_area[prefseedsites]

    # total area of seeded corals
    total_area_coral_TA_out = sum(area_TA)
    total_area_coral_CA_out = sum(area_CA)

    # absolute available area to seed
    available_area = available_space[prefseedsites] .* total_site_area[prefseedsites]

    # index of max proportion for available space
    max_ind_out = findfirst(item -> item == maximum(seed_dist[1]), seed_dist[1])
    max_ind = findfirst(item -> item == maximum(available_space[prefseedsites]), available_space[prefseedsites])

    # index of min proportion for available space
    min_ind_out = findfirst(item -> item == minimum(seed_dist[1]), seed_dist[1])
    min_ind = findfirst(item -> item == minimum(available_space[prefseedsites]), available_space[prefseedsites])

    # run tests
    @test ((total_area_coral_TA - total_area_coral_TA_out)<10^-5) && ((total_area_coral_CA - total_area_coral_CA_out)<10^-5) || "Area of corals seeded not equal to (colony area) * (number or corals)"
    @test all((seed_dist[1] .< 1) .&& (seed_dist[2] .< 1)) || "Some proportions of seeded corals greater than 1"
    @test all((seed_dist[1] .>= 0) .&& (seed_dist[2] .>= 0)) || "Some proportions of seeded corals less than zero"
    @test all((area_TA .<= available_area) .&& (area_CA .<= available_area)) || "Area seeded greater than available area"
    @test (max_ind_out == max_ind) || "Maximum distributed proportion of seeded coral not seeded in largest available area."
    @test (min_ind_out == min_ind) || "Minimum distributed proportion of seeded coral not seeded smallest available area."

end