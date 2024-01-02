using Test
using ADRIA: decision.constrain_reef_cluster
using LinearAlgebra
using Random
using Statistics

@testset "Constrain spatial groups" begin
	n_sites = 30
	site_ids = collect(1:n_sites)
	area_to_seed = 695.11
	orig_site_order = shuffle(Vector(1:n_sites))
	site_order = copy(orig_site_order)
	available_space = rand(Uniform(area_to_seed + 100.0, area_to_seed + 1000.0), 30)
	n_iv_locs = 5

	prefsites = site_order[1:5]
	reef_locs = [fill("1", 10)..., fill("2", 10)..., fill("3", 10)...]

	s_order = Union{Float64, Int64}[Int64.(site_order) rand(n_sites)]
	# All selected sites are in the same reef, so 2 should be replaced
	reef_locs[prefsites] .= "2"
	# Empty ranking just for testing
	rankings = Int64[site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

	new_prefsites, rankings = constrain_reef_cluster(
		reef_locs, 
		s_order, 
		rankings,
		area_to_seed,
		available_space,
		n_iv_locs,
		3,)

	num_reefs = [sum(reef_locs[new_prefsites] .== rr) for rr in unique(reef_locs)]
	l_diff_sites = length(setdiff(prefsites, new_prefsites))

	@test all(l_diff_sites == 2) ||
		  "Too few or too many sites have been removed when constraining spatial groups."
	@test all(prefsites[1:3] .== new_prefsites[1:3]) ||
		  "Some sites which should not have been replaced have been replaced when constraining spatial groups."
	@test all(num_reefs .<= 3) ||
		  "More sites than allowed by the spatial group constraint are in the same reef."

	# Set no 3 sites to be in the same reef and check none are replaced
	reef_locs[prefsites] .= ["1", "2", "3", "4", "1"]

	s_order = Union{Float64, Int64}[Int64.(orig_site_order) rand(n_sites)]
	rankings = Int64[site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

	new_prefsites, rankings = constrain_reef_cluster(
		reef_locs,
		s_order, 
		rankings, 
		area_to_seed, 
		available_space, 
		n_iv_locs, 
		3,
	)

	@test all(new_prefsites .== prefsites) ||
		  "All sites in different reefs but some were still replaced."

	# Make slected sites not have enough space to seed corals
	available_space[prefsites] .= (area_to_seed - 100.0) / n_iv_locs
	available_space[s_order[n_iv_locs+1, 1]] = area_to_seed

	s_order = Union{Float64, Int64}[Int64.(orig_site_order) rand(n_sites)]
	rankings = Int64[site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

	new_prefsites, rankings = constrain_reef_cluster(
		reef_locs, s_order, rankings, area_to_seed, available_space, n_iv_locs, 3,
	)
	@test length(new_prefsites) == (n_iv_locs + 1) ||
		  "Not enough sites were selected to fit the corals to be seeded."

end
