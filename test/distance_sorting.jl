using Test
using ADRIA: distance_sorting
using LinearAlgebra
using Random
using Statistics

@testset "Distance Sorting" begin

    n_locations = 30
    top_n = 15
    dist_thresh = 0.1
    location_order = shuffle(Vector(1:n_locations))
    preflocations = location_order[1:5]
    dists = rand(n_locations, n_locations)
    dists[diagind(dists)] .= NaN

    min_dist = median(dists[.!isnan.(dists)]) - dist_thresh * median(dists[.!isnan.(dists)])
    # only preflocations are < min_dist so only these should be replaced
    dists[dists.<min_dist] .= min_dist + 0.1
    dists[preflocations, preflocations] .= min_dist - 0.1
    pref_dists = findall(dists[preflocations, preflocations] .< min_dist)
    rep_locations = sort(unique(reinterpret(Int64, pref_dists)))

    new_preflocations = distance_sorting(preflocations, location_order, dists, dist_thresh, top_n)

    n_diff_locations = length(setdiff(preflocations, new_preflocations))
    @test all([in(new_preflocations[k], location_order[1:top_n]) for k = 1:5]) || "Not all locations are selected from the top_n."
    @test n_diff_locations == length(rep_locations) || "Some locations which should have been replaced have not been replaced."

    # set all locations above threshold and check none are replaced,
    dists[preflocations, preflocations] .= min_dist + 0.1
    pref_dists = findall(dists[preflocations, preflocations] .< min_dist)
    rep_locations = sort(unique(reinterpret(Int64, pref_dists)))

    new_preflocations = distance_sorting(preflocations, location_order, dists, dist_thresh, top_n)
    @test all(new_preflocations .== preflocations) || "All locations were sufficiently far apart but some were still replaced."

    # set all distances below threshold and check none are replaced.
    dists[dists.>min_dist] .= min_dist - 0.1
    pref_dists = findall(dists[preflocations, preflocations] .< min_dist)
    rep_locations = sort(unique(reinterpret(Int64, pref_dists)))
    new_preflocations = distance_sorting(preflocations, location_order, dists, dist_thresh, top_n)

    @test all(new_preflocations .== preflocations) || "All locations were too close according to the threshold but some were still replaced."
end