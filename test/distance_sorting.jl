using Test
using ADRIA: distance_sorting
using LinearAlgebra
using Random
using Statistics

@testset "Distance Sorting" begin

    n_locations = 30
    dist_thresh = 0.1
    location_order = Union{Float64,Int64}[shuffle(collect(1:n_locations)) collect(1.0:-1/30:0.0)[1:n_locations] 1:n_locations]

    preflocations = Int.(location_order[1:5, 1])
    new_locations = Int.(location_order[6:10, 1])
    dists = rand(n_locations, n_locations)
    dists[diagind(dists)] .= NaN

    min_dist = dist_thresh * median(dists[.!isnan.(dists)])

    # only preflocations are < min_dist so only these should be replaced
    dists .= min_dist + 0.1
    dists[preflocations, preflocations] .= min_dist - 0.1
    pref_dists = findall(dists[preflocations, preflocations] .< min_dist)
    rep_locations = sort(unique(reinterpret(Int64, pref_dists)))

    rankings = zeros(Int64, length(location_order[:, 3]), 2)
    rankings .= location_order[:, 3]
    new_preflocations, new_rankings = distance_sorting(preflocations, location_order, dists, dist_thresh, rankings, 2)

    n_diff_locations = length(setdiff(preflocations, new_preflocations))
    @test n_diff_locations == length(rep_locations) || "Some locations which should have been replaced have not been replaced."
    @test all(new_preflocations .== new_locations) || "Higher ranked locations satisfying the distance threshold were available but not used to replace unsuitable sites."

    # set all locations above threshold and check none are replaced,
    dists .= min_dist - 0.1
    dists[preflocations, preflocations] .= min_dist + 0.1
    pref_dists = findall(dists[preflocations, preflocations] .< min_dist)
    rep_locations = sort(unique(reinterpret(Int64, pref_dists)))

    new_preflocations, new_rankings = distance_sorting(preflocations, location_order, dists, dist_thresh, rankings, 2)
    @test all(new_preflocations .== preflocations) || "All locations were sufficiently far apart but some were still replaced."

    # set all distances below threshold and check none are replaced.
    dists .= min_dist - 0.1
    pref_dists = findall(dists[preflocations, preflocations] .< min_dist)
    rep_locations = sort(unique(reinterpret(Int64, pref_dists)))
    new_preflocations, new_rankings = distance_sorting(preflocations, location_order, dists, dist_thresh, rankings, 2)

    @test all(new_preflocations .== preflocations) || "All locations were too close according to the threshold but some were still replaced."
end