using Test
using ADRIA: distance_sorting
using LinearAlgebra
using Random
using Statistics

@testset "Distance Sorting" begin
    rng = MersenneTwister(1234)

    n_sites = 30
    top_n = 15
    dist_thresh = 0.1
    site_order = shuffle(rng, Vector(1:n_sites))
    prefsites = site_order[1:5]
    dists = rand(n_sites, n_sites)
    dists[diagind(dists)] .= NaN

    min_dist = median(dists[.!isnan.(dists)]) - dist_thresh * median(dists[.!isnan.(dists)])
    # only prefsites are < min_dist so only these should be replaced
    dists[dists.<min_dist] .= min_dist + 0.1
    dists[prefsites, prefsites] .= min_dist - 0.1
    pref_dists = findall(dists[prefsites, prefsites] .< min_dist)
    rep_sites = sort(unique(reinterpret(Int64, pref_dists)))

    new_prefsites = distance_sorting(prefsites, site_order, dists, dist_thresh, top_n)

    n_diff_sites = length(setdiff(prefsites, new_prefsites))
    @test all([in(new_prefsites[k], site_order[1:top_n]) for k = 1:5]) || "Not all sites are selected from the top_n."
    @test n_diff_sites == length(rep_sites) || "Some sites which should have been replaced have not been replaced."

    # set all sites above threshold and check none are replaced,
    dists[prefsites, prefsites] .= min_dist + 0.1
    pref_dists = findall(dists[prefsites, prefsites] .< min_dist)
    rep_sites = sort(unique(reinterpret(Int64, pref_dists)))

    new_prefsites = distance_sorting(prefsites, site_order, dists, dist_thresh, top_n)
    @test all(new_prefsites .== prefsites) || "All sites were sufficiently far apart but some were still replaced."

    # set all distances below threshold and check none are replaced.
    dists[dists.>min_dist] .= min_dist - 0.1
    pref_dists = findall(dists[prefsites, prefsites] .< min_dist)
    rep_sites = sort(unique(reinterpret(Int64, pref_dists)))
    new_prefsites = distance_sorting(prefsites, site_order, dists, dist_thresh, top_n)

    @test all(new_prefsites .== prefsites) || "All sites were too close according to the threshold but some were still replaced."
end