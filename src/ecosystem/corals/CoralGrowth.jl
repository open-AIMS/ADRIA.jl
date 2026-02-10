"""
    CoralGrowth(n_locs::Integer, n_groups::Integer, n_sizes::Integer, n_group_and_size::Integer)

Coral growth specification.
"""
struct CoralGrowth{A<:Integer}
    n_locs::A
    n_groups::A
    n_sizes::A
    n_group_and_size::A
end

"""
    CoralGrowth(n_locs)

Implements temporary hardcoded caches for a scenario with 35 'species' (split into 5 groups).
"""
function CoralGrowth(n_locs::Int64)
    n_groups::Int64, n_sizes::Int64 = size(linear_extensions())
    return CoralGrowth(n_locs, n_groups, n_sizes, n_groups * n_sizes)
end

#     @NamedTuple{
#         small::StaticArrays.SVector{5,Int64},  # indices for small size classes
#         mid::StaticArrays.SVector{25,Int64},   # indices for mid-size corals
#         large::StaticArrays.SVector{5,Int64}  # indices for large corals
#     }((
#     # Store specific indices for 5 functional groups and 7 size classes
#     @SVector [1, 8, 15, 22, 29],
#     SVector{25}(collect([2:6; 9:13; 16:20; 23:27; 30:34])),
#     @SVector [7, 14, 21, 28, 35]
# ))
