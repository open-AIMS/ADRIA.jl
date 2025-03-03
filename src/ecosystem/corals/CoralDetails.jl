"""
    CoralDetails{A<:Int64}

Convenience store of coral groups and size classes.
"""
struct CoralDetails{A<:Int64}
    n_groups::A
    n_sizes::A
    n_group_and_size::A
end

"""
    CoralDetails()::CoralDetails
    CoralDetails(n_groups::Int64, n_sizes::Int64)::CoralDetails

Convenience store of coral groups and size classes.
"""
function CoralDetails()
    n_groups::Int64, n_sizes::Int64 = size(linear_extensions())
    return CoralDetails(n_groups, n_sizes, n_groups * n_sizes)
end
function CoralDetails(n_groups::Int64, n_sizes::Int64)::CoralDetails
    n_group_and_size::Int64 = n_groups * n_sizes

    return CoralDetails(n_groups, n_sizes, n_group_and_size)
end
