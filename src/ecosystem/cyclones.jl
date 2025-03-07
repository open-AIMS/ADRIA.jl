"""
    cyclone_mortality!(coral_cover, coral_params, cyclone_mortality)::Nothing
    cyclone_mortality!(coral_cover::AbstractArray{Float64, 3}, cyclone_mortality::AbstractMatrix{Float64})::Nothing

Apply cyclone mortalities.

# Arguments
- `coral_cover` : Coral cover for current time step
- `coral_params` : Coral parameters indicating indices of small/mid/large size classes
- `cyclone_mortality` : Mortalities for each functional group and size class
"""
function cyclone_mortality!(coral_cover, coral_params, cyclone_mortality)::Nothing
    # TODO: Move to own file.

    # Small class coral mortality
    coral_deaths_small = coral_cover[:, 1, :] .* cyclone_mortality
    coral_cover[coral_cover, 1, :] -= coral_deaths_small

    # Mid class coral mortality
    coral_mid = hcat(
        collect(Iterators.partition(coral_params.mid, length(coral_params.small)))...
    )
    for i in axes(coral_mid, 1)
        coral_deaths_mid = coral_cover[coral_mid[i, :], :] .* cyclone_mortality
        coral_cover[coral_mid[i, :], :] -= coral_deaths_mid
    end

    # Large class coral mortality
    coral_deaths_large = coral_cover[coral_params.large, :] .* cyclone_mortality
    coral_cover[coral_params.large, :] -= coral_deaths_large

    # Ensure no negative values
    clamp!(coral_cover, 0.0, 1.0)

    return nothing
end
function cyclone_mortality!(
    coral_cover::AbstractArray{Float64,3}, cyclone_mortality::AbstractMatrix{Float64}
)::Nothing
    n_groups, n_locs = size(cyclone_mortality)
    coral_cover .= coral_cover .* (1 .- reshape(cyclone_mortality, (n_groups, 1, n_locs)))
    clamp!(coral_cover, 0.0, 1.0)
    return nothing
end