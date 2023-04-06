using Statistics, OnlineStats


"""
    per_location(metric, data::NamedDimsArray)

Get metric results applied to the location-level at indicated time (or across timesteps).

# Arguments
- metric : Any function from the Statistics package to be applied to `data`
- data : Data set to apply metric to
- timestep : Target time step, or time frame

# Returns
Vector of N elements, where N is the number of locations.
"""
function per_location(metric, data::NamedDimsArray)
    num_locations = size(data, dim(data, :locations))
    met = zeros(num_locations)

    @inbounds for i in 1:num_locations
        met[i] = metric(data[locations=i])
    end

    return met
end