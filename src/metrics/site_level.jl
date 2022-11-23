using Statistics, OnlineStats


"""
    per_site(metric, data::NamedDimsArray)

Get metric results applied to the site-level at indicated time (or across timesteps).

# Arguments
- metric : Any function from the Statistics package to be applied to `data`
- data : Data set to apply metric to
- timestep : Target time step, or time frame

# Returns
Vector of N elements, where N is the number of sites.
"""
function per_site(metric, data::NamedDimsArray)
    num_sites = size(data, dim(data, :sites))
    met = zeros(num_sites)

    @inbounds for i in 1:num_sites
        met[i] = metric(data[sites=i])
    end

    return met
end