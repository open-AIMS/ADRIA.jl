using NamedArrays
using ADRIA: ResultSet


"""
    seeded_sites_frequency(rs::ResultSet,scens::NamedTuple)::NamedTuple

Count frequency of seeded sites for scenarios satisfying a condition.

# Arguments
- 'rs' : ResultSet
- `scens` : contains scenario ids for scenarios satisfying the condition of interest.

# Returns
NamedArray, where each entry relates to an RCP of interest, e.g., `[RCP45=[... frequency of selection for each site ...]; 
RCP60=[ ... frequency of selection for each site ...]]`

"""
function seeded_sites_frequency(rs::ResultSet, scens::NamedTuple)::NamedArray

    # retrieve RCPs
    rcps = keys(scens)
    # create frequencies storage container
    seeded_sites_store = NamedArray(zeros(length(rcps), size(rs.site_data, 1)))
    idx_rows = [String(i) for i = rcps]
    setnames!(seeded_sites_store, idx_rows, 1)
    for rcp in rcps
        ind_cond_temp = scens[rcp]

        # select scenarios satisfying condition and sum up selection tally for each site
        seed_log = dropdims(sum(rs.seed_log[:, :, :, ind_cond_temp], dims=2), dims=2)
        seeded_sites_store[String(rcp), :] .= vec(dropdims(sum(seed_log .> 0, dims=[1, 3]), dims=3))
    end

    return seeded_sites_store
end
