using Combinatorics


"""
    seeded_sites_robust_scens(rs::ResultSet, y::AbstractArray, rcps::Vector{Int} offset::Int=0)    

Count frequency of seeded sites for robust and non-robust scenarios where robustness is defined by a 
    level of output metric(s) achieved.

# Arguments
- `rs` : ResultSet
- `y` : scenario outcomes
- `rcps` : RCP ids as integers (e.g., 45, 60, 85)
- `offset` : include scenarios that are some distance from pareto front, where 0 refers to
             the pareto front itself.

# Returns
NamedTuple, where each entry relates to an RCP of interest, e.g., `(RCP45=[... scenario ids ...], RCP60=[ ...scenario_ids ...])`

"""
function seeded_sites_robust_scens(scens::NamedTuple, nscens::Int64, rcps::Vector{Int})::NamedTuple

    for rcp in rcps
        ind_metrics_robust = scens["RCP$rcp"]
        ind_metrics_not_robust = setdiff(collect(Int64, 1:nscens), ind_metrics_robust)

        seed_log = dropdims(sum(rs.seed_log, dims=2), dims=2)
        seeded_sites = zeros(Int64, (size(seed_log)[1], 5, size(seed_log)[3]))

        for k in 1:size(seed_log)[1]
            for j in 1:size(seed_log)[3]
                sites = findall(seed_log[k, :, j] .> 0)
                if !isempty(sites)
                    seeded_sites[k, 1:length(sites), j] .= sites
                end
            end
        end
    end

end
