function scenario_clusters(clusters::BitVector)::Dict{Symbol,BitVector}
    return Dict(:target => clusters, :non_target => .!clusters)
end
function scenario_clusters(clusters::Vector{Int64})::Dict{Symbol,BitVector}
    return Dict(Symbol("Cluster_$(c)") => clusters .== c for c in unique(clusters))
end
