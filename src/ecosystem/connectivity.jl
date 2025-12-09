"""
location_connectivity(file_loc::String, unique_loc_ids::Vector{String};
                  con_cutoff::Float64=1e-6, agg_func::Function=mean, swap::Bool=false)::NamedTuple

Create transitional probability matrix indicating connectivity between
sites, level of centrality, and the strongest predecessor for each site.

NOTE: Transposes transitional probability matrix if `swap == true`
      If multiple files are read in, this assumes all file rows/cols
      follow the same order as the first file read in.

# Examples
```julia
    location_connectivity("MooreTPmean.csv", site_order)
    location_connectivity("MooreTPmean.csv", site_order; con_cutoff=0.02, agg_func=mean, swap=true)
```

# Arguments
- `file_path` : Path to data file (or datasets) to load
               If a folder, searches subfolders as well
- `loc_ids` : Unique site ids in their expected order
- `conn_cutoff` : Percent thresholds of max for weak connections in
                 network (defined by user or defaults in `SimConstants`)
- `agg_func` : Summary statistic to take (defaults to `mean`)
- `swap` : Whether to transpose data (defaults to `false`)

# Returns
NamedTuple:
- `conn` : Matrix, containing the connectivity for all locations.
           Accounts for larvae which do not settle, so rows are not required to sum to 1
- `truncated` : ID of locations removed
- `loc_ids` : ID of locations kept
"""
function location_connectivity(
    file_path::String,
    loc_ids::Vector{String};
    conn_cutoff::Float64=1e-6,
    agg_func::Function=mean,
    swap::Bool=false
)::NamedTuple
    if !isdir(file_path) && !isfile(file_path)
        error("Could not find location: $(file_path)")
    end

    local extracted_TP::Matrix{Float64}
    if isfile(file_path)
        conn_files::Vector{String} = String[file_path]
        first_file = conn_files[1]
    elseif isdir(file_path)
        conn_fns = readdir(file_path)
        conn_fns = String[fn for fn in conn_fns if endswith(fn, ".csv")]

        first_file = joinpath.(file_path, conn_fns[1])

        # Assume years are always in the second position
        years::Vector{String} = unique(getindex.(split.(conn_fns, "_"), 2))

        # Organize files by their connectivity years
        year_conn_fns = NamedTuple{Tuple(Symbol.("year_" .* years))}(
        [filter(x -> occursin(yr, x), joinpath.(file_path, conn_fns)) for yr in years]
)

        # Create store for each year
        tmp_store::Vector{Matrix{Float64}} = Matrix{Float64}[]
        for yr in years
            assoc_files = getfield(year_conn_fns, Symbol.("year_" * yr))
            conn_data::Vector{Matrix{Float64}} = Matrix{Float64}[
                Matrix(
                    CSV.read(
                        fn,
                        DataFrame;
                        comment="#",
                        missingstring="NA",
                        transpose=swap,
                        types=Float64,
                        drop=[1]
                    )
                ) for fn in assoc_files
            ]

            push!(tmp_store, agg_func(conn_data))
        end

        # Mean across all years
        extracted_conn = agg_func(tmp_store)
    end

    # Get location ids from first file
    conn_file1 = CSV.read(
        first_file,
        DataFrame;
        comment="#",
        missingstring="NA",
        transpose=swap,
        types=Float64,
        drop=[1]
    )

    conn_loc_ids::Vector{String} = names(conn_file1)
    if isfile(file_path)
        extracted_conn = Matrix{Float64}(conn_file1)
    end

    # Identify locations that are not found in either conn_loc_ids or loc_ids
    # Get IDs that are in conn_loc_ids but not in loc_ids
    invalid_ids::Vector{String} = setdiff(conn_loc_ids, loc_ids)

    # Add locations that in loc_ids, but not in conn_loc_ids
    append!(invalid_ids, setdiff(loc_ids, conn_loc_ids))

    # Identify indices of IDs that do not appear in `invalid_ids`
    valid_ids::Vector{String} = setdiff(loc_ids, invalid_ids)
    valid_idx = findall(x -> x in loc_ids, valid_ids)

    # Align IDs
    loc_ids::Vector{String} = coalesce(loc_ids[valid_idx])
    loc_order = [findfirst(c_id .== conn_loc_ids) for c_id in loc_ids]

    if length(invalid_ids) > 0
        if length(invalid_ids) >= length(conn_loc_ids)
            error("All sites appear to be missing from data set. Aborting.")
        end

        @warn "The following sites (n=$(length(invalid_ids))) were not found in `loc_ids` and were removed:\n$(invalid_ids)"
    end

    # Reorder all data into expected form
    extracted_conn = extracted_conn[loc_order, loc_order]
    if conn_cutoff > 0.0
        extracted_conn[extracted_conn .< conn_cutoff] .= 0.0
    end

    conn = DataCube(extracted_conn; Source=loc_ids, Sink=loc_ids)

    @assert all(0.0 .<= conn .<= 1.0) "Connectivity data not scaled between 0 - 1"

    return (conn=conn, truncated=invalid_ids, loc_ids=loc_ids)
end
function location_connectivity(
    file_loc::String,
    unique_loc_ids::Vector{Union{Missing,String}};
    con_cutoff::Float64=1e-6,
    agg_func::Function=mean,
    swap::Bool=false
)::NamedTuple

    # Remove any row marked as missing
    if any(ismissing.(unique_loc_ids))
        @warn "Removing entries marked as `missing` from provided list of sites."
        unique_loc_ids::Vector{String} =
            String.(unique_loc_ids[.!ismissing.(unique_loc_ids)])
    else
        unique_loc_ids = String.(unique_loc_ids)
    end

    return location_connectivity(
        file_loc, unique_loc_ids; con_cutoff=con_cutoff, agg_func=agg_func, swap=swap
    )
end

function custom_indegree_centrality(conn_matrix::AbstractMatrix)
    # Count non-zero elements in each column (incoming connections)
    # This works for both dense and sparse matrices efficiently.
    indegrees = vec(sum(.!iszero.(conn_matrix); dims=1)) .|> Float64

    n = size(conn_matrix, 1)
    # Avoid division by zero if N=1
    if n > 1
        # Normalization factor for degree centrality is usually N-1
        # (the maximum number of possible connections from other nodes)
        indegrees ./= (n - 1)
    end

    return indegrees
end

# Overload for Graphs.SimpleWeightedDiGraph
function custom_indegree_centrality(g::SimpleWeightedDiGraph)
    # Extract the sparse adjacency matrix from the graph
    adj_matrix = Graphs.adjacency_matrix(g)
    # Call the matrix version of the function
    return custom_indegree_centrality(adj_matrix)
end

"""
    connectivity_strength(conn::AbstractArray)::NamedTuple

Determine in/out degree centralities for all nodes.

# Arguments
- `conn` : Base connectivity matrix to create Directed Graph from
- `in_method` : Centrality method used to determine incoming connectivity
- `out_method` : Centrality method used to determine outgoing connectivity

# Returns
NamedTuple:
- `in_conn` : Locations ranked by incoming connectivity
- `out_conn` : Locations ranked by outgoing connectivity
- `network` : Generated graph network
"""
function connectivity_strength(
    conn::AbstractMatrix{<:Union{Float32,Float64}};
    in_method=custom_indegree_centrality,
    out_method=eigenvector_centrality
)::NamedTuple
    g = SimpleWeightedDiGraph(conn)

    # Measure centrality based on number of incoming connections
    # Note: Graphs.indegree_centrality is unweighted.
    # If a weighted metric is desired, use `custom_indegree_centrality`.
    C1 = in_method(g)

    C2 = out_method(g)

    return (in_conn=C1, out_conn=C2, network=g)
end

"""
    connectivity_strength(area_weighted_conn::AbstractMatrix{Float64}, cover::Vector{Float64}, conn_cache::AbstractMatrix{Float64})::NamedTuple

Determine in/out degree centralities for all nodes.

# Arguments
- `area_weighted_conn` : Transfer probability matrix weighted by location `k` area
- `cover` : Total relative coral cover at location
- `conn_cache` : Cache matrix of same size as TP_base to hold intermediate values
- `in_method` : Centrality method used to determine incoming connectivity
- `out_method` : Centrality method used to determine outgoing connectivity

# Returns
NamedTuple:
- `in_conn` : Locations ranked by incoming connectivity
- `out_conn` : Locations ranked by outgoing connectivity
- `network` : Generated graph network
"""
function connectivity_strength(
    area_weighted_conn::AbstractMatrix{<:Union{Float32,Float64}},
    cover::Vector{<:Union{Float32,Float64}},
    conn_cache::AbstractMatrix{Float64};
    in_method=custom_indegree_centrality,
    out_method=eigenvector_centrality
)::NamedTuple
    # Accounts for cases where there is no coral cover
    conn_cache .= (area_weighted_conn .* cover)
    max_conn = maximum(conn_cache)
    if max_conn > 0.0
        conn_cache .= conn_cache ./ max_conn
    end

    return connectivity_strength(conn_cache; in_method, out_method)
end

"""
    strongest_source(g, indegree)::Vector{Int64}

For each node in the graph, identify the neighboring node with the strongest edge weight.
In other words, find the strongest source node.

# Arguments
- `g` : Graph network to assess
- `indegree` : vector of indegree connections for each node/vertex (taken from indegree_centrality)

# Returns
Vector of node ids indicating the strongest source for each node.
"""
function strongest_source(g, indegree)::Vector{Int64}
    # For each node, find strongly connected predecessor (by number of connections)
    strong_src = zeros(Int64, size(indegree)...)
    for v_id in vertices(g)
        incoming = inneighbors(g, v_id)

        if length(incoming) > 0
            # For each incoming connection, find the one with most "in"
            # connections themselves
            in_conns = Int64[length(inneighbors(g, in_id)) for in_id in incoming]

            # Find index of predecessor with most connections
            # (use `first` to get the first match in case of a tie)
            most_conns = maximum(in_conns)
            idx = first(findall(in_conns .== most_conns))
            strong_src[v_id] = incoming[idx]
        else
            strong_src[v_id] = 0
        end
    end

    return strong_src
end
