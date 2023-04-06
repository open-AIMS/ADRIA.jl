"""
location_connectivity(file_loc::String, unique_location_ids::Vector{String};
                  con_cutoff::Float64=1e-6, agg_func::Function=mean, swap::Bool=false)::NamedTuple

Create transitional probability matrix indicating connectivity between
locations, level of centrality, and the strongest predecessor for each location.

NOTE: Transposes transitional probability matrix if `swap == true`
      If multiple files are read in, this assumes all file rows/cols
      follow the same order as the first file read in.

# Examples
```julia
    location_connectivity("MooreTPmean.csv", location_order)
    location_connectivity("MooreTPmean.csv", location_order; con_cutoff=0.02, agg_func=mean, swap=true)
```

# Arguments
- `file_loc` : Path to data file (or datasets) to load
               If a folder, searches subfolders as well
- `unique_location_ids` : Unique location ids in their expected order
- `con_cutoff` : Percent thresholds of max for weak connections in
                 network (defined by user or defaults in `SimConstants`)
- `agg_func` : Summary statistic to take (defaults to `mean`)
- `swap` : Whether to transpose data (defaults to `false`)

# Returns
NamedTuple:
- `TP_data` : Matrix, containing the transition probability for all locations
- `truncated` : ID of locations removed
- `location_ids` : ID of locations kept
"""
function location_connectivity(file_loc::String, unique_location_ids::Vector{String};
    con_cutoff::Float64=1e-6, agg_func::Function=mean, swap::Bool=false)::NamedTuple

    if !isdir(file_loc) && !isfile(file_loc)
        error("Could not find location: $(file_loc)")
    end

    local extracted_TP::Matrix{Float64}
    if isfile(file_loc)
        con_files::Vector{String} = String[file_loc]
    elseif isdir(file_loc)
        # Get connectivity years available in data store
        years::Vector{String} = getindex(first(walkdir(file_loc)), 2)
        year_conn_fns = NamedTuple{Tuple(Symbol.(years))}(
            [[joinpath.(first(fl), last(fl))
              for fl in walkdir(joinpath(file_loc, yr))][1]
             for yr in years]
        )

        con_files = vcat([x for x in values(year_conn_fns)]...)

        # Pre-allocate store
        tmp_store::Vector{Matrix{Float64}} = Vector{Matrix{Float64}}(undef, length(years))

        # Get average connectivity for each represented year
        for (i, yr) in enumerate(Symbol.(years))
            conn_data::Vector{Matrix{Float64}} = Matrix{Float64}[
                Matrix(CSV.read(fn, DataFrame, comment="#", missingstring="NA", transpose=swap, types=Float64, drop=[1]))
                for fn in year_conn_fns[yr]
            ]

            tmp_store[i] = agg_func(conn_data)
        end

        # Mean across all years
        extracted_TP = agg_func(tmp_store)
    end

    # Get location ids from first file
    con_file1::DataFrame = CSV.read(con_files[1], DataFrame, comment="#", missingstring="NA", transpose=swap, types=Float64, drop=[1])
    con_location_ids::Vector{String} = String[x[1] for x in split.(names(con_file1), "_v"; limit=2)]

    if isfile(file_loc)
        extracted_TP = Matrix{Float64}(con_file1)
    end

    # Get IDs missing in con_location_ids
    invalid_ids::Vector{String} = setdiff(con_location_ids, unique_location_ids)

    # Get IDs missing in location_order
    append!(invalid_ids, setdiff(unique_location_ids, con_location_ids))

    # Identify IDs that do not appear in `invalid_ids`
    valid_ids::Vector{String} = [x ∉ invalid_ids ? x : missing for x in unique_location_ids]
    valid_idx = .!ismissing.(valid_ids)

    # Align IDs
    unique_location_ids::Vector{String} = coalesce(unique_location_ids[valid_idx])
    location_order = [findfirst(c_id .== con_location_ids) for c_id in unique_location_ids]

    if length(invalid_ids) > 0
        if length(invalid_ids) >= length(con_location_ids)
            error("All locations appear to be missing from data set. Aborting.")
        end

        @warn "The following locations (n=$(length(invalid_ids))) were not found in location_ids and were removed:\n$(invalid_ids)"
    end

    # Reorder all data into expected form
    extracted_TP = extracted_TP[location_order, location_order]

    if con_cutoff > 0.0
        extracted_TP[extracted_TP.<con_cutoff] .= 0.0
    end

    TP_base = NamedDimsArray(sparse(extracted_TP), Source=unique_location_ids, Receiving=unique_location_ids)
    @assert all(0.0 .<= TP_base .<= 1.0) "Connectivity data not scaled between 0 - 1"

    return (TP_base=TP_base, truncated=invalid_ids, location_ids=unique_location_ids)
end
function location_connectivity(file_loc::String, unique_location_ids::Vector{Union{Missing,String}};
    con_cutoff::Float64=1e-6, agg_func::Function=mean, swap::Bool=false)::NamedTuple

    # Remove any row marked as missing
    if any(ismissing.(unique_location_ids))
        @warn "Removing entries marked as `missing` from provided list of locations."
        unique_location_ids::Vector{String} = String.(unique_location_ids[.!ismissing.(unique_location_ids)])
    else
        unique_location_ids = String.(unique_location_ids)
    end

    return location_connectivity(file_loc, unique_location_ids;
        con_cutoff=con_cutoff, agg_func=agg_func, swap=swap)
end


"""
    connectivity_strength(TP_base::AbstractArray)::NamedTuple
    connectivity_strength(area_weighted_TP::AbstractMatrix{Float64}, cover::Vector{Float64})::NamedTuple

Create in/out degree centralities for all nodes, and vector of their strongest predecessors.

# Arguments
- `TP_base` : Base transfer probability matrix to create Directed Graph from.
- `area_weighted_TP` : Transfer probability matrix weighted by location `k` area
- `cover` : total coral cover at location

# Returns
NamedTuple:
- `in_conn` : locations ranked by incoming connectivity
- `out_conn` : locations ranked by outgoing connectivity
- `strongest_predecessor` : strongest predecessor for each location
"""
function connectivity_strength(TP_base::AbstractMatrix{Float64})::NamedTuple

    g = SimpleDiGraph(TP_base)

    # Measure centrality based on number of incoming connections
    C1 = indegree_centrality(g)
    C2 = outdegree_centrality(g)

    # For each node, find strongly connected predecessor (by number of connections)
    strong_pred = zeros(Int64, size(C1)...)
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
            strong_pred[v_id] = incoming[idx]
        else
            strong_pred[v_id] = 0
        end
    end

    return (in_conn=C1, out_conn=C2, strongest_predecessor=strong_pred)
end
function connectivity_strength(area_weighted_TP::AbstractMatrix{Float64}, cover::Vector{Float64})::NamedTuple

    # Accounts for cases where there is no coral cover
    tp = (area_weighted_TP .* cover)
    tp .= tp ./ maximum(tp)

    return connectivity_strength(tp)
end
