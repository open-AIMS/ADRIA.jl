"""
site_connectivity(file_loc::String, unique_site_ids::Vector{String};
                  con_cutoff::Float64=1e-6, agg_func::Function=mean, swap::Bool=false)::NamedTuple

Create transitional probability matrix indicating connectivity between
sites, level of centrality, and the strongest predecessor for each site.

NOTE: Transposes transitional probability matrix if `swap == true`
      If multiple files are read in, this assumes all file rows/cols
      follow the same order as the first file read in.

# Examples
```julia
    site_connectivity("MooreTPmean.csv", site_order)
    site_connectivity("MooreTPmean.csv", site_order; con_cutoff=0.02, agg_func=mean, swap=true)
```

# Arguments
- `file_loc` : Path to data file (or datasets) to load
               If a folder, searches subfolders as well
- `unique_site_ids` : Unique site ids in their expected order
- `con_cutoff` : Percent thresholds of max for weak connections in
                 network (defined by user or defaults in `SimConstants`)
- `agg_func` : Summary statistic to take (defaults to `mean`)
- `swap` : Whether to transpose data (defaults to `false`)

# Returns
NamedTuple:
- `TP_data` : Matrix, containing the transition probability for all sites
- `truncated` : ID of sites removed
- `site_ids` : ID of sites kept
"""
function site_connectivity(file_loc::String, unique_site_ids::Vector{String};
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

    # Get site ids from first file
    con_file1::DataFrame = CSV.read(con_files[1], DataFrame, comment="#", missingstring="NA", transpose=swap, types=Float64, drop=[1])
    con_site_ids::Vector{String} = String[x[1] for x in split.(names(con_file1), "_v"; limit=2)]

    if isfile(file_loc)
        extracted_TP = Matrix{Float64}(con_file1)
    end

    # Get IDs missing in con_site_ids
    invalid_ids::Vector{String} = setdiff(con_site_ids, unique_site_ids)

    # Get IDs missing in site_order
    append!(invalid_ids, setdiff(unique_site_ids, con_site_ids))

    # Identify IDs that do not appear in `invalid_ids`
    valid_ids::Vector{String} = [x ∉ invalid_ids ? x : missing for x in unique_site_ids]
    valid_idx = .!ismissing.(valid_ids)

    # Align IDs
    unique_site_ids::Vector{String} = coalesce(unique_site_ids[valid_idx])
    site_order = [findfirst(c_id .== con_site_ids) for c_id in unique_site_ids]

    if length(invalid_ids) > 0
        if length(invalid_ids) >= length(con_site_ids)
            error("All sites appear to be missing from data set. Aborting.")
        end

        @warn "The following sites (n=$(length(invalid_ids))) were not found in site_ids and were removed:\n$(invalid_ids)"
    end

    # Reorder all data into expected form
    extracted_TP = extracted_TP[site_order, site_order]

    if con_cutoff > 0.0
        extracted_TP[extracted_TP.<con_cutoff] .= 0.0
    end

    TP_base = NamedDimsArray(sparse(extracted_TP), Source=unique_site_ids, Receiving=unique_site_ids)
    @assert all(0.0 .<= TP_base .<= 1.0) "Connectivity data not scaled between 0 - 1"

    return (TP_base=TP_base, truncated=invalid_ids, site_ids=unique_site_ids)
end
function site_connectivity(file_loc::String, unique_site_ids::Vector{Union{Missing,String}};
    con_cutoff::Float64=1e-6, agg_func::Function=mean, swap::Bool=false)::NamedTuple

    # Remove any row marked as missing
    if any(ismissing.(unique_site_ids))
        @warn "Removing entries marked as `missing` from provided list of sites."
        unique_site_ids::Vector{String} = String.(unique_site_ids[.!ismissing.(unique_site_ids)])
    else
        unique_site_ids = String.(unique_site_ids)
    end

    return site_connectivity(file_loc, unique_site_ids;
        con_cutoff=con_cutoff, agg_func=agg_func, swap=swap)
end


"""
    connectivity_strength(TP_base::AbstractArray)::NamedTuple

Generate array of outdegree connectivity strength for each node and its
strongest predecessor.

# Returns
NamedTuple:
- `in_conn` : sites ranked by incoming connectivity
- `out_conn` : sites ranked by outgoing connectivity
- `strongest_predecessor` : strongest predecessor for each site
"""
function connectivity_strength(TP_base::AbstractArray)::NamedTuple

    g = SimpleDiGraph(TP_base)

    # ew_base = weights(g)  # commented out ew_base are all equally weighted anyway...

    # Measure centrality based on number of incoming connections
    # C1 = indegree_centrality(g)
    # C2 = outdegree_centrality(g)
    C1 = betweenness_centrality(g)
    C2 = stress_centrality(g)

    # strong_pred = closeness_centrality(g)

    # For each edge, find strongly connected predecessor (by number of connections)
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
