
function seed_ranks(rs::ResultSet; kwargs...)
    selected = slice_results(rs.ranks[intervention=1]; kwargs...)
    nsteps, nsites = size(selected)

    @assert length(timesteps) == nsteps

    r_ids = rs.site_data.reef_siteid
    min_rank = length(r_ids) + 1
    if haskey(kwargs, :sites)
        r_ids = r_ids[kwargs[:sites]]
    end

    @assert length(r_ids) == nsites

    return min_rank .- NamedArray(selected, (timesteps(rs), r_ids))

    # n = NamedArray([1 3; 2 4], ( OrderedDict("A"=>1, "B"=>2), OrderedDict("C"=>1, "D"=>2) ),
    #            ("Rows", "Cols"))
    # @show n;

    # n = 2×2 Named Array{Int64,2}
    # Rows ╲ Cols │ C  D
    # ────────────┼─────
    # A           │ 1  3
    # B           │ 2  4

    # n = NamedArray([1 2 3; 4 5 6], (["one", "two"], [:a, :b, :c]))

    # return partialsortperm(vec(mean(rs.ranks[intervention=1, scenarios=5], dims=1)), 1:10)
end

function top_n_seeded_sites(rs::ResultSet; kwargs...)
    ranked_sites = seed_ranks(rs; kwargs...)
    
end