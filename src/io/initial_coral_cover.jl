using Distributions,
    NetCDF,
    YAXArrays

"""
    load_initial_cover(data_fn::String)::YAXArray
    load_initial_cover(n_species::Int64, n_locs::Int64)::YAXArray

Load initial coral cover data from netCDF.
"""
function load_initial_cover(data_fn::String)::YAXArray
    _dim_names_replace = [:covers => :species, :reef_siteid => :sites]
    return _split_cover(
        load_nc_data(data_fn, "covers"; dim_names_replace=_dim_names_replace)
    )
end
function load_initial_cover(n_group_and_size::Int64, n_locs::Int64)::YAXArray
    @warn "Using random initial coral cover"
    random_cover_data = rand(Float32, n_group_and_size, n_locs)
    return DataCube(random_cover_data; species=1:n_group_and_size, sites=1:n_locs)
end

"""
    function split_cover(cover::YAXArray)::YAXArray

Split initial cover for each functional group into size classes assuming diameter follows
a LogNormal distribution.
"""
function _split_cover(cover::YAXArray)::YAXArray
    # Initial cover fraction for each functional group and size class
    init_cover_weights::Matrix{Float64} = _init_cover_weights()

    # Select only functional groups that ADRIA is currently working with
    selected_groups::BitVector = _select_functional_groups(cover)

    # Build initial cover final data
    n_groups::Int64, n_sizes::Int64 = size(init_cover_weights)
    n_groups_sizes::Int64 = n_groups * n_sizes
    cover_data::Matrix{Float64} =
        reshape(init_cover_weights', n_groups_sizes) .*
        repeat(cover[selected_groups, :]; inner=(n_sizes, 1))

    return DataCube(cover_data; _cover_labels(cover, n_sizes, selected_groups)...)
end

function _cover_labels(
    cover::YAXArray, n_sizes::Int64, selected_groups::BitVector
)::@NamedTuple{species::Vector{String}, locations::Vector{String}}
    taxa_names::Vector{String} = cover.species.val.data[selected_groups]
    tn::Vector{String} = repeat(taxa_names; inner=n_sizes)
    taxa_id::Vector{Int64} = repeat(1:n_sizes; inner=n_sizes)
    size_id::Vector{Int64} = repeat(1:n_sizes, n_sizes)
    species_labels::Vector{String} = String[join(x, "_") for x in zip(tn, taxa_id, size_id)]

    return (species=species_labels, locations=cover.locations.val.data)
end

function _reef_mod_area_dist(bin_edges::Matrix{Float64})::Vector{LogNormal{Float64}}
    # Proportionally adjust dist mean and std for each functional group based on old values
    base_mean::Float64 = 700.0
    base_std::Float64 = 4.0
    base_final_bin_edge::Float64 = 80.0
    scale_factor = bin_edges[:, end] ./ base_final_bin_edge

    # The mean scale factor is squared because area is proportional to the diameter squared
    dist_mean::Vector{Float64} = base_mean .* (scale_factor .^ 2)
    dist_std::Vector{Float64} = base_std .* scale_factor

    return LogNormal.(log.(dist_mean), log.(dist_std))
end

"""
    _init_cover_weights()::Matrix{Float64}

Initial cover fraction for each functional group and size class
"""
function _init_cover_weights()::Matrix{Float64}
    _bin_edges::Matrix{Float64} = bin_edges(; unit=:cm)

    # Build distinct LogNormal dists for each functional group
    reef_mod_area_dist::Vector{LogNormal{Float64}} = _reef_mod_area_dist(_bin_edges)

    bin_edges_area::Matrix{Float64} = colony_mean_area(_bin_edges)
    cdf_integral::Matrix{Float64} = cdf.(reef_mod_area_dist, bin_edges_area)
    init_cover_fracs::Matrix{Float64} = (
        cdf_integral[:, 2:end] .- cdf_integral[:, 1:(end - 1)]
    )
    init_cover_fracs = init_cover_fracs ./ sum(init_cover_fracs; dims=2)
    return replace!(init_cover_fracs, NaN => 0.0)
end

"""
    _select_functional_groups(cover::YAXArray)::BitVector

Mask to select only functional groups that ADRIA is currently working with
"""
function _select_functional_groups(cover::YAXArray)::BitVector
    groups::Vector{String} = replace.(lowercase.(cover.species.val.data), (' ', '-') => '_')
    considered_groups::Vector{String} = lowercase.(string.(functional_group_names()))
    return groups .∈ [considered_groups]
end
