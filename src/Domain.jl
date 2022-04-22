"""

Core ADRIA domain. Represents study area.
"""
struct Domain{M, I, D, S, V, T, X}
    # Matrix{Float64, 2}, Vector{Int}, DataFrame, String, Vector{Float64}, Vector{String}, Matrix{Float64, 3}

    TP_data::D     # site connectivity data
    site_ranks::V  # site rank
    strongpred::I  # strongest predecessor
    site_data::D   # table of site data (depth, carrying capacity, etc)
    site_id_col::S  # column to use as site ids
    init_coral_cover::M  # initial coral cover dataset
    coral_growth::CoralGrowth  # coral
    connectivity_site_ids::T  # Site IDs as specified by the connectivity dataset (indicates order of `TP_data`)
    removed_sites::T  # indices of sites that were removed. Used to align site_data, DHW, connectivity, etc.
    dhw_scens::X  # DHW scenarios
    wave_scens::X # wave scenarios

    # Parameters
    intervention::Intervention
    criteria::Criteria
    coral::Coral
    sim_constants::SimConstants
end


"""
Barrier function to create Domain struct without specifying Intervention/Criteria/Coral/SimConstant parameters.
"""
function Domain(TP_base, site_ranks, strongest_predecessor, 
                  site_data, site_id_col, init_coral_cover, coral_growth, site_ids, removed_sites, DHWs, waves)
    
    intervention = Intervention()
    criteria = Criteria()
    coral = Coral()
    sim_constants = SimConstants()
    return Domain(TP_base, site_ranks, strongest_predecessor, site_data, site_id_col,
                  init_coral_cover, coral_growth, site_ids, removed_sites, DHWs, waves,
                  intervention, criteria, coral, sim_constants);
end


"""
    Domain(site_data_fn::String, site_id_col::String, init_coral_fn::String, conn_path::String, dhw_fn::String, wave_fn::String)

Convenience constructor for Domain
"""
function Domain(site_data_fn::String, site_id_col::String, init_coral_fn::String, 
                conn_path::String, dhw_fn::String, wave_fn::String)::Domain

    site_data = GeoDataFrames.read(site_data_fn)
    tmp_site_ids = site_data[:, site_id_col]

    site_conn = site_connectivity(conn_path, tmp_site_ids)
    conns = connectivity_strength(site_conn.TP_base)

    # Filter out missing entries
    site_data = site_data[coalesce.(in.(tmp_site_ids, [site_conn.site_ids]), false), :]

    # Sort data to maintain consistent order
    sort!(site_data, [Symbol(site_id_col)]);

    coral_growth = CoralGrowth(nrow(site_data))
    n_sites = coral_growth.n_sites

    # TODO: Load DHW/Wave data from
    # TODO: Clean these repetitive lines up
    if endswith(dhw_fn, ".mat")
        dhw = matread(dhw_fn)["dhw"]

        if size(dhw, 2) != n_sites
            @warn "Mismatch in DHW data. Truncating so that data size matches!"
            dhw = dhw[:, 1:n_sites, :];
        end

    else
        @warn "Using empty DHW data"
        dhw = zeros(74, n_sites, 50)
    end

    if endswith(wave_fn, ".mat")
        waves = matread(wave_fn)["wave"]

        if size(waves, 2) != n_sites
            @warn "Mismatch in wave data. Truncating so that data size matches!"
            waves = waves[:, 1:n_sites, :];
        end
    else
        @warn "Using empty wave data"
        waves = zeros(74, n_sites, 50)
    end


    if endswith(init_coral_fn, ".mat")
        coral_cover = matread(init_coral_fn)["covers"]

        if !isempty(site_conn.truncated)
            @warn "Mismatch in coral cover data. Truncating so that data size matches!"
            coral_cover = coral_cover[:, 1:n_sites];
        end
    else
        @warn "Using random initial coral cover"
        coral_cover = rand(coral_growth.n_species, n_sites)
    end

    return Domain(site_conn.TP_base, conns.site_ranks, conns.strongest_predecessor, 
                  site_data, site_id_col, coral_cover, coral_growth, site_conn.site_ids, site_conn.truncated, dhw, waves)
end


function Base.getproperty(o::Domain, s::Symbol)
    if s == :coral_params
        return to_spec(o.coral)
    end

    if hasfield(typeof(o), s)
        return getfield(o, s)
    end

    return o.super.s
end
