"""

Core ADRIA domain. Represents study area.
"""
struct Domain{M,I,D,S,V,T,X}
    # Matrix{Float64, 2}, Vector{Int}, DataFrame, String, Vector{Float64}, Vector{String}, Matrix{Float64, 3}

    name::S           # human-readable name
    scenario_invoke_time::S  # time latest set of scenarios were run
    TP_data::D     # site connectivity data
    site_ranks::V  # site rank
    strongpred::I  # strongest predecessor
    site_data::D   # table of site data (depth, carrying capacity, etc)
    site_id_col::S  # column to use as site ids, also used by the connectivity dataset (indicates order of `TP_data`)
    unique_site_id_col::S  # column of unique site ids
    init_coral_cover::M  # initial coral cover dataset
    coral_growth::CoralGrowth  # coral
    site_ids::T  # Site IDs that are represented (i.e., subset of site_data[:, site_id_col], after missing sites are filtered)
    removed_sites::T  # indices of sites that were removed. Used to align site_data, DHW, connectivity, etc.
    dhw_scens::X  # DHW scenarios
    wave_scens::X # wave scenarios

    # Parameters
    model::Model
    # intervention::Intervention
    # criteria::Criteria
    # coral::Coral
    sim_constants::SimConstants
end


"""
Barrier function to create Domain struct without specifying Intervention/Criteria/Coral/SimConstant parameters.
"""
function Domain(name, TP_base, site_ranks, strongest_predecessor,
    site_data, site_id_col, unique_site_id_col, init_coral_cover, coral_growth,
    site_ids, removed_sites, DHWs, waves)::Domain

    # intervention = Intervention()
    # criteria = Criteria()
    # coral = Coral()
    model = Model((Intervention(), Criteria(), Coral()))
    sim_constants = SimConstants()
    return Domain(name, "", TP_base, site_ranks, strongest_predecessor, site_data, site_id_col, unique_site_id_col,
        init_coral_cover, coral_growth, site_ids, removed_sites, DHWs, waves,
        model, sim_constants)
end


"""
    Domain(site_data_fn::String, site_id_col::String, unique_site_id_col::String, init_coral_fn::String, conn_path::String, dhw_fn::String, wave_fn::String)

Convenience constructor for Domain
"""
function Domain(name::String, site_data_fn::String, site_id_col::String, unique_site_id_col::String, init_coral_fn::String,
    conn_path::String, dhw_fn::String, wave_fn::String)::Domain

    site_data = GeoDataFrames.read(site_data_fn)

    # Sort data to maintain consistent order
    sort!(site_data, [Symbol(site_id_col)])

    site_data.row_id = 1:nrow(site_data)
    site_data._siteref_id = groupindices(groupby(site_data, Symbol(site_id_col)))

    # tmp_site_ids = site_data[:, [:_siteref_id, Symbol(site_id_col)]]
    conn_ids = site_data[:, site_id_col]
    site_conn = site_connectivity(conn_path, conn_ids, site_data[:, unique_site_id_col], site_data._siteref_id)
    conns = connectivity_strength(site_conn.TP_base)

    # Filter out missing entries
    site_data = site_data[coalesce.(in.(conn_ids, [site_conn.site_ids]), false), :]

    coral_growth = CoralGrowth(nrow(site_data))
    n_sites = coral_growth.n_sites

    # conn_site_names = names(site_conn.TP_base)

    # TODO: Load DHW/Wave data from
    # TODO: Clean these repetitive lines up
    if endswith(dhw_fn, ".mat")
        dhw = matread(dhw_fn)["dhw"]

        if size(dhw, 2) != n_sites
            @warn "Mismatch in DHW data. Truncating so that data size matches!"
            dhw = dhw[:, 1:n_sites, :]
        end

    else
        @warn "Using empty DHW data"
        dhw = zeros(74, n_sites, 50)
    end

    if endswith(wave_fn, ".mat")
        waves = matread(wave_fn)["wave"]

        if size(waves, 2) != n_sites
            @warn "Mismatch in wave data. Truncating so that data size matches!"
            waves = waves[:, 1:n_sites, :]
        end
    else
        @warn "Using empty wave data"
        waves = zeros(74, n_sites, 50)
    end


    if endswith(init_coral_fn, ".mat")
        coral_cover = matread(init_coral_fn)["covers"]

        if !isempty(site_conn.truncated)
            @warn "Mismatch in coral cover data. Truncating so that data size matches!"
            coral_cover = coral_cover[:, 1:n_sites]
        end
    else
        @warn "Using random initial coral cover"
        coral_cover = rand(coral_growth.n_species, n_sites)
    end

    return Domain(name, site_conn.TP_base, conns.site_ranks, conns.strongest_predecessor,
        site_data, site_id_col, unique_site_id_col, coral_cover, coral_growth,
        site_conn.site_ids, site_conn.truncated, dhw, waves)
end


function Base.getproperty(o::Domain, s::Symbol)
    if s == :unique_site_ids
        return o.site_data[:, o.unique_site_id_col]
    end

    # if s == :coral_params
    #     return to_spec(o.coral)
    # end

    if hasfield(typeof(o), s)
        return getfield(o, s)
    end

    return o.super.s
end


function param_table(d::Domain)::DataFrame
    f_names = collect(d.model[:fieldname])
    vals = collect(d.model[:val])
    p_df = DataFrame(OrderedDict(k => v for (k, v) in zip(f_names, vals)))

    return p_df
end


function update_domain!(d::Domain, params::DataFrameRow)
    p_df = DataFrame(d.model)[:, [:fieldname, :val, :ptype, :bounds]]
    p_df[!, :val] = collect(params)

    to_floor = (p_df.ptype .== "integer") .& .!isinteger.(p_df.val)
    if any(to_floor .> 0)
        v = p_df[to_floor, :val]

        # Floor values, capping to maximum bound value
        v .= min.(floor.(v), [b[2] for b in p_df[to_floor, :bounds]])
        p_df[to_floor, :val] = Int.(v)
    end

    # update with new parameters
    update!(d.model, p_df)
end


"""
Extract parameters for a specific model component
"""
function component_params(m::Model, component::Type)::DataFrame
    df::DataFrame = DataFrame(m)
    return df[df.component.==component, :]
end


"""

# Returns
Matrix : n_reps * sites * 3

last dimension indicates: site_id, seeding rank, shading rank
"""
function site_selection(domain::Domain, criteria::DataFrame, ts::Int, n_reps::Int, alg_ind::Int)
    # Site Data
    site_d = domain.site_data
    sr = domain.site_ranks
    area = site_d.area

    # Weights for connectivity , waves (ww), high cover (whc) and low
    wtwaves = criteria.wave_stress           # weight of wave damage in MCDA
    wtheat = criteria.heat_stress            # weight of heat damage in MCDA
    wtconshade = criteria.shade_connectivity # weight of connectivity for shading in MCDA
    wtconseed = criteria.seed_connectivity   # weight of connectivity for seeding in MCDA
    wthicover = criteria.coral_cover_high    # weight of high coral cover in MCDA (high cover gives preference for seeding corals but high for SRM)
    wtlocover = criteria.coral_cover_low     # weight of low coral cover in MCDA (low cover gives preference for seeding corals but high for SRM)
    wtpredecseed = criteria.seed_priority    # weight for the importance of seeding sites that are predecessors of priority reefs
    wtpredecshade = criteria.shade_priority  # weight for the importance of shading sites that are predecessors of priority reefs
    risktol = criteria.deployed_coral_risk_tol # risk tolerance
    depth_min = criteria.depth_min
    depth_offset = criteria.depth_offset

    # Filter out sites outside of desired depth range
    max_depth = depth_min + depth_offset
    depth_criteria = (site_d.sitedepth .> -max_depth) .& (site_d.sitedepth .< -depth_min)

    depth_priority = collect(1:nrow(site_d))[depth_criteria]

    max_cover = site_d.k / 100.0  # Max coral cover at each site

    n_siteint = domain.sim_constants.nsiteint
    w_scens = domain.wave_scens
    dhw_scen = domain.dhw_scens

    sumcover = sum(domain.init_coral_cover, dims=2)
    sumcover = sumcover / 100.0

    ranks = zeros(n_reps, length(depth_priority), 3)

    for i = 1:n_reps
        # site_id, seeding rank, shading rank
        rankingsin = [depth_priority zeros(length(depth_priority), 1) zeros(length(depth_priority), 1)]
        prefseedsites = zeros(1, n_siteint)
        prefshadesites = zeros(1, n_siteint)
        dhw_step = dhw_scen[ts, :, i]
        heatstressprob = dhw_step

        w_step = w_scens[ts, :, i]
        damprob = w_step

        mcda_vars = DMCDA_vars(
            depth_priority,
            n_siteint,
            domain.sim_constants.prioritysites,
            domain.strongpred,
            sr,  # sr.C1
            damprob,
            heatstressprob,
            sumcover,
            max_cover,
            area,
            risktol,
            wtconseed,
            wtconshade,
            wtwaves,
            wtheat,
            wthicover,
            wtlocover,
            wtpredecseed,
            wtpredecshade
        )

        # dMCDA(d_vars, alg_ind, log_seed, log_shade, prefseedsites, prefshadesites, rankingsin)
        (_, _, _, _, rankings) = dMCDA(mcda_vars, alg_ind, false, false, prefseedsites, prefshadesites, rankingsin)
        ranks[i, :, :] = rankings
    end

    return ranks
end
