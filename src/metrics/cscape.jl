using ADRIA: CScapeResultSet, _load_variable!

"""
    _total_internal_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}

Get the number of internal received larvae across the entire domain. Calculated from the
NetCDF variable "internal_received_larvae".
"""
function _total_internal_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :total_internal_larvae
    if outcome_name in keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Expected dimensions after aggregation excluding dimensions
    out_dims::Tuple = (:year,)
    agg_f = x -> dropdims(sum(x, dims=(2, 3, 4, 5)), dims=(2, 3, 4, 5))

    # name of variable to be used for calculation
    input_var_name::Symbol = :internal_received_larvae
    return _load_variable!(
        rs, input_var_name, out_dims, agg_f, outcome_name; show_progress=show_progress
    )
end
total_internal_larvae = Metric(
    _total_internal_larvae, 
    (:timesteps, :scenarios), 
    "Total Internal Larvae", 
    IS_NOT_RELATIVE,
    "count"
)

"""
    _loc_internal_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}

Get the number of internal received larvae for each location. Calculated from the NetCDF 
variable "internal_received_larvae".
"""
function _loc_internal_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :loc_internal_larvae
    if outcome_name in keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Expected dimensions after aggregation excluding dimensions
    out_dims::Tuple = (:year, :reef_sites)
    agg_f = x -> dropdims(sum(x, dims=(3, 4, 5)), dims=(3, 4, 5))

    # name of variable to be used for calculation
    input_var_name::Symbol = :internal_received_larvae
    return _load_variable!(
        rs, input_var_name, out_dims, agg_f, outcome_name; show_progress=show_progress
    )
end
loc_internal_larvae = Metric(
    _loc_internal_larvae, 
    (:timesteps, :locations, :scenarios), 
    "Location Internal Larvae", 
    IS_NOT_RELATIVE,
    "count"
)

"""
    _total_external_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}

Get the number of external larvae used by cscape across the entire domain. Calculated from
the NetCDF "external_larvae".
"""
function _total_external_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :total_external_larvae
    if outcome_name in keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Expected dimensions after aggregation excluding dimensions
    out_dims::Tuple = (:year,)
    agg_f = x -> dropdims(sum(x, dims=(2, 3, 4, 5)), dims=(2, 3, 4, 5))

    # name of variable to be used for calculation
    input_var_name::Symbol = :external_larvae
    return _load_variable!(
        rs, input_var_name, out_dims, agg_f, outcome_name; show_progress=show_progress
    )
end
total_external_larvae = Metric(
    _total_external_larvae, 
    (:timesteps, :scenarios), 
    "Total External Larvae", 
    IS_NOT_RELATIVE,
    "count"
)

"""
    _loc_external_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}

Get the number of external larvae used in the cscape model for each location. Calculated
from the NetCDF variable "external_larvae".
"""
function _loc_external_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :loc_external_larvae
    if outcome_name in keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Expected dimensions after aggregation excluding dimensions
    out_dims::Tuple = (:year, :reef_sites)
    agg_f = x -> dropdims(sum(x, dims=(3, 4, 5)), dims=(3, 4, 5))

    # name of variable to be used for calculation
    input_var_name::Symbol = :external_larvae
    return _load_variable!(
        rs, input_var_name, out_dims, agg_f, outcome_name; show_progress=show_progress
    )
end
loc_external_larvae = Metric(
    _loc_external_larvae, 
    (:timesteps, :locations, :scenarios), 
    "Location External Larvae", 
    IS_NOT_RELATIVE,
    "count"
)

"""
    _total_eggs_produced(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}

Get the number of eggs produced across the entire domain. Calculated from the NetCDF variable
"eggs".
"""
function _total_eggs_produced(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :total_eggs_produced
    if outcome_name in keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Expected dimensions after aggregation excluding dimensions
    out_dims::Tuple = (:year,)
    agg_f = x -> dropdims(sum(x, dims=(2, 3, 4, 5)), dims=(2, 3, 4, 5))

    # name of variable to be used for calculation
    input_var_name::Symbol = :eggs
    return _load_variable!(
        rs, input_var_name, out_dims, agg_f, outcome_name; show_progress=show_progress
    )
end
total_eggs_produced = Metric(
    _total_eggs_produced, 
    (:timesteps, :scenarios), 
    "Total Eggs Produced", 
    IS_NOT_RELATIVE,
    "count"
)

"""
    _loc_eggs_produced(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}

Get the number of eggs produced at each location of the domain. Calculated from the NetCDF 
variable "eggs".
"""
function _loc_eggs_produced(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :loc_eggs_produced
    if outcome_name in keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Expected dimensions after aggregation excluding dimensions
    out_dims::Tuple = (:year, :reef_sites)
    agg_f = x -> dropdims(sum(x, dims=(3, 4, 5)), dims=(3, 4, 5))

    # name of variable to be used for calculation
    input_var_name::Symbol = :eggs
    return _load_variable!(
        rs, input_var_name, out_dims, agg_f, outcome_name; show_progress=show_progress
    )
end
loc_eggs_produced = Metric(
    _loc_eggs_produced, 
    (:timesteps, :locations, :scenarios), 
    "Location Eggs Produced", 
    IS_NOT_RELATIVE,
    "count"
)

"""
    _total_settlers(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}

Get the number of coral settlers across the entire domain. Calculated from 
NetCDF variable "settlers".
"""
function _total_settlers(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :total_settlers
    if outcome_name in keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Expected dimensions after aggregation excluding dimensions
    out_dims::Tuple = (:year,)
    agg_f = x -> dropdims(sum(x, dims=(2, 3, 4, 5)), dims=(2, 3, 4, 5))

    # name of variable to be used for calculation
    input_var_name::Symbol = :settlers
    return _load_variable!(
        rs, input_var_name, out_dims, agg_f, outcome_name; show_progress=show_progress
    )
end
total_settlers = Metric(
    _total_settlers, 
    (:timesteps, :scenarios), 
    "Total Settlers", 
    IS_NOT_RELATIVE,
    "count"
)

"""
    _loc_settlers(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}

Get the number of coral settlers at each location. Calculated from NetCDF variable 
"settlers".
"""
function _loc_settlers(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :loc_settlers
    if outcome_name in keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Expected dimensions after aggregation excluding dimensions
    out_dims::Tuple = (:year, :reef_sites)
    agg_f = x -> dropdims(sum(x, dims=(3, 4, 5)), dims=(3, 4, 5))

    # name of variable to be used for calculation
    input_var_name::Symbol = :settlers
    return _load_variable!(
        rs, input_var_name, out_dims, agg_f, outcome_name; show_progress=show_progress
    )
end
loc_settlers = Metric(
    _loc_settlers, 
    (:timesteps, :locations, :scenarios), 
    "Location Settlers", 
    IS_NOT_RELATIVE,
    "count"
)

function _relative_loc_taxa_cover(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :relative_loc_taxa_cover
    if outcome_name in keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Expected dimensions after aggregation excluding dimensions
    out_dims::Tuple = (:year, :reef_sites, :ft)
    agg_f = x -> dropdims(sum(x, dims=(4, 5)), dims=(4, 5))

    # name of variable to be used for calculation
    input_var_name::Symbol = :cover
    return _load_variable!(
        rs, input_var_name, out_dims, agg_f, outcome_name;
        use_combined_cover=true, show_progress=show_progress
    )
end

function _relative_taxa_cover(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :relative_taxa_cover
    if outcome_name in keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    _site_k_area = reshape(site_k_area(rs), (1, :, 1, 1))
    loc_taxa_cover::YAXArray = relative_loc_taxa_cover(rs; show_progress=show_progress)
    taxa_cover = dropdims(sum(
        loc_taxa_cover .* _site_k_area, dims=:locations
    ) ./ sum(_site_k_area), dims=:locations)
    rs.outcomes[outcome_name] = taxa_cover

    return taxa_cover
end

function _relative_cover(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :relative_cover
    if outcome_name in keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Clarify to the user why relative species cover is being calculated to prevent
    # confusion.
    @info "Calculating relative species cover for relative cover."
    rel_taxa_loc_cover = relative_loc_taxa_cover(rs; show_progress=show_progress)
    rel_cover = dropdims(sum(rel_taxa_loc_cover, dims=:species), dims=:species)
    rs.outcomes[outcome_name] = rel_cover

    return rs.outcomes[outcome_name]
end
