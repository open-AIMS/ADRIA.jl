using ADRIA: CScapeResultSet, _load_variable!

# =============================================================================
# Larvae / reproduction / recruitment totals
#
# Each raw NetCDF variable below is reduced to either a domain total (`total_*`,
# summed over every dimension except year) or a per-location series (`loc_*`,
# summed over every dimension except year and location). Only the source variable
# and cached outcome name differ between metrics, so the aggregation lives in the
# two helpers here.
# =============================================================================

"""
    _cscape_domain_sum(rs::CScapeResultSet, var::Symbol, out_name::Symbol; show_progress=true)::YAXArray{<:Real}

Domain-wide sum of NetCDF variable `var` over every dimension except year.
"""
function _cscape_domain_sum(
    rs::CScapeResultSet, var::Symbol, out_name::Symbol; show_progress=true
)::YAXArray{<:Real}
    haskey(rs.outcomes, out_name) && return rs.outcomes[out_name]
    agg_f = x -> dropdims(sum(x; dims=(2, 3, 4, 5)); dims=(2, 3, 4, 5))
    return _load_variable!(rs, var, (:year,), agg_f, out_name; show_progress=show_progress)
end

"""
    _cscape_loc_sum(rs::CScapeResultSet, var::Symbol, out_name::Symbol; show_progress=true)::YAXArray{<:Real}

Per-location sum of NetCDF variable `var` over every dimension except year and location.
"""
function _cscape_loc_sum(
    rs::CScapeResultSet, var::Symbol, out_name::Symbol; show_progress=true
)::YAXArray{<:Real}
    haskey(rs.outcomes, out_name) && return rs.outcomes[out_name]
    agg_f = x -> dropdims(sum(x; dims=(3, 4, 5)); dims=(3, 4, 5))
    return _load_variable!(
        rs, var, (:year, :reef_sites), agg_f, out_name; show_progress=show_progress
    )
end

_cscape_total_metric(f, feature) = Metric(
    f, (:timesteps, :scenarios), (:timesteps, :scenarios), feature, IS_NOT_RELATIVE,
    "count"
)
_cscape_loc_metric(f, feature) = Metric(
    f,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :locations, :scenarios),
    feature,
    IS_NOT_RELATIVE,
    "count"
)

"""
    _total_internal_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    _loc_internal_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    _total_external_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    _loc_external_larvae(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    _total_eggs_produced(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    _loc_eggs_produced(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    _total_settlers(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    _loc_settlers(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}

Domain-total and per-location sums of the C~scape larval (`internal_received_larvae`,
`external_larvae`), reproductive (`eggs`) and recruitment (`settlers`) NetCDF variables.
"""
_total_internal_larvae(rs::CScapeResultSet; show_progress=true) = _cscape_domain_sum(
    rs, :internal_received_larvae, :total_internal_larvae; show_progress=show_progress
)
_loc_internal_larvae(rs::CScapeResultSet; show_progress=true) = _cscape_loc_sum(
    rs, :internal_received_larvae, :loc_internal_larvae; show_progress=show_progress
)
_total_external_larvae(rs::CScapeResultSet; show_progress=true) = _cscape_domain_sum(
    rs, :external_larvae, :total_external_larvae; show_progress=show_progress
)
_loc_external_larvae(rs::CScapeResultSet; show_progress=true) = _cscape_loc_sum(
    rs, :external_larvae, :loc_external_larvae; show_progress=show_progress
)
_total_eggs_produced(rs::CScapeResultSet; show_progress=true) = _cscape_domain_sum(
    rs, :eggs, :total_eggs_produced; show_progress=show_progress
)
_loc_eggs_produced(rs::CScapeResultSet; show_progress=true) = _cscape_loc_sum(
    rs, :eggs, :loc_eggs_produced; show_progress=show_progress
)
_total_settlers(rs::CScapeResultSet; show_progress=true) = _cscape_domain_sum(
    rs, :settlers, :total_settlers; show_progress=show_progress
)
_loc_settlers(rs::CScapeResultSet; show_progress=true) = _cscape_loc_sum(
    rs, :settlers, :loc_settlers; show_progress=show_progress
)

total_internal_larvae = _cscape_total_metric(
    _total_internal_larvae, "Total Internal Larvae"
)
loc_internal_larvae = _cscape_loc_metric(_loc_internal_larvae, "Location Internal Larvae")
total_external_larvae = _cscape_total_metric(
    _total_external_larvae, "Total External Larvae"
)
loc_external_larvae = _cscape_loc_metric(_loc_external_larvae, "Location External Larvae")
total_eggs_produced = _cscape_total_metric(_total_eggs_produced, "Total Eggs Produced")
loc_eggs_produced = _cscape_loc_metric(_loc_eggs_produced, "Location Eggs Produced")
total_settlers = _cscape_total_metric(_total_settlers, "Total Settlers")
loc_settlers = _cscape_loc_metric(_loc_settlers, "Location Settlers")

# =============================================================================
# Coral cover
# =============================================================================

function _relative_loc_taxa_cover(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :relative_loc_taxa_cover
    if outcome_name ∈ keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Expected dimensions after aggregation excluding dimensions
    out_dims::Tuple = (:year, :reef_sites, :ft)
    agg_f = x -> dropdims(sum(x; dims=(4, 5)); dims=(4, 5))

    # name of variable to be used for calculation
    input_var_name::Symbol = :cover
    return _load_variable!(
        rs, input_var_name, out_dims, agg_f, outcome_name;
        use_combined_cover=true, show_progress=show_progress
    )
end

function _relative_taxa_cover(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :relative_taxa_cover
    if outcome_name ∈ keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    _site_k_area = reshape(loc_k_area(rs), (1, 1, :, 1))
    loc_taxa_cover::YAXArray = relative_loc_taxa_cover(rs; show_progress=show_progress)
    taxa_cover = dropdims(
        sum(
            loc_taxa_cover .* _site_k_area; dims=:locations
        ) ./ sum(_site_k_area); dims=:locations)
    rs.outcomes[outcome_name] = taxa_cover

    return taxa_cover
end

function _relative_cover(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :relative_cover
    if outcome_name ∈ keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Clarify to the user why relative taxa cover is being calculated to prevent
    # confusion.
    @info "Calculating relative taxa cover for relative cover."
    rel_taxa_loc_cover = relative_loc_taxa_cover(rs; show_progress=show_progress)
    rel_cover = dropdims(sum(rel_taxa_loc_cover; dims=:groups); dims=:groups)
    rs.outcomes[outcome_name] = rel_cover

    return rs.outcomes[outcome_name]
end

function _coral_evenness(rs::CScapeResultSet; show_progress=true)::YAXArray{<:Real}
    outcome_name::Symbol = :coral_evenness
    if outcome_name ∈ keys(rs.outcomes)
        return rs.outcomes[outcome_name]
    end

    # Warn user about different metric calculation, otherwise progress bar may be mistaken
    # for a bug.
    if :relative_taxa_cover ∉ keys(rs.outcomes)
        @info "Calculating relative location species cover for coral evenness index."
    end
    loc_taxa_cover = relative_loc_taxa_cover(rs; show_progress=show_progress)

    _cor_evenness = YAXArray(
        (loc_taxa_cover.timesteps, loc_taxa_cover.locations, loc_taxa_cover.scenarios),
        zeros(
            Float64,
            length(loc_taxa_cover.timesteps),
            length(loc_taxa_cover.locations),
            length(loc_taxa_cover.scenarios)
        ),
        Dict{Symbol,Any}()
    )

    for scen_idx in 1:size(_cor_evenness, :scenarios)
        _cor_evenness[:, :, scen_idx] .= coral_evenness(
            @view(loc_taxa_cover[:, :, :, scen_idx])
        )
    end
    return _cor_evenness
end
