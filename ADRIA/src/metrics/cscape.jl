using ADRIA: CScapeResultSet

# =============================================================================
# Outcome loading
#
# Read a NetCDF variable from every scenario in the result set, aggregate each
# with a caller-supplied function, and assemble one ADRIA-oriented cube. Only
# used by the metric definitions below.
# =============================================================================

"""
    _n_scenarios(nc_handle::NcFile)::Int64

Number of draws (scenarios) stored in a C~scape NetCDF; 1 when there is no `draws` dimension.
"""
function _n_scenarios(nc_handle::NcFile)::Int64
    if "draws" in keys(nc_handle.dim)
        return nc_handle.dim["draws"].dimlen
    end
    return 1
end

"""
    _combine_intervention_sites(nc_handle::NcFile, scenario_idx::Union{Int64,Nothing}=nothing)::Array{<:Real}

Area-weighted combination of the intervened and non-intervened `cover` slices. Pass
`scenario_idx` to select a single draw when the variable carries a leading `draws` dimension.
"""
function _combine_intervention_sites(
    nc_handle::NcFile, scenario_idx::Union{Int64,Nothing}=nothing
)::Array{<:Real}
    cover::Array = NetCDF.readvar(nc_handle["cover"]) ./ 100
    area::Array{Float64,2} = NetCDF.readvar(nc_handle["area"])

    lead = isnothing(scenario_idx) ? () : (scenario_idx,)
    area_shape = isnothing(scenario_idx) ? (1, :, 1, 1, 1) : (1, 1, :, 1, 1, 1)
    non_int_idx = (lead..., :, :, 1, :, :)
    int_idx = (lead..., :, :, 2, :, :)
    return (
        @view(cover[non_int_idx...]) .* reshape(@view(area[1, :]), area_shape) .+
        @view(cover[int_idx...]) .* reshape(@view(area[2, :]), area_shape)
    ) ./ reshape(sum(area; dims=1), area_shape)
end

"""
    _read_scenario_slice(nc_handle::NcFile, var_name::String, use_combined_cover::Bool; draw::Union{Int64,Nothing}=nothing)

Read one scenario's data for `var_name` from `nc_handle`. `draw` selects a single index of a
leading `draws` dimension; `nothing` reads the whole variable. When `use_combined_cover` is
set, the area-weighted intervention/counterfactual `cover` combination is returned instead.
"""
function _read_scenario_slice(
    nc_handle::NcFile, var_name::String, use_combined_cover::Bool;
    draw::Union{Int64,Nothing}=nothing
)
    use_combined_cover && return _combine_intervention_sites(nc_handle, draw)
    raw = NetCDF.readvar(nc_handle[var_name])
    isnothing(draw) && return raw
    return raw[draw, ntuple(_ -> Colon(), ndims(raw) - 1)...]
end

"""
    _calculate_first_scenario(nc_handle::NcFile, var_name::String, scenario_func, use_combined_cover::Bool)::Array

Apply `scenario_func` to the first scenario of `var_name`, selecting the first draw when the
variable carries a `draws` dimension. Used to determine the shape of the aggregated output.
"""
function _calculate_first_scenario(
    nc_handle::NcFile,
    var_name::String,
    scenario_func,
    use_combined_cover::Bool
)::Array
    draw = _n_scenarios(nc_handle) > 1 ? 1 : nothing
    return scenario_func(
        _read_scenario_slice(nc_handle, var_name, use_combined_cover; draw=draw)
    )
end

"""
    _load_variable!(rs::CScapeResultSet, variable_name::Symbol, out_dims::Tuple, scenario_func, out_name::Symbol; use_combined_cover=false, show_progress=true)::YAXArray

Aggregate NetCDF variable `variable_name` across every dataset in `rs` into one outcome
cube, store it in `rs.outcomes[out_name]` and return it.

`scenario_func` is applied to one scenario's variable — dimensions
`(year, reef_sites, intervened, ft, thermal_tolerance)`, or
`(year, reef_sites, ft, thermal_tolerance)` when `use_combined_cover` collapses the
`intervened` axis — and must reduce it to exactly the axes named in `out_dims`, with `year`
first. For example `out_dims = (:year,)` with
`scenario_func = x -> dropdims(sum(x; dims=(2,3,4,5)); dims=(2,3,4,5))` for a domain total,
or `out_dims = (:year, :reef_sites)` to keep the location axis. The output axes are then
renamed and reordered to ADRIA's convention by `_reformat_cube`.

`use_combined_cover` reads the area-weighted intervention/counterfactual `cover` combination
instead of `variable_name`.
"""
function _load_variable!(
    rs::CScapeResultSet,
    variable_name::Symbol,
    out_dims::Tuple,
    scenario_func,
    out_name::Symbol;
    use_combined_cover=false,
    show_progress=true
)::YAXArray
    var_name_str::String = String(variable_name)

    # Calculate the shape and dimensions of the new array given the aggregation function
    first_comp::Array = _calculate_first_scenario(
        rs.raw_data[1], var_name_str, scenario_func, use_combined_cover
    )

    # Use the shape of the first calculation to preallocate the space for the rest
    non_scenario_dims = Tuple(
        Dim{dim_name}(1:dim_size) for
        (dim_name, dim_size) in zip(out_dims, size(first_comp))
    )
    n_scens::Vector{Int64} = _n_scenarios.(rs.raw_data)
    dims = (
        Dim{:draws}(1:sum(n_scens)),
        non_scenario_dims...
    )
    out_shape = Tuple(length(d) for d in dims)
    output_variable = YAXArray(
        dims,
        zeros(eltype(rs.raw_data[1][var_name_str]), out_shape...),
        Dict{Symbol,Any}()
    )

    cur_indx = 1
    @showprogress desc = "Calculating $(out_name)" enabled = show_progress for (
        n_sc, nc_handle
    ) in zip(
        n_scens, rs.raw_data
    )
        if n_sc == 1
            output_variable[draws=cur_indx] .= scenario_func(
                _read_scenario_slice(nc_handle, var_name_str, use_combined_cover)
            )
        else
            for j in 0:(n_sc - 1)
                output_variable[draws=cur_indx + j] .= scenario_func(
                    _read_scenario_slice(
                        nc_handle, var_name_str, use_combined_cover; draw=j
                    )
                )
            end
        end
        cur_indx += n_sc
    end
    # Reformat cube to use ADRIA expected dimension names
    output_variable = _reformat_cube(CScapeResultSet, output_variable)

    rs.outcomes[out_name] = output_variable
    return output_variable
end

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
    f, (:timesteps, :scenarios), (:timesteps, :scenarios), feature, IS_NOT_RELATIVE, "count"
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

total_internal_larvae = _cscape_total_metric(_total_internal_larvae, "Total Internal Larvae")
loc_internal_larvae = _cscape_loc_metric(_loc_internal_larvae, "Location Internal Larvae")
total_external_larvae = _cscape_total_metric(_total_external_larvae, "Total External Larvae")
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
