using YAXArrays

# =============================================================================
# C~scape outcome loading
#
# Read a NetCDF variable from every scenario in the result set, aggregate each
# with a caller-supplied function, and assemble one ADRIA-oriented cube. Only
# used by the C~scape metric definitions in metrics/cscape.jl.
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
