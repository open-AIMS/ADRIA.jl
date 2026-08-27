"""
    _coral_calib_overrides(nc_ds)::Dict{String,Float64}

Read calibrated coral parameter values from a NetCDF dataset (produced by CoralBlox
calibration) and return a Dict mapping Coral struct field names to their calibrated values.
Covers `linear_extension`, `mb_rate`, `dist_mean`, `dist_std`, `linear_extension_scale`,
and `mb_rate_scale`. Any field not present in the dataset falls back to the ADRIA default.

Variables are normally `(n_groups, n_sizes)`. A 1-D `(n_groups,)` variable is accepted and
broadcast across every size class: `dist_std` is calibrated per functional group, matching
ADRIA's own default of `repeat(<5 group values>; inner=n_sizes)`.
"""
function _coral_calib_overrides(nc_ds)::Dict{String,Float64}
    overrides = Dict{String,Float64}()
    fg_names = string.(functional_group_names())
    n_sizes = size(bin_widths(), 2)

    for (param_name, varname) in (
        ("linear_extension", "linear_extension"),
        ("mb_rate", "mb_rate"),
        ("dist_mean", "dist_mean"),
        ("dist_std", "dist_std")
    )
        # Older calibration outputs predate calibrated `dist_std` entirely
        Symbol(varname) in keys(nc_ds.cubes) || continue

        data = Array(nc_ds[varname])
        # Per-group (n_groups,) → broadcast to (n_groups, n_sizes). Files written before
        # dist_std became per-group carry the full matrix and are used as-is.
        vals = ndims(data) == 1 ? repeat(data, 1, n_sizes) : data

        for (fg_idx, fg) in enumerate(fg_names), sc = 1:size(vals, 2)
            overrides["$(fg)_$(fg_idx)_$(sc)_$(param_name)"] = vals[fg_idx, sc]
        end
    end

    le_scale_da = nc_ds["linear_extension_scale"]
    mb_scale_da = nc_ds["mb_rate_scale"]
    bg_ids = collect(DimensionalData.lookup(le_scale_da, :cb_calib_group))
    le_scale = Array(le_scale_da)
    mb_scale = Array(mb_scale_da)

    for (fg_idx, fg) in enumerate(fg_names), (bg_idx, bg) in enumerate(bg_ids)
        overrides["linear_extension_scale_cb_group_$(bg)_$(fg)"] = le_scale[fg_idx, bg_idx]
        overrides["mb_rate_scale_cb_group_$(bg)_$(fg)"] = mb_scale[fg_idx, bg_idx]
    end

    return overrides
end

"""
    _growth_accel_calib_overrides(nc_ds)::Dict{String,Float64}

Read calibrated growth acceleration parameter values from a NetCDF dataset and return a
Dict mapping GrowthAcceleration struct field names to their calibrated values.
"""
function _growth_accel_calib_overrides(nc_ds)::Dict{String,Float64}
    overrides = Dict{String,Float64}()
    ga_da = nc_ds["growth_acceleration"]  # (cb_calib_group=12, accel_param=3)
    bg_ids = collect(DimensionalData.lookup(ga_da, :cb_calib_group))
    ap_vals = collect(DimensionalData.lookup(ga_da, :accel_param))
    ga = Array(ga_da)

    for (bg_idx, bg) in enumerate(bg_ids), (ap_idx, ap) in enumerate(ap_vals)
        overrides["growth_acceleration_cb_group_$(bg)_$(ap)"] = ga[bg_idx, ap_idx]
    end

    return overrides
end

"""
    _depth_attenuation_calib_overrides(calib_dataset::Dataset)::Dict{String,Float64}

Read calibrated depth-attenuation values from a NetCDF dataset and return a Dict mapping
`DepthAttenuation` struct field names to their calibrated values.

Calibration outputs written before `depth_attenuation` entered the search space have no such
variable, in which case an empty Dict is returned and the ADRIA defaults stand. Values are
looked up by their `depth_atten_param` label rather than by position, so the two factors
cannot be silently transposed.

# Arguments
- `calib_dataset` : Calibration parameters NetCDF, opened with `open_dataset` (produced by
    CoralBlox calibration; see `load_calib_params`)
"""
function _depth_attenuation_calib_overrides(calib_dataset::Dataset)::Dict{String,Float64}
    overrides = Dict{String,Float64}()
    :depth_attenuation in keys(calib_dataset.cubes) || return overrides

    da = calib_dataset["depth_attenuation"]
    p_names = string.(collect(DimensionalData.lookup(da, :depth_atten_param)))
    vals = Array(da)

    for fn in string.(fieldnames(DepthAttenuation))
        idx = findfirst(==(fn), p_names)
        isnothing(idx) || (overrides[fn] = vals[idx])
    end

    return overrides
end

"""
    load_calib_params(calib_params_fn::String)

Build the `Coral`, `GrowthAcceleration` and `DepthAttenuation` component instances for a
domain.

If `calib_params_fn` points to an existing CoralBlox calibration NetCDF, its parameter
values are applied as overrides; otherwise the ADRIA defaults are used.

# Returns
Tuple of (Coral, GrowthAcceleration, DepthAttenuation)
"""
function load_calib_params(calib_params_fn::String)
    if isempty(calib_params_fn) || !isfile(calib_params_fn)
        return Coral(), GrowthAcceleration(), DepthAttenuation()
    end

    @info "Loading calibrated coral parameters from $(calib_params_fn)"
    nc_ds = open_dataset(calib_params_fn)

    return (
        create_coral_instance(; overrides=_coral_calib_overrides(nc_ds)),
        create_growth_acceleration_instance(;
            overrides=_growth_accel_calib_overrides(nc_ds)
        ),
        create_depth_attenuation_instance(;
            overrides=_depth_attenuation_calib_overrides(nc_ds)
        )
    )
end
