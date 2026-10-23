"""
Model factors for an optional reduction of bleaching mortality at low heat stress, applied
in `bleaching_mortality!`.

Where a location's effective (depth-attenuated) DHW is below `low_dhw_bleach_threshold`,
the bleaching mortality fraction is multiplied by `1 - low_dhw_bleach_reduction`.

The defaults (reduction 0, threshold 0 DHW) leave bleaching mortality unchanged, so the
established model is reproduced unless the factors are set. Both defaults are zero so
that adding the factors does not change the per-scenario seed that ADRIA derives from
the sum of parameter values (`set_random_seed`).

The bounds are broad placeholders, not evidence-derived: the reduction is a fraction in
[0, 1] and the threshold spans 0-12 DHW. The factors are fixed (not sampled) by default,
like `DepthAttenuation`, and can be calibrated or swept as a scenario axis, or set from a
calibration NetCDF via `load_calib_params`.
"""
Base.@kwdef struct LowDHWBleaching <: EcoModel
    low_dhw_bleach_reduction::Param = _low_dhw_bleach_factor(:low_dhw_bleach_reduction)
    low_dhw_bleach_threshold::Param = _low_dhw_bleach_factor(:low_dhw_bleach_threshold)
end

const LOW_DHW_BLEACH_DEFAULTS = (
    low_dhw_bleach_reduction=0.0,
    low_dhw_bleach_threshold=0.0
)

const LOW_DHW_BLEACH_BOUNDS = (
    low_dhw_bleach_reduction=(0.0, 1.0),
    low_dhw_bleach_threshold=(0.0, 12.0)
)

const LOW_DHW_BLEACH_NAMES = (
    low_dhw_bleach_reduction="Low-DHW Bleaching Reduction",
    low_dhw_bleach_threshold="Low-DHW Bleaching Threshold"
)

const LOW_DHW_BLEACH_DESCRIPTIONS = (
    low_dhw_bleach_reduction="Fractional reduction of bleaching mortality where effective \
                              DHW is below the low-DHW threshold. Zero leaves bleaching \
                              mortality unchanged.",
    low_dhw_bleach_threshold="Effective DHW (DHW-weeks) below which the low-DHW bleaching \
                              reduction applies. Zero disables the reduction."
)

function _low_dhw_bleach_factor(fieldname::Symbol, val::Float64)::Param
    return Factor(
        val;
        ptype="continuous",
        dist=Uniform,
        dist_params=getproperty(LOW_DHW_BLEACH_BOUNDS, fieldname),
        name=getproperty(LOW_DHW_BLEACH_NAMES, fieldname),
        description=getproperty(LOW_DHW_BLEACH_DESCRIPTIONS, fieldname)
    )
end
function _low_dhw_bleach_factor(fieldname::Symbol)::Param
    return _low_dhw_bleach_factor(fieldname, getproperty(LOW_DHW_BLEACH_DEFAULTS, fieldname))
end

"""
    create_low_dhw_bleaching_instance(; overrides::Dict{String,Float64}=Dict{String,Float64}())::LowDHWBleaching

Build a `LowDHWBleaching` instance, replacing default factor values with any matching entry
in `overrides` (keyed by field name). Bounds and metadata are unaffected.
"""
function create_low_dhw_bleaching_instance(;
    overrides::Dict{String,Float64}=Dict{String,Float64}()
)::LowDHWBleaching
    return LowDHWBleaching(;
        (
            fn => _low_dhw_bleach_factor(
                fn, get(overrides, string(fn), getproperty(LOW_DHW_BLEACH_DEFAULTS, fn))
            )
            for fn in fieldnames(LowDHWBleaching)
        )...
    )
end
