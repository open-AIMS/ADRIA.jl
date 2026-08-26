"""
Model factors of [`effective_dhw_at_depth`](@ref), which attenuates surface DHW before it
reaches the thermal tolerance distributions used by `bleaching_mortality!`.

`eff_dhw_base` and `eff_dhw_mix` supply that function's `κ_base` and `mixing_scale`
arguments respectively.

Both are free parameters — they are unconstrained by the Baird et al. (2018) depth-refuge
data cited in `effective_dhw_at_depth` (see that docstring for why). Exposing them as model
factors lets them be calibrated or swept as a scenario axis instead of being pinned to the
defaults below.

Defaults are the calibrated values from `calib_20260811` (`C:/AIMS/calib_results/calib_20260811/params/calibrated_params.nc`),
with bounds set to +/- 20% of each nominal value.
"""

Base.@kwdef struct DepthAttenuation <: EcoModel
    eff_dhw_base::Param = _eff_dhw_factor(:eff_dhw_base)
    eff_dhw_mix::Param = _eff_dhw_factor(:eff_dhw_mix)
end

const EFF_DHW_DEFAULTS = (
    eff_dhw_base=0.6339028960956812,
    eff_dhw_mix=0.7015827131973273
)

const EFF_DHW_BOUNDS = NamedTuple{keys(EFF_DHW_DEFAULTS)}(
    Tuple((v * 0.8, v * 1.2) for v in EFF_DHW_DEFAULTS)
)

const EFF_DHW_NAMES = (
    eff_dhw_base="Effective DHW Attenuation Base",
    eff_dhw_mix="Effective DHW Mixing Scale"
)

const EFF_DHW_DESCRIPTIONS = (
    eff_dhw_base="Base DHW attenuation coefficient (m⁻¹) in the surface DHW → 0 limit \
                  (`κ_base`). Zero removes the depth refuge entirely.",
    eff_dhw_mix="Surface DHW at which mixing halves the effective attenuation coefficient \
                 (`mixing_scale`). Zero represents permanently well-mixed conditions, in \
                 which depth provides no refuge at any DHW."
)

function _eff_dhw_factor(fieldname::Symbol, val::Float64)::Param
    return Factor(
        val;
        ptype="continuous",
        dist=Uniform,
        dist_params=getproperty(EFF_DHW_BOUNDS, fieldname),
        name=getproperty(EFF_DHW_NAMES, fieldname),
        description=getproperty(EFF_DHW_DESCRIPTIONS, fieldname)
    )
end
function _eff_dhw_factor(fieldname::Symbol)::Param
    return _eff_dhw_factor(fieldname, getproperty(EFF_DHW_DEFAULTS, fieldname))
end

"""
    create_depth_attenuation_instance(; overrides::Dict{String,Float64}=Dict{String,Float64}())::DepthAttenuation

Build a `DepthAttenuation` instance, replacing default factor values with any matching entry
in `overrides` (keyed by field name). Bounds and metadata are unaffected.
"""
function create_depth_attenuation_instance(;
    overrides::Dict{String,Float64}=Dict{String,Float64}()
)::DepthAttenuation
    return DepthAttenuation(;
        (
            fn => _eff_dhw_factor(
                fn, get(overrides, string(fn), getproperty(EFF_DHW_DEFAULTS, fn))
            )
            for fn in fieldnames(DepthAttenuation)
        )...
    )
end
