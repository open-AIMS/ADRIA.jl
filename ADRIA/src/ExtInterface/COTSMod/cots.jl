# ADRIA adapter for the crown-of-thorns starfish (COTS) submodel in COTSMod.jl.
#
# ADRIA owns the parameter table, domain wiring, coral-cover tensors, result logging and
# scenario timing. COTSMod owns the COTS state transition, predation, larval dispersal and
# initialisation. Units follow COTSMod: COTS densities are COTS ha^-1.

"""
    cots_runtime_params(param_set::YAXArray)::COTSMod.COTSParams

Build COTSMod's runtime parameters from a scenario's `CotsParams` factors, following
`COTS_PARAM_MAP`. Mechanisms that are not exposed as factors keep COTSMod's defaults (off).
"""
function cots_runtime_params(param_set::YAXArray)::COTSMod.COTSParams
    kwargs = (
        field => _cots_param_value(factor, param_set[At(String(factor))])
        for (factor, field) in pairs(COTS_PARAM_MAP)
    )
    return COTSMod.COTSParams(; kwargs...)
end

# Convert a factor value to the type COTSMod expects, from the factor's `ptype`
function _cots_param_value(factor::Symbol, val)
    ptype = getproperty(COTS_FACTOR_SPECS, factor).ptype
    ptype == "unordered categorical" && return val == 1
    ptype == "ordered discrete" && return round(Int, val)
    return Float64(val)
end
