"""Scenario grouping helper functions used by AnnotatedOutcomes."""

"""
    _CAq_off_grp(scenarios::DataFrame)::BitVector

Identify scenarios where no corals were deployed via coral aquaculture (all `N_CAq_*` columns are zero).
"""
function _CAq_off_grp(scenarios::DataFrame)::BitVector
    return dropdims(
        sum(Matrix(scenarios[:, contains.(names(scenarios), "iv_CAq_N")]); dims=2); dims=2
    ) .== 0
end

"""
    _counterfactual_grp(scenarios::DataFrame)::BitVector

Identify counterfactual scenarios: no coral aquaculture, no fogging, no shading, no larval
methods deployment, and no marine cloud brightening.
"""
function _counterfactual_grp(scenarios::DataFrame)::BitVector
    iv_CAq_off = _CAq_off_grp(scenarios)
    no_Fog = scenarios.iv_Fog .== 0
    no_Shd = scenarios.iv_Shd .== 0
    no_LvM = scenarios.iv_LvM_N_settlers .== 0
    no_MCB = scenarios.iv_MCB_duration .== 0
    return iv_CAq_off .& no_Fog .& no_Shd .& no_LvM .& no_MCB
end

"""
    _unguided_grp(scenarios::DataFrame)::BitVector

Identify unguided intervention scenarios: at least one intervention is active but
`guided == 0` (i.e., site selection is random rather than MCDA-driven).
"""
function _unguided_grp(scenarios::DataFrame)::BitVector
    iv_CAq_on = .!_CAq_off_grp(scenarios)
    has_Shd = (scenarios.iv_Fog .> 0) .| (scenarios.iv_Shd .> 0)
    has_LvM = scenarios.iv_LvM_N_settlers .> 0
    has_MCB = scenarios.iv_MCB_duration .> 0
    return (scenarios.guided .== 0) .& (iv_CAq_on .| has_Shd .| has_LvM .| has_MCB)
end

"""
    _guided_grp(scenarios::DataFrame)::BitVector

Identify guided intervention scenarios: all scenarios that are neither counterfactual
nor unguided.
"""
function _guided_grp(scenarios::DataFrame)::BitVector
    return .!(_counterfactual_grp(scenarios) .| _unguided_grp(scenarios))
end

const _SCENARIO_TYPE_FNS = (
    counterfactual=_counterfactual_grp,
    unguided=_unguided_grp,
    guided=_guided_grp
)
const _SCENARIO_TYPE_KEYS = [:counterfactual, :unguided, :guided]

"""
    scenario_rcps(scenarios::DataFrame)::Dict{Symbol,BitVector}

Return a `Dict` mapping each RCP label (e.g. `:RCP45`) to a `BitVector` identifying
which rows of `scenarios` belong to that RCP.
"""
function scenario_rcps(scenarios::DataFrame)::Dict{Symbol,BitVector}
    rcps::Vector{Symbol} = Symbol.(:RCP, Int64.(scenarios[:, :RCP]))
    return Dict(rcp => rcps .== rcp for rcp in unique(rcps))
end

"""
    scenario_types(scenarios::DataFrame)::Dict{Symbol,BitVector}

Return a `Dict` mapping each non-empty scenario type (`:counterfactual`, `:unguided`,
`:guided`) to a `BitVector` identifying which rows of `scenarios` belong to that type.
Types with no matching scenarios are omitted from the result.
"""
function scenario_types(scenarios::DataFrame)::Dict{Symbol,BitVector}
    return Dict(
        type => _SCENARIO_TYPE_FNS[type](scenarios) for
        type in _SCENARIO_TYPE_KEYS if count(_SCENARIO_TYPE_FNS[type](scenarios)) != 0
    )
end
