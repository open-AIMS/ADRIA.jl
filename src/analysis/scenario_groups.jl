"""Scenario grouping helper functions used by AnnotatedOutcomes."""

"""
    _no_seed_grp(scenarios::DataFrame)::BitVector

Identify scenarios where no corals were seeded (all `N_seed_*` columns are zero).
"""
function _no_seed_grp(scenarios::DataFrame)::BitVector
    return dropdims(
        sum(Matrix(scenarios[:, contains.(names(scenarios), "N_seed")]); dims=2); dims=2
    ) .== 0
end

"""
    _counterfactual_grp(scenarios::DataFrame)::BitVector

Identify counterfactual scenarios: no seeding, no fogging, no SRM, no macro-colonization,
and no marine cloud brightening.
"""
function _counterfactual_grp(scenarios::DataFrame)::BitVector
    no_seed = _no_seed_grp(scenarios)
    no_fog = scenarios.fogging .== 0
    no_SRM = scenarios.SRM .== 0
    no_mc = scenarios.N_mc_settlers .== 0
    no_mcb = scenarios.mcb_duration .== 0
    return no_seed .& no_fog .& no_SRM .& no_mc .& no_mcb
end

"""
    _unguided_grp(scenarios::DataFrame)::BitVector

Identify unguided intervention scenarios: at least one intervention is active but
`guided == 0` (i.e., site selection is random rather than MCDA-driven).
"""
function _unguided_grp(scenarios::DataFrame)::BitVector
    has_seed = .!_no_seed_grp(scenarios)
    has_shade = (scenarios.fogging .> 0) .| (scenarios.SRM .> 0)
    has_mc_corals = scenarios.N_mc_settlers .> 0
    has_mcb = scenarios.mcb_duration .> 0
    return (scenarios.guided .== 0) .& (has_seed .| has_shade .| has_mc_corals .| has_mcb)
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
    guided=_guided_grp,
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
