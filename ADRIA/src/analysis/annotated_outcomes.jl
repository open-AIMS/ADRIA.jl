"""
    AnnotatedOutcomes

A `YAXArray` of simulation outcomes paired with scenario-grouping metadata.

# Fields
- `data`: the raw outcome array (scenarios as one of its dimensions).
- `metadata`: a `Dict{Symbol,Any}` whose recognised keys are
  - `:scenario_type_groups` — `Dict{Symbol,BitVector}` mapping scenario-type
    labels (`:counterfactual`, `:unguided`, `:guided`) to boolean masks over
    the scenario dimension of `data`.
  - `:scenario_rcp_groups` — `Dict{Symbol,BitVector}` mapping RCP labels
    (e.g. `:RCP45`) to boolean masks, or `nothing` when RCP information is
    unavailable (e.g. for `RMEResultSet`).

Use `attach_scenario_metadata` to construct an `AnnotatedOutcomes` from a
`ResultSet` or `RMEResultSet`; construct directly only when building metadata
by hand.
"""
struct AnnotatedOutcomes
    data::YAXArray
    metadata::Dict{Symbol,Any}
end

"""
    attach_scenario_metadata(outcomes::YAXArray, rs::ResultSet)::AnnotatedOutcomes

Construct an `AnnotatedOutcomes` by attaching scenario-type and RCP grouping
metadata derived from `rs.inputs` to `outcomes`.

# Arguments
- `outcomes`: a `YAXArray` whose `scenarios` dimension holds integer scenario
  indices that are valid row indices into `rs.inputs`.
- `rs`: an `ADRIAResultSet` (or any `ResultSet`) whose `inputs` field is a
  `DataFrame` with at least the columns required by `scenario_types` and
  `scenario_rcps`.

# Returns
An `AnnotatedOutcomes` with metadata keys `:scenario_type_groups` and
`:scenario_rcp_groups`, both populated as `Dict{Symbol,BitVector}`.
"""
function attach_scenario_metadata(outcomes::YAXArray, rs::ResultSet)::AnnotatedOutcomes
    scenarios_df = rs.inputs[collect(outcomes.scenarios), :]
    return AnnotatedOutcomes(
        outcomes,
        Dict{Symbol,Any}(
            :scenario_type_groups => scenario_types(scenarios_df),
            :scenario_rcp_groups => scenario_rcps(scenarios_df)
        )
    )
end

"""
    attach_scenario_metadata(outcomes::YAXArray, rs::RMEResultSet)::AnnotatedOutcomes

Construct an `AnnotatedOutcomes` by attaching scenario-type grouping metadata
from `rs.scenario_groups` to `outcomes`.

# Arguments
- `outcomes`: a `YAXArray` whose `scenarios` dimension corresponds to the
  scenarios stored in `rs`.
- `rs`: an `RMEResultSet` whose `scenario_groups` field provides pre-computed
  scenario-type masks.

# Returns
An `AnnotatedOutcomes` with `:scenario_type_groups` populated from
`rs.scenario_groups` and `:scenario_rcp_groups` set to `nothing` (RCP
information is not available for `RMEResultSet`).
"""
function attach_scenario_metadata(outcomes::YAXArray, rs::RMEResultSet)::AnnotatedOutcomes
    return AnnotatedOutcomes(
        outcomes,
        Dict{Symbol,Any}(
            :scenario_type_groups => rs.scenario_groups,
            :scenario_rcp_groups => nothing
        )
    )
end
