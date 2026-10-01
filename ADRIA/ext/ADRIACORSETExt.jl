#=
ADRIACORSETExt -- weak-dependency extension (CscapeIntegration branch).

Provides the `:corset` branch of `ADRIA.run_scenarios(dom, scens, RCP, engine; ...)` (declared as a
zero-method stub, `_run_scenarios_corset`, in `src/scenario.jl`). Only active when `CORSETInterface`
is loaded (`using CORSETInterface`), the same pattern as `ADRIAKoraExt`/`ADRIACoralPopExt`.

CORSET's domain translation needs several inputs an ADRIA `Domain` has no source for at all (consumer
biomass above all -- see `CORSETInterface.jl/src/domain_translator.jl`'s file header), so unlike the
other three engines' dispatch targets, this one does NOT translate the domain itself; it forwards
`dom` and the required extra keywords straight to `CORSETInterface.domain_to_corset_inputs`, which
raises a clear `ArgumentError` naming whichever required keyword (`initial_herbivore_biomass`, etc.)
is missing, rather than this extension silently defaulting it.
=#
module ADRIACORSETExt

using ADRIA
import ADRIA: _run_scenarios_corset, Domain
using CORSETInterface
using DataFrames: DataFrame

"""
    ADRIA._run_scenarios_corset(dom, scens, fpath; dhw_idx=0, wave_idx=0, cyclone_idx=0,
                                base_scenario_id=1, n_workers=0, keep_raw=false, spec=nothing,
                                translator_kwargs=NamedTuple(), kwargs...)

Translate `dom` into per-site CORSET inputs (forwarding `translator_kwargs` --
`initial_herbivore_biomass`, `initial_small_piscivore_biomass`, `initial_large_piscivore_biomass`,
`initial_urchin_biomass` are REQUIRED among them, see `domain_translator.jl` -- plus any of
`brooder_fraction`/`turf_fraction`/`algae_fraction`/`other_cover_fraction`/`fish_conn`/`urchin_conn`/
`use_cyclones_as_hurricanes`/`use_dhw_as_bleaching`), run every row of the sampled `corset_spec` table
`scens` through CORSET, and load the result back as a `CORSETInterface.CorsetResultSet`.

Mirrors the manual pipeline in `CORSETInterface.jl`'s README: `domain_to_corset_inputs` ->
`corset_sample_to_scenarios` -> `run_corset_scenarios` -> `load_corset_results`, invoked through
`ADRIA.run_scenarios`.

Keywords
- `dhw_idx`/`cyclone_idx`: forwarded to the translator as the member used for the optional
  DHW-as-bleaching / cyclone-as-hurricane substitution (off by default; see `translator_kwargs`'s
  `use_dhw_as_bleaching`/`use_cyclones_as_hurricanes`). `wave_idx` is accepted for a uniform signature
  and unused (CORSET has no wave-mortality pathway distinct from its hurricane forcing).
- `translator_kwargs`: a `NamedTuple` forwarded to `domain_to_corset_inputs` -- this is where the
  REQUIRED consumer-biomass keywords and any benthic-composition/connectivity overrides go.
- `base_scenario_id`/`n_workers`/`keep_raw`/`spec`: same meanings as `ADRIAKoraExt`'s.

A failing scenario is reported with a warning and is absent from the returned result set.
"""
function ADRIA._run_scenarios_corset(
    dom::Domain, scens::DataFrame, fpath::String;
    dhw_idx::Int=0, wave_idx::Int=0, cyclone_idx::Int=0,
    base_scenario_id::Int=1, n_workers::Int=0, keep_raw::Bool=false,
    spec::Union{DataFrame,Nothing}=nothing, translator_kwargs::NamedTuple=NamedTuple()
)
    inp = CORSETInterface.domain_to_corset_inputs(dom; dhw_idx, cyclone_idx, translator_kwargs...)
    scenarios, ids, _ = CORSETInterface.corset_sample_to_scenarios(scens, base_scenario_id)

    _, failures = CORSETInterface.run_corset_scenarios(
        inp, scenarios, fpath; n_workers, calc_indicators=true, keep_raw
    )
    isempty(failures) || @warn "Some CORSET scenarios failed" ids = [f.id for f in failures]

    return CORSETInterface.load_corset_results(
        fpath; scenarios, dom, spec, scenario_ids=ids, show_progress=false
    )
end

end
