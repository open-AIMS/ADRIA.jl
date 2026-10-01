#=
ADRIAKoraExt -- weak-dependency extension (CscapeIntegration branch).

Provides the `:kora` branch of `ADRIA.run_scenarios(dom, scens, RCP, engine; ...)` (declared as a
zero-method stub, `_run_scenarios_kora`, in `src/scenario.jl`). Only active when `KoraInterface` is
loaded (`using KoraInterface`) -- CoralBlox-only users never pull in Kora's dependency tree.

Like `ADRIACoralPopExt`, it depends on the interface package rather than the underlying model: the
domain translator (`domain_to_kora_inputs`), the scenario runner (`run_kora_scenarios`) and the
result-set loader (`load_kora_results`) already exist and are tested there, so this file only wires
them behind `ADRIA.run_scenarios`.

Kora has no project directory of its own (no external demography data, unlike C~scape) and no
PythonCall bridge (unlike CoralPop), so the ADRIA `Domain` is ALWAYS the source of its spatial,
initial-cover and forcing data, and `n_workers` is a `Threads.@threads` concurrency cap rather than a
Distributed worker count -- see `KoraInterface.run_kora_scenarios`.
=#
module ADRIAKoraExt

using ADRIA
import ADRIA: _run_scenarios_kora, Domain
using KoraInterface
using DataFrames: DataFrame

"""
    ADRIA._run_scenarios_kora(dom, scens, fpath; dhw_idx=0, wave_idx=0, cyclone_idx=0,
                              base_scenario_id=1, n_workers=0, keep_raw=false,
                              growth_models=nothing, survival_models=nothing, spec=nothing, kwargs...)

Translate `dom` into per-site Kora inputs, run every row of the sampled `kora_spec` table `scens`
through the Kora engine, and load the result back as a `KoraInterface.KoraResultSet`.

Mirrors the manual pipeline in `KoraInterface.jl/README.md`: `domain_to_kora_inputs` ->
`kora_sample_to_scenarios` -> `run_kora_scenarios` -> `load_kora_results`, just invoked through
`ADRIA.run_scenarios`.

Keywords
- `dhw_idx`: which member of `dom`'s DHW data to use (0 = none); pass the same value to every engine
  being compared. `wave_idx`/`cyclone_idx` are accepted for a uniform signature and ignored (Kora has
  no wave or cyclone mortality pathway).
- `base_scenario_id`: scenario IDs are `base_scenario_id * 100_000 + i`.
- `n_workers`: `0` (default) runs sequentially; `n > 1` runs scenarios across up to that many Julia
  threads via `Threads.@threads` (start Julia with `-t N` for this to do anything -- unlike
  `:coralpop`/`:cscape`, this never spawns separate OS processes).
- `keep_raw`: also keep the per-scenario raw output (default: indicators only).
- `growth_models`/`survival_models`: override Kora's bundled "offshore_north" fits.
- `spec`: the `kora_spec` table `scens` was sampled from, so the result set's `model_spec` carries the
  real factor bounds (optional).
- any other keyword is an error (nothing is silently dropped).

A failing scenario is reported with a warning and is absent from the returned result set.
"""
function ADRIA._run_scenarios_kora(
    dom::Domain, scens::DataFrame, fpath::String;
    dhw_idx::Int=0, wave_idx::Int=0, cyclone_idx::Int=0,
    base_scenario_id::Int=1, n_workers::Int=0, keep_raw::Bool=false,
    growth_models=nothing, survival_models=nothing,
    spec::Union{DataFrame,Nothing}=nothing
)
    inp = KoraInterface.domain_to_kora_inputs(dom; dhw_idx, wave_idx, cyclone_idx)
    scenarios, ids, _ = KoraInterface.kora_sample_to_scenarios(scens, base_scenario_id)

    ## indicators are always computed: the result set is built from them
    _, failures = KoraInterface.run_kora_scenarios(
        inp, scenarios, fpath; growth_models, survival_models, n_workers,
        calc_indicators=true, keep_raw
    )
    isempty(failures) || @warn "Some Kora scenarios failed" ids = [f.id for f in failures]

    return KoraInterface.load_kora_results(
        fpath; scenarios, dom, spec, scenario_ids=ids, show_progress=false
    )
end

end
