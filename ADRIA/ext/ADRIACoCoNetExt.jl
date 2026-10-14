#=
ADRIACoCoNetExt -- weak-dependency extension (CscapeIntegration branch).

Provides the `:coconet` branch of `ADRIA.run_scenarios(dom, scens, RCP, engine; ...)` (declared as a
zero-method stub, `_run_scenarios_coconet`, in `src/scenario.jl`). Only active when `CoCoNetInterface`
is loaded (`using CoCoNetInterface`) -- CoralBlox-only users never pull in its dependencies (PythonCall
and a Python interpreter with numpy/pandas/pyyaml, plus the local CoCoNet clone).

Like `ADRIACoralPopExt`, it depends on the interface package, not the model: the domain translator,
the scenario runner and the result-set loader exist and are tested there, so this file only wires them
behind `ADRIA.run_scenarios`. CoCoNet has no project directory of its own, so the ADRIA `Domain` is
always the source of its reefs, initial cover, connectivity and forcing.
=#
module ADRIACoCoNetExt

using ADRIA
import ADRIA: _run_scenarios_coconet, Domain
using CoCoNetInterface
using DataFrames: DataFrame

"""
    ADRIA._run_scenarios_coconet(dom, scens, fpath; dhw_idx=0, wave_idx=0, cyclone_idx=0,
                                 base_scenario_id=1, n_workers=0, exeflags=String[], keep_raw=false,
                                 spec=nothing, translator_kwargs=NamedTuple())

Translate `dom` into a CoCoNet reef network, run every row of the sampled `coconet_spec` table `scens`
through CoCoNet, and load the result back as a `CoCoNetInterface.CoCoNetResultSet`.

The manual pipeline in `CoCoNetInterface.jl/README.md` -- `domain_to_coconet_inputs` ->
`coconet_sample_to_scenarios` -> `run_coconet_scenarios` -> `load_coconet_results` -- invoked through
`ADRIA.run_scenarios`.

Keywords
- `dhw_idx`/`cyclone_idx`: which member of `dom`'s DHW / cyclone-mortality data to use (0 = none); pass
  the same values to every engine being compared. `wave_idx` is accepted for a uniform signature and
  ignored (CoCoNet has no wave pathway).
- `translator_kwargs`: a `NamedTuple` forwarded to `domain_to_coconet_inputs` (`km_offshore`,
  `initial_rubble`, `ocean_acidification`, `cots`, `connectivity_scale`, `max_sites`, ...).
- `base_scenario_id`: scenario IDs are `base_scenario_id * 100_000 + i`.
- `n_workers`/`exeflags`: `0` runs in this process; `n > 0` starts that many Distributed workers (each
  with its own Python interpreter). Workers inherit this process's environment, so
  `JULIA_PYTHONCALL_EXE`/`COCONET_DIR` must be set before `using CoCoNetInterface`.
- `keep_raw`: also keep each scenario's raw output (default: indicators only).
- `spec`: the `coconet_spec` table `scens` was sampled from, for an accurate `model_spec` (optional).
- any other keyword is an error (nothing is silently dropped).

A failing scenario is reported with a warning and is absent from the returned result set.
"""
function ADRIA._run_scenarios_coconet(
    dom::Domain, scens::DataFrame, fpath::String;
    dhw_idx::Int=0, wave_idx::Int=0, cyclone_idx::Int=0,
    base_scenario_id::Int=1, n_workers::Int=0, exeflags::Vector{String}=String[],
    keep_raw::Bool=false, spec::Union{DataFrame,Nothing}=nothing,
    translator_kwargs::NamedTuple=NamedTuple()
)
    inp = CoCoNetInterface.domain_to_coconet_inputs(dom; dhw_idx, wave_idx, cyclone_idx, translator_kwargs...)
    scenarios, ids, _ = CoCoNetInterface.coconet_sample_to_scenarios(scens, base_scenario_id)

    ## indicators are always computed: the result set is built from them
    _, failures = CoCoNetInterface.run_coconet_scenarios(
        inp, scenarios, fpath; n_workers, exeflags, calc_indicators=true, keep_raw
    )
    isempty(failures) || @warn "Some CoCoNet scenarios failed" ids = [f.id for f in failures]

    return CoCoNetInterface.load_coconet_results(
        fpath; scenarios, dom, spec, scenario_ids=ids, show_progress=false
    )
end

end
