#=
ADRIACoralPopExt -- weak-dependency extension (CscapeIntegration branch).

Provides the `:coralpop` branch of `ADRIA.run_scenarios(dom, scens, RCP, engine; ...)`
(declared as a zero-method stub, `_run_scenarios_coralpop`, in `src/scenario.jl`). Only active
when `CoralPopInterface` is loaded (`using CoralPopInterface`) -- CoralBlox-only users never pull
in its dependencies (PythonCall and a Python interpreter with numpy/pandas).

Like `ADRIACscapeExt`, it depends on the interface package rather than the underlying model: the
domain translator (`domain_to_coralpop_inputs`), the scenario runner (`run_coralpop_scenarios`) and
the result-set loader (`load_coralpop_results`) already exist and are tested there, so this file
only wires them behind `ADRIA.run_scenarios`.

One difference from `:cscape`: CoralPop has no project directory of its own (no ScenarioID.xlsx,
no data files), so the ADRIA `Domain` is ALWAYS the source of its spatial, initial-cover,
connectivity and forcing data. There is no `use_domain_translation` switch.
=#
module ADRIACoralPopExt

using ADRIA
import ADRIA: _run_scenarios_coralpop, Domain
using CoralPopInterface
using DataFrames: DataFrame

"""
    ADRIA._run_scenarios_coralpop(dom, scens, fpath; dhw_idx=0, wave_idx=0, cyclone_idx=0,
                                  base_scenario_id=1, n_workers=0, exeflags=String[], coupled=true,
                                  draw_params=false, keep_raw=false, spec=nothing, kwargs...)

Translate `dom` into per-site CoralPop inputs, run every row of the sampled `coralpop_spec` table
`scens` through the CoralPop engine (larvae exchanged between sites through `dom.conn`), and load the
result back as a `CoralPopInterface.CoralPopResultSet`.

Mirrors the manual pipeline in `CoralPopInterface.jl/README.md`:
`domain_to_coralpop_inputs` -> `coralpop_sample_to_scenarios` -> `run_coralpop_scenarios` ->
`load_coralpop_results`, just invoked through `ADRIA.run_scenarios`.

Keywords
- `dhw_idx`/`cyclone_idx`: which member of `dom`'s DHW / cyclone-mortality data to use (0 = none);
  pass the same values to every engine being compared. `wave_idx` is accepted for a uniform
  signature and ignored (CoralPop has no wave mortality).
- `base_scenario_id`: scenario IDs are `base_scenario_id * 100_000 + i`.
- `n_workers`/`exeflags`: `0` runs in this process; `n > 0` starts that many Distributed workers
  (each with its own Python interpreter). Workers inherit this process's environment, so
  `JULIA_PYTHONCALL_EXE` etc. must be set before `using CoralPopInterface`.
- `coupled`: exchange larvae between sites (default) or run sites independently.
- `keep_raw`: also keep the ~15 MB raw output per scenario (default: indicators only, ~4 MB).
- `draw_params`: redraw CoralPop's native random parameters from each scenario seed.
- `spec`: the `coralpop_spec` table `scens` was sampled from, so the result set's `model_spec`
  carries the real factor bounds (optional).
- any other keyword is an error (nothing is silently dropped).

A failing scenario is reported with a warning and is absent from the returned result set.
"""
function ADRIA._run_scenarios_coralpop(
    dom::Domain, scens::DataFrame, fpath::String;
    dhw_idx::Int=0, wave_idx::Int=0, cyclone_idx::Int=0,
    base_scenario_id::Int=1, n_workers::Int=0, exeflags::Vector{String}=String[],
    coupled::Bool=true, draw_params::Bool=false, keep_raw::Bool=false,
    spec::Union{DataFrame,Nothing}=nothing
)
    inp = CoralPopInterface.domain_to_coralpop_inputs(dom; dhw_idx, wave_idx, cyclone_idx)
    scenarios, ids, _ = CoralPopInterface.coralpop_sample_to_scenarios(scens, base_scenario_id)

    ## indicators are always computed: the result set is built from them
    _, failures = CoralPopInterface.run_coralpop_scenarios(
        inp, scenarios, fpath; n_workers, exeflags, coupled, draw_params,
        calc_indicators=true, keep_raw
    )
    isempty(failures) || @warn "Some CoralPop scenarios failed" ids = [f.id for f in failures]

    return CoralPopInterface.load_coralpop_results(
        fpath; scenarios, dom, spec, scenario_ids=ids, show_progress=false
    )
end

end
