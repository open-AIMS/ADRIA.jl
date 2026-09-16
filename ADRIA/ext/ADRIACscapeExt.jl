#=
ADRIACscapeExt -- Phase C weak-dependency extension (CscapeIntegration branch).

Provides the `:cscape` branch of `ADRIA.run_scenarios(dom, scens, RCP, engine; ...)`
(declared as a zero-method stub, `_run_scenarios_cscape`, in `src/scenario.jl`). Only active
when `CscapeInterface` is loaded (`using CscapeInterface`) -- CoralBlox-only users never pull
in `CscapeInterface`'s dependency tree (RCall, WGLMakie, ModelParameters, ...).

Depends on `CscapeInterface` rather than the lower-level `CScape.jl` directly: the
sampled-row translator (`cscape_sample_to_scenarios`), the parallel run driver
(`run_cscape_native_parallel`), and the `CScapeResultSet` loader (`load_results`) were all
already built and verified in Phase A/B against `CscapeInterface`'s API -- reimplementing
that orchestration layer against raw `CScape.jl` here would duplicate it for no benefit. This
is a deliberate deviation from `00_PLAN.md`'s literal wording (which named `CScape.jl` as the
weak dependency); see `03_PHASE_C.md` for the rationale.
=#
module ADRIACscapeExt

using ADRIA
import ADRIA: _run_scenarios_cscape
using CscapeInterface
using DataFrames: DataFrame

"""
    ADRIA._run_scenarios_cscape(scens, fpath, base_scenario_id; kwargs...)

Translate a sampled `cscape_spec` table into runnable scenario rows, run them through native
`CScape.jl` in parallel, and load the result back into a `CscapeInterface.CScapeResultSet`.
Mirrors `CscapeInterface.jl/sandbox/ADRIAIntegration/02_PHASE_B.md`'s manual Phase B
pipeline (`cscape_sample_to_scenarios` → `run_cscape_native_parallel` →
`load_results(CScapeResultSet, ...)`), just invoked through `ADRIA.run_scenarios` instead of
by hand.

`kwargs` are forwarded to `run_cscape_native_parallel` (e.g. `n_workers`, `export_adria`,
`calc_indicators`).
"""
function ADRIA._run_scenarios_cscape(
    scens::DataFrame, fpath::String, base_scenario_id::Int; kwargs...
)
    scenarios, ids, draws = CscapeInterface.cscape_sample_to_scenarios(
        fpath, base_scenario_id, scens
    )
    failed = CscapeInterface.run_cscape_native_parallel(
        ids, fpath; scenarios=scenarios, draws=draws, kwargs...
    )
    isempty(failed) || @warn "Some C~scape scenarios failed" failed
    return CscapeInterface.load_results(
        CscapeInterface.CScapeResultSet, fpath; scenarios=scenarios
    )
end

end
