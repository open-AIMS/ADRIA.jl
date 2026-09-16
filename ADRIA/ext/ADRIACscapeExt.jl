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
import ADRIA: _run_scenarios_cscape, Domain
using CscapeInterface
using DataFrames: DataFrame

"""
    ADRIA._run_scenarios_cscape(dom, scens, fpath, base_scenario_id;
                                 use_domain_translation=false, dhw_idx=0, wave_idx=0,
                                 cyclone_idx=0, kwargs...)

Translate a sampled `cscape_spec` table into runnable scenario rows, run them through native
`CScape.jl` in parallel, and load the result back into a `CscapeInterface.CScapeResultSet`.
Mirrors `CscapeInterface.jl/sandbox/ADRIAIntegration/02_PHASE_B.md`'s manual Phase B
pipeline (`cscape_sample_to_scenarios` → `run_cscape_native_parallel` →
`load_results(CScapeResultSet, ...)`), just invoked through `ADRIA.run_scenarios` instead of
by hand.

`use_domain_translation=true` additionally builds `domain_inputs` (via
`CscapeInterface.domain_to_cscape_inputs(dom; dhw_idx, wave_idx, cyclone_idx)`) and a
`demography_override` (via `CScape.build_demography`, reusing the translated spatial table
and `fpath`'s own real `Growth_Surv_file`s for `base_scenario_id`), and passes both through
to every scenario in the batch -- so the whole sweep runs against `dom`'s own spatial/
connectivity/disturbance/cyclone data instead of `fpath`'s `data/` files for those, while
still using `fpath`'s real demography files (ADRIA has no equivalent). `false` (the default)
is unchanged from before this capability was added.

`kwargs` are forwarded to `run_cscape_native_parallel` (e.g. `n_workers`, `export_adria`,
`calc_indicators`).
"""
function ADRIA._run_scenarios_cscape(
    dom::Domain, scens::DataFrame, fpath::String, base_scenario_id::Int;
    use_domain_translation::Bool=false, dhw_idx::Int=0, wave_idx::Int=0, cyclone_idx::Int=0,
    kwargs...
)
    scenarios, ids, draws = CscapeInterface.cscape_sample_to_scenarios(
        fpath, base_scenario_id, scens
    )

    extra_kwargs = if use_domain_translation
        base_row = CscapeInterface.CScape.load_scenario_table(fpath)
        row = base_row[base_row.ID .== base_scenario_id, :]
        fts = String.(split(row.fts[1], "/"))
        demog_files = replace.(
            String.(split(row.Growth_Surv_file[1], "/")), "demog_inputs_" => "demog_julia_"
        )
        domain_inputs = CscapeInterface.domain_to_cscape_inputs(dom; dhw_idx, wave_idx, cyclone_idx)
        demography = CscapeInterface.CScape.build_demography(
            fts, demog_files, fpath, String(row.Region[1]), domain_inputs.spatial_override
        )
        (; domain_inputs, demography_override=demography)
    else
        (;)
    end

    failed = CscapeInterface.run_cscape_native_parallel(
        ids, fpath; scenarios=scenarios, draws=draws, extra_kwargs..., kwargs...
    )
    isempty(failed) || @warn "Some C~scape scenarios failed" failed
    return CscapeInterface.load_results(
        CscapeInterface.CScapeResultSet, fpath; scenarios=scenarios
    )
end

end
