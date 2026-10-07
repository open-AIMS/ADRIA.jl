#=
ADRIAReefModExt -- weak-dependency extension (CscapeIntegration branch).

Provides the `:reefmod` branch of `ADRIA.run_scenarios(dom, scens, RCP, engine; ...)` (declared as a
zero-method stub, `_run_scenarios_reefmod`, in `src/scenario.jl`). Only active when `ReefModInterface`
is loaded (`using ReefModInterface`), the same pattern as `ADRIAKoraExt`/`ADRIACoralPopExt`/
`ADRIACORSETExt`.

Unlike those three, ReefMod is real MATLAB (`ymbozec/REEFMOD.7.3_GBR`), run as one `matlab -batch ...`
subprocess per scenario rather than an in-process call -- see `ReefModInterface.jl/src/reefmod_bridge.jl`
and its local `ReefMod-GBR` clone's `LOCAL_CHANGES.md`. This needs a real MATLAB on `PATH` (or
`ENV["REEFMOD_MATLAB_EXE"]`); it is designed to run the same way on a workstation with a MATLAB
license or inside an Apptainer image that bundles Julia + MATLAB + the ReefMod-GBR clone together with
the other five engines.
=#
module ADRIAReefModExt

using ADRIA
import ADRIA: _run_scenarios_reefmod, Domain
using ReefModInterface
using DataFrames: DataFrame

"""
    ADRIA._run_scenarios_reefmod(dom, scens, fpath; dhw_idx=0, wave_idx=0, cyclone_idx=0,
                                 base_scenario_id=1, n_workers=0, keep_raw=true,
                                 calc_indicators=false, matlab_exe=nothing, reefmod_gbr_dir=nothing,
                                 spec=nothing, kwargs...)

Translate `dom` into per-site ReefMod inputs, run every row of the sampled `reefmod_spec` table
`scens` through ReefMod-GBR (one `matlab -batch` subprocess per scenario, up to `n_workers`
concurrently), and load the result back as a `ReefModInterface.ReefModResultSet`.

Mirrors the manual pipeline in `ReefModInterface.jl`'s README: `domain_to_reefmod_inputs` ->
`reefmod_sample_to_scenarios` -> `run_reefmod_scenarios` -> `load_reefmod_results`, invoked through
`ADRIA.run_scenarios`.

Keywords
- `dhw_idx`: which member of `dom`'s DHW data to use (0 = none: ReefMod then gets an all-zero DHW,
  i.e. no bleaching, see `domain_translator.jl`'s `_translate_dhw`); pass the same value every other
  engine being compared uses. `wave_idx`/`cyclone_idx` are accepted for a uniform signature and
  ignored (ReefMod has no ADRIA-forced wave pathway, and its own cyclone pathway is switched off for
  this comparison, not driven by `cyclone_idx` -- see `domain_translator.jl`'s header).
- `base_scenario_id`: scenario IDs are `base_scenario_id * 100_000 + i`.
- `n_workers`: `0` (default) runs scenarios sequentially; `n > 0` runs up to that many `matlab`
  subprocesses concurrently (each its own OS process and MATLAB license checkout -- like
  `:coralpop`'s Distributed workers, unlike `:kora`'s `Threads.@threads` cap).
- `keep_raw`: keep each scenario's raw `.mat` (default `true`; set `false` once `calc_indicators=true`
  is confirmed working, to save disk).
- `calc_indicators`: also compute `reefmod_indicators` right after each run (default `false` --
  `load_reefmod_results` can recompute them from the raw `.mat` files either way, see that package's
  README for which indicators are implemented so far).
- `matlab_exe`/`reefmod_gbr_dir`: override `ENV["REEFMOD_MATLAB_EXE"]`/`ENV["REEFMOD_GBR_DIR"]`
  (defaults: `matlab` on `PATH`; the sibling `ReefMod-GBR` folder).
- `spec`: the `reefmod_spec` table `scens` was sampled from, so the result set's `model_spec` carries
  the real factor bounds (optional).
- any other keyword is an error (nothing is silently dropped).

A failing scenario is reported with a warning and is absent from the returned result set.
"""
function ADRIA._run_scenarios_reefmod(
    dom::Domain, scens::DataFrame, fpath::String;
    dhw_idx::Int=0, wave_idx::Int=0, cyclone_idx::Int=0,
    base_scenario_id::Int=1, n_workers::Int=0, keep_raw::Bool=true, calc_indicators::Bool=false,
    matlab_exe::Union{String,Nothing}=nothing, reefmod_gbr_dir::Union{String,Nothing}=nothing,
    spec::Union{DataFrame,Nothing}=nothing
)
    inp = ReefModInterface.domain_to_reefmod_inputs(dom; dhw_idx)
    scenarios, ids, _ = ReefModInterface.reefmod_sample_to_scenarios(scens, base_scenario_id)

    _, failures = ReefModInterface.run_reefmod_scenarios(
        inp, scenarios, fpath; n_workers, keep_raw, calc_indicators, matlab_exe, reefmod_gbr_dir
    )
    isempty(failures) || @warn "Some ReefMod scenarios failed" ids = [f.id for f in failures]

    return ReefModInterface.load_reefmod_results(
        fpath; scenarios, dom, spec, show_progress=false
    )
end

end
