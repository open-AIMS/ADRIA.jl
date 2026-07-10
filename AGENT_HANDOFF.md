# ADRIA COTS Submodel & Lizard Island Calibration: Agent Handoff Document

> [!IMPORTANT]
> **To Future AI Agents (Codex, Antigravity, etc.) and Developers:**
> This document (`AGENT_HANDOFF.md`) is the authoritative **Point of Truth** for the Crown-of-Thorns Starfish (COTS) submodel calibration and integration within `ADRIA.jl`. Please read this entire document carefully before making code edits or running optimization sweeps.

---

## 1. Executive Summary & Objective
Our objective is to calibrate and refine the ADRIA COTS submodel on the **Lizard Island Cluster domain** (`Lizard_Historical_v0.1`, covering 1985–2024) to accurately reproduce historical ~15-year predator-prey boom-bust cycles and empirical observations (`reef_cots.csv` and `reef_manta.csv`).

We have successfully resolved previous lag and flatlining issues by combining:
1. A **4-parameter Latin Hypercube Sweep (LHS)** of core ecological parameters (`a_F`, `a_S`, `IMM`, `seed_mult`).
2. The **Cairns Initiation Box** external larval pulse boundary condition, which injects incoming larvae from upstream northern GBR outbreaks every 15 years starting at `t = 1` (1985).

---

## 2. Directory Structure & Key Files (`sandbox/`)
To maintain clean separation and avoid ad-hoc script sprawl, all COTS calibration and visualization workflows live inside organized subdirectories under `sandbox/`:

```
sandbox/
├── calibration/
│   ├── formal_calibration_sweep.jl    # Sequential/Parallel LHS sweep generating parameter ensembles
│   ├── analyze_calibration.jl         # Post-sweep analysis extracting top-N parameter candidates
│   ├── simulate_best_stochastic.jl    # Runs N stochastic simulations for a single fixed parameter candidate
│   └── save_best_calibration.jl       # Utility to run & export CSVs for single parameter candidates
├── domain_building/
│   ├── build_lizard_domain.jl         # Builds Lizard_Historical_v0.1 domain from RME spatial/connectivity data
│   └── build_historical_dhw.jl        # Prepares NetCDF DHW historical cubes (dhw_RCPhistorical.nc / dhw_RCP45.nc)
├── plotting/
│   ├── plot_ensemble.py               # Python Matplotlib script for multi-candidate or multi-run plotting
│   ├── plot_best_stochastic.py        # Python script plotting median + 10th-90th %ile uncertainty bands
│   ├── plot_calibration.py            # Plots single trajectory validation curves vs empirical data
│   └── plot_disturbances.py           # Plots DHW thermal stress and cyclone disturbance timelines
├── data/
│   ├── Lizard_Historical_v0.1/        # The compiled historical domain data package (spatial, connectivity, DHWs)
│   ├── reef_cots.csv                  # Empirical COTS manta/tow survey observations (1985–2024)
│   ├── reef_manta.csv                 # Empirical coral cover survey observations (1985–2024)
│   ├── formal_calibration_ensemble.csv# Full table of LHS sweep parameter candidates and their loss scores
│   ├── top5_trajectories.csv          # Simulated COTS & coral trajectories for the top 5 candidates
│   └── best_stochastic_trajectories.csv # Simulated trajectories from N stochastic runs of the single best candidate
└── AGENT_HANDOFF.md                   # Reference copy of this handoff document
```

---

## 3. Core Architecture & Integration Points (`src/`)

### A. Relative Cover Translation Layer (`src/scenario.jl` & `src/ecosystem/cots.jl`)
* **How COTS Predation Operates:** COTS predation operates on **relative coral cover (`C_cover_t`)**, not absolute cover directly. In `scenario.jl`, environmental stressors (bleaching/DHW, cyclones) are applied first to `C_cover_t`, and `cots_mortality!(C_cover_t, cots_models, cots_prey_map)` runs *afterward*. This ensures COTS correctly sense and respond to food-limited environments following bleaching or storms.
* **Initial Density Injection:** `cots.jl` reads `ENV["COTS_INITIAL_MULTIPLIER"]` during model initialization to scale the baseline density (`N[1]`) at `t = 1` across all sites.

### B. The Cairns Initiation Box (`src/scenario.jl`)
To model external larval ingress riding the East Australian Current southward from the Cairns initiation zone every ~15 years, `scenario.jl` checks `ENV` variables during every time step:
```julia
# Mimic incoming larvae/recruitment from upstream outbreak (Cairns Initiation Box)
_cots_pulse_enabled = cots_enabled && get(ENV, "COTS_EXTERNAL_PULSE", "false") == "true"
if _cots_pulse_enabled && !isnothing(_cots_seed_locs)
    _pulse_period = parse(Int, get(ENV, "COTS_PULSE_PERIOD", "15"))
    _pulse_offset = parse(Int, get(ENV, "COTS_PULSE_OFFSET", "1"))
    _pulse_val = parse(Float64, get(ENV, "COTS_PULSE_VAL", "2.0"))
    
    if mod(tstep - _pulse_offset, _pulse_period) == 0
        inject_upstream_pulse!(cots_models, _cots_seed_locs; pulse_val=_pulse_val)
    end
end
```
* **Required ENV Setup for Phase 9b & Best Candidate:**
  ```julia
  ENV["COTS_EXTERNAL_PULSE"] = "true"
  ENV["COTS_PULSE_PERIOD"] = "15"
  ENV["COTS_PULSE_OFFSET"] = "1"
  ENV["COTS_PULSE_VAL"] = "1.5"
  ```

---

## 4. Current Best-Fit Parameterization (Phase 9b Results)

Our lowest loss achieved to date is **`Loss = 23.0897`** (Candidate **`run_id = 2`** from the 250-scenario LHS sweep).

| Parameter | Symbol | Best Candidate (`run_id = 2`) Value | LHS Search Bounds | Description |
| :--- | :---: | :---: | :---: | :--- |
| **Fast-Coral Starvation Threshold** | `a_F` | **`1.34378`** | `[0.1, 2.0]` | Critical consumption rate / threshold for Acropora / fast growers |
| **Slow-Coral Starvation Threshold** | `a_S` | **`0.0957831`** | `[0.01, 0.9]` | Consumption rate / threshold for massive / slow growers |
| **Background Immigration Rate** | `IMM` | **`0.0947791`** | `[0.0, 0.1]` | Baseline low-level background immigration across reefs |
| **Initial Seed Multiplier** | `seed_mult` | **`2.23695`** | `[0.5, 3.0]` | Multiplier on starting COTS population density via `COTS_INITIAL_MULTIPLIER` |
| **Cairns Pulse Amplitude** | `pulse_val` | **`1.5`** | Fixed / ENV | Added to `Age-0` (`N[1]`) at `t=1, 16, 31` (1985, 2000, 2015) |

> [!TIP]
> **Why `run_id = 2` works so well:** The relatively high `a_F` (`1.34`) combined with low `a_S` (`0.096`) allows COTS to rapidly consume fast-growing corals during an outbreak and then experience sharp, deep starvation collapses (**busts**) down to ~7% normalized density (`2000–2008`). The external Cairns pulse arriving around `t=31` (`2015`) perfectly ignites the second historical outbreak wave, peaking around `2018–2019`.

---

## 5. Known Gotchas & Troubleshooting Guide

### Current Codex Phase: Stochastic Best-Fit Diagnostics
Codex is actively expanding `sandbox/calibration/simulate_best_stochastic.jl` and
`sandbox/plotting/plot_best_stochastic.py` to separate three different sources of
spread:
1. **Across-run environmental spread** from `ADRIA.sample`.
2. **Within-reef site spread** from site-level COTS/coral trajectories.
3. **Controlled COTS demographic spread** around the best-fit candidate.

Completed so far:
- `simulate_best_stochastic.jl` now writes reef-level trajectories with
  `sim_cots_adult` plus shared reef-level normalization, and also writes
  `sandbox/data/best_stochastic_site_trajectories.csv`.
- `plot_best_stochastic.py` now reads the site-level CSV, plots site p10-p90
  envelopes, fixes the hardcoded run-count label, and uses the same COTS
  normalization scale for site bands and reef medians.
- The external Cairns pulse now has explicit diagnostic seed locations via
  `ENV["ADRIA_DEBUG_SEED_FIRST_N"] = "10"` in the stochastic simulation script.

Observed diagnostic result:
- Across-run stochastic spread remains very small when only environmental/sample
  fields vary.
- Within-reef site spread is substantial, especially for reefs with many sites
  (Lizard Island Reef and Eyrie Reef).

Completed implementation step:
- `simulate_best_stochastic.jl` now supports `COTS_STOCHASTIC_MODE=environmental|demographic`.
- Default mode is `demographic`: `a_F`, `a_S`, and `IMM` stay locked to the best candidate, while other COTS demographic parameters, initial seed multiplier, and Cairns pulse amplitude get controlled lognormal jitter.
- Per-run values are written to `sandbox/data/best_stochastic_metadata.csv`.

In progress:
- The first demographic jitter pass (`COTS_DEMOGRAPHIC_CV=0.15`, seed/pulse CV `0.10`) created visible across-run spread but scattered peak timing too widely.
- Defaults were tightened to `COTS_DEMOGRAPHIC_CV=0.08`, `COTS_SEED_MULT_CV=0.05`, and `COTS_PULSE_VAL_CV=0.05` before regenerating the baseline plot.


Validation metrics completed:
- Pearson correlation, Spearman rank correlation, RMSE, and Percent Bias have been added to the calibration plotting workflow.
- Metrics are calculated on matched observation years only.
- `plot_best_stochastic.py` compares observations against the simulated median trajectory, displays the COTS metrics in each panel subtitle, and writes `sandbox/data/best_stochastic_validation_metrics.csv` for both COTS and coral.
- `plot_calibration.py` compares observations against the single simulated trajectory, displays the COTS metrics in each panel subtitle, and writes `sandbox/data/calibration_validation_metrics.csv` for both COTS and coral.
- Regenerated figures: `sandbox/best_stochastic_plot.png` and `sandbox/calibration_plot.png`.
- Current stochastic COTS metrics are weak on several reefs, so use the new metrics as calibration diagnostics rather than as evidence that the current stochastic envelope is final.

Pulse +15 comparison completed:
- Baseline periodic pulse already fires at model timesteps `1, 16, 31` because `COTS_PULSE_OFFSET=1` and `COTS_PULSE_PERIOD=15`.
- `src/scenario.jl` now supports optional extra one-off pulses through `COTS_EXTRA_PULSE_TIMESTEPS`.
- `simulate_best_stochastic.jl` and `plot_best_stochastic.py` support tagged outputs via `COTS_OUTPUT_TAG`.
- Generated a `plus15` diagnostic with `COTS_OUTPUT_TAG=plus15` and `COTS_EXTRA_PULSE_TIMESTEPS=16`, which adds an extra pulse at timestep `16` (year 2000) on top of the baseline schedule.
- Outputs: `sandbox/data/best_stochastic_plus15_trajectories.csv`, `sandbox/data/best_stochastic_plus15_site_trajectories.csv`, `sandbox/data/best_stochastic_plus15_metadata.csv`, `sandbox/data/best_stochastic_plus15_validation_metrics.csv`, and `sandbox/best_stochastic_plus15_plot.png`.
- Result: the extra timestep-16 pulse did **not** materially change validation metrics or median peak years. Second-peak years on the median trajectories remained Lizard `2019`, MacGillivray `2022`, North Direction `2017`, and Eyrie `2016`.
- Interpretation: because the baseline already includes a periodic timestep-16 pulse, adding an additional pulse at the same timestep is not enough to bring the second peak forward. A more direct next test would shift or add the later pulse around timestep `29` or `30` (2013-2014) rather than reinforcing timestep `16`.

No-pulse comparison completed:
- Goal was to test whether the external Cairns pulse schedule materially affects the current stochastic best-fit trajectories.
- `simulate_best_stochastic.jl` now respects `COTS_EXTERNAL_PULSE=false` and records `external_pulse_enabled` in metadata.
- Generated a `nopulse` diagnostic with `COTS_OUTPUT_TAG=nopulse`, `COTS_EXTERNAL_PULSE=false`, and no extra pulse timesteps.
- Outputs: `sandbox/data/best_stochastic_nopulse_trajectories.csv`, `sandbox/data/best_stochastic_nopulse_site_trajectories.csv`, `sandbox/data/best_stochastic_nopulse_metadata.csv`, `sandbox/data/best_stochastic_nopulse_validation_metrics.csv`, and `sandbox/best_stochastic_nopulse_plot.png`.
- Result: baseline, `plus15`, and `nopulse` had the same overall median peak years on the plotted reefs: Lizard `1993`, MacGillivray `1993`, North Direction `1994`, Eyrie `1993`.
- Result: baseline, `plus15`, and `nopulse` also had the same second-peak years (`>=2010`) on the plotted median trajectories: Lizard `2019`, MacGillivray `2022`, North Direction `2017`, Eyrie `2016`.
- Validation metrics were effectively unchanged across all three variants. Direct trajectory comparison showed a small local pulse effect (`max adult COTS abs diff ~0.325`, mean abs diff ~0.00032), but not enough to affect reef-level median timing or fit metrics.
- Interpretation: in the current stochastic best-fit workflow, the external pulses are not the dominant driver of the plotted reef-level timing. The second peak appears to be governed more by internal COTS/coral dynamics, demographic perturbations, environmental fields, and/or normalization than by the pulse schedule.

Pulse mechanism redesign completed:
- External COTS pulses are now **off by default** (`COTS_EXTERNAL_PULSE=false`).
- New ENV controls are `COTS_EXTERNAL_PULSE`, `COTS_PULSE_START`, `COTS_PULSE_DURATION`, `COTS_PULSE_REPEAT_INTERVAL`, and `COTS_PULSE_RELATIVE_MAGNITUDE`.
- Relative magnitude scales pulse input to each location as a proportion of that location's running maximum internally generated age-0 larval supply, giving a tunable proxy for external upstream supply until a real external source layer is implemented.
- `simulate_best_stochastic.jl` records the pulse controls in metadata (`external_pulse_enabled`, `pulse_start`, `pulse_duration`, `pulse_repeat_interval`, `pulse_relative_magnitude`).
- `docs/src/concepts/cots_submodel.md` now documents the new pulse controls and replaces the old absolute `pulse_val` description.
- Smoke test with `COTS_EXTERNAL_PULSE=false`, `COTS_N_STOCHASTIC_SCENS=2`, `COTS_OUTPUT_TAG=pulse_smoke` completed successfully; metadata showed pulses disabled and relative magnitude `0.0`.
- Smoke test with `COTS_EXTERNAL_PULSE=true`, `COTS_PULSE_START=16`, `COTS_PULSE_DURATION=1`, `COTS_PULSE_REPEAT_INTERVAL=15`, `COTS_PULSE_RELATIVE_MAGNITUDE=0.5`, `COTS_N_STOCHASTIC_SCENS=2`, `COTS_OUTPUT_TAG=pulse_smoke_on` completed successfully; metadata recorded those settings and trajectories differed from the no-pulse smoke test (`max adult COTS abs diff ~9.48`, mean abs diff ~0.688).
- Example diagnostic command:
  `$env:COTS_OUTPUT_TAG='pulse_timing_test'; $env:COTS_EXTERNAL_PULSE='true'; $env:COTS_PULSE_START='29'; $env:COTS_PULSE_DURATION='1'; $env:COTS_PULSE_REPEAT_INTERVAL='15'; $env:COTS_PULSE_RELATIVE_MAGNITUDE='0.5'; julia sandbox\calibration\simulate_best_stochastic.jl`

Pulse scenario plotting and pulse-only calibration sweep completed:
- Requested scenario plots used `COTS_PULSE_START=20` (~2004/2005 depending timestep convention), `COTS_PULSE_RELATIVE_MAGNITUDE=0.5`, and durations `1`, `2`, and `3` years.
- Added `sandbox/calibration/run_pulse_start20_scenarios.ps1` to generate the three requested start-20 duration variants and their individual plots.
- Generated tagged outputs for `best_stochastic_pulse_start20_dur1`, `best_stochastic_pulse_start20_dur2`, and `best_stochastic_pulse_start20_dur3`.
- Added and ran `sandbox/plotting/plot_pulse_scenario_comparison.py`; outputs are `sandbox/pulse_start20_scenario_comparison.png` and `sandbox/data/pulse_start20_scenario_comparison_metrics.csv`.
- Result: start-20 pulses with duration 1/2/3 and relative magnitude 0.5 barely changed matched-year metrics or peak years compared with `nopulse`; duration 1 and 3 had a tiny mean absolute percent-bias improvement, but no meaningful timing shift.
- Added `sandbox/calibration/pulse_calibration_sweep.jl` to keep current biological parameters fixed and vary only pulse controls (`start`, `duration`, `repeat_interval`, `relative_magnitude`) to rank pulse timing/magnitude fits using matched-year validation metrics.
- Smoke-tested the sweep with `PULSE_SWEEP_STARTS=20`, `PULSE_SWEEP_DURATIONS=1`, `PULSE_SWEEP_REPEAT_INTERVALS=0`, and `PULSE_SWEEP_RELATIVE_MAGNITUDES=0.0,0.5`; it wrote `sandbox/data/pulse_calibration_sweep_by_reef.csv` and `sandbox/data/pulse_calibration_sweep_summary.csv` successfully.

COTS documentation update completed:
- `docs/src/concepts/cots_submodel.md` now documents how the COTS submodel hooks into ADRIA's main `scenario.jl` ecosystem loop.
- Added details on how COTS reads and mutates CoralBlox-managed coral cover through `C_cover_t`, including absolute-vs-relative cover conventions.
- Added state ownership and side-effect boundaries for `cots_timestep!`, `cots_mortality!`, larval dispersal, external supply, and logging.
- Added a forward-looking design section for a possible `CotsBlox.jl` package, including proposed API boundaries, data contracts with CoralBlox, responsibilities ADRIA should retain, migration steps, and risks.
### A. `switch_RCPs!` on `LizardDomain`
When running `ADRIA.run_scenarios(dom, scens, "45")` multi-threaded, internal ADRIA logic calls `switch_RCPs!(dom, rcp)`. We have explicitly added the method to `src/ExtInterface/Lizard/LizardDomain.jl`:
```julia
function ADRIA.switch_RCPs!(d::LizardDomain, RCP::String)::LizardDomain
    return d
end
```
*Never remove this method from `LizardDomain.jl` unless `LizardDomain` is refactored to inherit dynamic multi-RCP paths.*

### B. Sampling Counterfactual Stochastic Scenarios (`guided = 0`)
If you call `scens = ADRIA.sample(dom, N)` to generate stochastic runs across environmental/demographic variables, `sample()` will by default perturb **Intervention / Seeding** parameters (`N_seed_TA`, `N_seed_CA`, `guided`, `fogging`, etc.).
If left active, those seeded corals will flood the domain and cause:
`Warning: Number of seeding devices exceeds available space. Excluding devices...`
* **Fix:** Always zero out intervention columns when testing counterfactual COTS dynamics:
  ```julia
  scens.guided .= 0
  for col in names(scens)
      if startswith(col, "N_seed_") || col in ["fogging", "SRM", "a_adapt"]
          scens[!, col] .= 0.0
      end
  end
  ```

### C. Serial vs Parallel Execution & Memory
* **`run_scenario(dom, p)`** runs a single scenario sequentially in ~10–15 seconds. For small ensembles (`N <= 20`), a simple `for` loop in `simulate_best_stochastic.jl` (`20 * 12s ≈ 4 minutes`) is completely robust and avoids Zarr batch I/O overhead.
* **`run_scenarios(dom, scens, "45")`** runs multi-threaded and uses Zarr chunked execution. Use this when sampling `N >= 64` scenarios.

---

## 6. Next Phases & Action Roadmap

### Phase 9a/c: COTS-Specific Larval Dispersal Matrix
Currently, COTS larval dispersal (`disperse_cots_larvae!`) uses coral connectivity as a proxy.
* **Action:** Source or generate an empirical hydrodynamic COTS connectivity matrix for the Lizard Island cluster (`Lizard_COTS_Connectivity.csv`). Update `build_lizard_domain.jl` and `scenario.jl` to pass this COTS-specific connectivity into `disperse_cots_larvae!` separately from coral connectivity.

### Phase 10: Formal Global Sensitivity Analysis (Sobol / PAWN)
* **Action:** Execute a formal Sobol or PAWN sensitivity sweep using `ADRIA.sensitivity` across the `[a_F, a_S, IMM, seed_mult, pulse_val]` parameter space around the `run_id = 2` regime to quantify the exact first-order and total-order sensitivity indices driving cycle period and trough depth.

### Phase 11+: Regional & Full-Scale GBR Cross-Validation
* **Action:** Validate the calibrated COTS parameter set against empirical AIMS tow datasets across broader Moore/Cairns regional domains (`RMEDomain`) to confirm generalizability outside the Lizard Island cluster.

Post-rebase smoke check completed (2026-07-10):
- Restored active-package COTS includes in `ADRIA/src/ADRIA.jl` (`ecosystem/cots_factors.jl`, `ecosystem/cots.jl`). Without these, `CotsParams` and COTS runtime types were missing after the rebase.
- Migrated the existing Lizard custom loader from old root `src/ExtInterface/Lizard/LizardDomain.jl` into active package path `ADRIA/src/ExtInterface/Lizard/LizardDomain.jl`, included/exported `LizardDomain`, and kept `switch_RCPs!` for the historical Lizard workflow.
- Fixed `ADRIA/src/io/ResultSet.jl` rebase duplication that left an invalid/incomplete `cots_condition_log` ternary and duplicated COTS log loading block.
- Adjusted `ADRIA/src/scenario.jl` to use the installed CoralBlox API: removed `LinearExtensionCache` and pass `_bin_edges` directly to `linear_extension_scale_factors`.
- Fixed `sandbox/Manifest.toml` conflict markers and pointed its ADRIA path at `../ADRIA` instead of stale root `..`; `sandbox/Project.toml` now declares direct script deps `CSV` and `DataFrames`.
- Updated `sandbox/calibration/simulate_best_stochastic.jl` and `sandbox/calibration/pulse_calibration_sweep.jl` to activate `sandbox/Project.toml` based on `@__DIR__` and `cd` to repo root, so they are not sensitive to launch directory.
- Verification passed: `julia --project=sandbox -e "using ADRIA"` precompiled/loaded successfully.
- Verification passed: two-scenario no-pulse stochastic smoke completed with `COTS_N_STOCHASTIC_SCENS=2`, `COTS_OUTPUT_TAG=post_rebase_smoke`, `COTS_EXTERNAL_PULSE=false`; outputs written under `sandbox/data/best_stochastic_post_rebase_smoke_*`.
- Verification passed: plot generation completed with `COTS_OUTPUT_TAG=post_rebase_smoke`; output `sandbox/best_stochastic_post_rebase_smoke_plot.png` and validation metrics CSV were written.
- Note: `sandbox/Manifest.toml` is currently staged from conflict resolution in Git; other repair files are unstaged. Review staging before commit.