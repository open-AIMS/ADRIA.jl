# COTS calibration: current state and handoff (updated 2026-10-13)

This is the single entry point for resuming work, in this chat or a new one. It records:
- where the work stands;
- how to run things;
- what is decided and what is open;
- what to do next.

Nothing here is promoted to the ADRIA default model.

**Read in this order:**
1. this page;
2. repository-root `AGENTS.md` (the rules);
3. the end of [crash_mechanism_protocol.md](crash_mechanism_protocol.md) (latest preregistrations and decisions);
4. the latest rows of [MODEL_LOG.md](MODEL_LOG.md);
5. [README.md](README.md) (reproduce commands per screen).

## Working agreements with the project lead

- **Workflow.** Every screen follows the same steps:
  1. diagnose with data;
  2. compare with ReefMod, CoCoNet and MICE;
  3. **preregister** in the protocol before implementation;
  4. add an opt-in switch whose default reproduces current behaviour, plus tests;
  5. run a small factorial with an **exact-replay reference arm**;
  6. audit;
  7. log in MODEL_LOG, the protocol, README and this page.
- **Decisions belong to the project lead.** Ask before changing the frozen observation model, retuning failed parameters or adopting inputs.
- **The lead commits.** Do not commit unless asked. Report uncommitted state in both repos.
- **Coral growth principle (2026-10-13).** Keep CoralBlox calibrated growth rates untouched as far as possible; they were calibrated to LTMP data.
- **Data restrictions.** `SpatialData/GBR-HabitatMapping_INTERNALONLY` is internal only, and `GBRMPA_RRN_2025_AIMS.xlsx` is under a data licence. Derived outputs stay in the git-ignored `sandbox/data/`. Never publish or commit them.
- **Run naming.** Runs are named `runs/YYYYMMDDT_NN_<experiment>_1991_vN`, sequential per date. The runner auto-numbers IDs given as `YYYYMMDDT_<name>`. Older runs were renamed by `renumber_runs.py`; old → new names are in `runs/RUN_RENAMES.csv`.

## Where the code lives

| Piece | Location |
| --- | --- |
| COTS equations, switches, pure tests | **Separate repo** `..\COTSMod.jl` (`src/COTSMod.jl`, `test/runtests.jl`; 363 tests). Loaded by local path from `sandbox/Manifest.toml`, so uncommitted edits affect runs. Committed by the lead through `0a45f58`. |
| ADRIA adapter | `ADRIA/src/ecosystem/cots.jl`: `_cots_runtime_params`, `load_cots_forced_adults`, diagnostics. |
| Env switches and logs | `ADRIA/src/scenario.jl`: `COTS_*` env vars. `cots_mechanism_log` has 10 channels: food memory, food survival, senescent adults, senescent deaths, diet share, juvenile predation, maturation fraction, outbreak clock, collapse removals, low-food years. |
| Experiment runner | `sandbox/calibration/run_1991_initialization_test.jl`, one `OWEN_1991_EXPERIMENT` per screen. |
| Audits | `sandbox/calibration/audit_1991_<experiment>.py`. |
| Derived inputs | See the table below. |
| Run index | `python sandbox/calibration/build_model_run_inventory.py` |
| Adapter tests | `julia --project=sandbox sandbox/calibration/test_1991_initial_state.jl`; parse check `test_calibration_scripts_parse.jl`. |

**Derived inputs** (all checked-in builders):

| Builder | Output |
| --- | --- |
| `sandbox/domain_building/build_site_habitable_fraction.py` | site habitable fraction |
| `build_cyclone_hours_rrn.py` | cyclone hours |
| `build_lizard_domain_v3.py` | the v0.3 domain |
| `sandbox/calibration/build_observed_cots_forcing.py` | observed COTS forcing |
| `export_mellin_disturbance.R` + `build_mellin_background_mortality.jl` | Mellin background |

**Tools:**
- **Julia:** `julia --project=sandbox`.
- **Python:** base miniconda (`C:\Users\smatthew\AppData\Local\miniconda3\python.exe`: pandas, geopandas, scipy, openpyxl). The user Python has rasterio (GDAL 3.9) for the FileGDB rasters.
- **R:** `C:\Users\smatthew\AppData\Local\Programs\R\R-4.5.0\bin\x64\Rscript.exe`, with base packages only.
- **Shell:** Bash heredocs with apostrophes, and Python f-strings containing `{...}`, have both broken patch scripts. Write patch scripts with the Write tool.

**Running a screen.** Set the environment, then run the runner:
```powershell
$env:OWEN_1991_EXPERIMENT = '<experiment>'
$env:OWEN_1991_RUN_ID = '<YYYYMMDD>T_<experiment>_1991_v1'
julia --project=sandbox sandbox/calibration/run_1991_initialization_test.jl
```
- **Arm filter.** `OWEN_1991_ARM_FILTER` (regex on arm names) runs a subset, for example to complete an interrupted grid in a new run directory.
- **Timing output.** Each arm prints its time split into `run_scenario`, setup checkpoints and post-processing.
- **Duration.** About 1.5 minutes for domain load and compile, then about 10–20 s per arm.

## Model state

### Leading COTS candidate (COTS dynamics)

Run `20261012T_15_food_grace_1991_v1` arm `acrpoc_g0`, identical to run 14's `obs_m100_q30`.

**Working COTS model** ("condition, no switching", the lead's preference):

| Component | Setting |
| --- | --- |
| Grazing | per-capita 10 m²/adult/yr, preferred-first, ÷ site `k` |
| Starvation (H1) | preferred-prey starvation, θ 0.05 on remembered fast cover, memory 0.5 yr |
| Senescence (H2) | 6+ plus group at 0.8/yr |
| Condition | preferred-food condition, saturation 0.10, memory 1 yr |

**Labelled assumptions:**
- **External ×5** on the evidence-derived outbreak boundary production;
- **predator pit `Q = 30`** COTS/ha/yr (tuned);
- **gradual Allee** `exp(-(3/A)^0.5)`.

**Conventions:** `q = 0.015` COTS/tow per adult/ha; `k_unit = 19.61`; Allee scale 3; `ADRIA_RNG_SEED=90972102`; Owen annual connectivity seed 20260930.

**Inputs:**
- `Lizard_Historical_v0.3`: geomorphic-infilled site `k`; RRN cyclones at κ = 1;
- CoralBlox **group 5** via per-reef RME groups, relabelled contiguously;
- fast-growth multiplier 1.0;
- **observed 1993–95 photo-transect composition** (*Acropora* + Pocilloporidae share: Lizard 0.21, MacGillivray 0.13, North Direction 0.25, others 0.19);
- inherited 2023 total cover.

**Performance:**
- **Score.** Summed frozen loss 12.18 (best to date); 5/7 peaks matched.
- **First peaks.** 1997–99, against observed 1996–98, but too low on Lizard and MacGillivray (about 0.45 vs 1.92 / 0.86 COTS/tow).
- **Troughs.** 2–3 adults/ha, inside the observed 0–3.
- **Second waves.** Small: Lizard 0.13–0.17 vs 0.74; North Direction 0.45 vs 0.71; Eyrie 0.64 vs 0.85.

### Coral side (COTS forced from observed CPUE / q)

These runs isolate coral dynamics from the COTS model.

- **Under observed COTS, the group-5 CoralBlox model overgrows both coral types on these reefs:**
  - slow coral to about 0.5 vs observed 0.16–0.24;
  - fast coral to 0.14–0.24 vs observed 0.04–0.11 in 2005–13;
  - recovery after 2015–17 is too fast.
- **Why.** CoralBlox was calibrated over 2008–22 with *observed* LTMP storm/COTS losses imposed, whereas these runs use weaker modelled disturbances.
- **Calibration-group provenance check.** Group membership is identical to the canonical file used in the CoralBlox calibration. The ADRIA factors come from commit `f931cfa42`.
- **Joint fit** (growth × cyclone κ × background μ): `J` 0.54 → 0.24, but only by cutting growth (s 0.55, f 0.8) with μ = 6%/yr, both at bounds.
- **Mellin background** (lead's COTS_Model: disease/unknown log-cover effects, 0.39%/yr plus recorded LTMP events) helps in every pair but is too small.
- **Bleaching sensitivity with growth untouched:**
  - `J` reaches 0.284 at Mellin + 3%/yr, κ 2.5, bleaching ×2.5;
  - but ×2.5 kills fast coral in 1994 and stops post-2020 recovery;
  - a ×1–1.5 scalar is more plausible;
  - bleaching cannot fix the 2005–13 fast excess, because these reefs have no heat stress over 1995–2015.

**Recommended working coral configuration (project lead to confirm):**

| Component | Setting |
| --- | --- |
| Growth | CoralBlox untouched (`s = f = 1`) |
| Background | Mellin table + 3%/yr chronic loss |
| Cyclone strength | κ = 5 |
| Bleaching scalar | β = 1–1.5 |

In the forced runs this gives `J` about 0.32–0.34 and a 2017–19 trough of about 0.13–0.16 (observed 0.07).

**Residual.** Fast coral is about 2× too high in 2005–13 under observed COTS.

## Where the COTS work got to: open problems

1. **Re-baseline needed.** The leading COTS candidate has **not** been run with the recommended coral configuration (κ 5, Mellin + 3%, β 1–1.5, growth untouched). Coral supply, and so outbreak size and second waves, will change. Re-check `Q`.
2. **First-wave height** is too low on Lizard and MacGillivray with the observed composition. The composition is uncertain:
   - photo transects are at 9 m and miss crest *Acropora*;
   - "Other" hides *Montipora* / *Isopora*;
   - the plausible range runs up to Acropora + Pocilloporidae + Other, mean 0.36.
3. **Second waves are small.** Fast coral on Lizard and MacGillivray stays at or below the starvation threshold θ 0.05 in the dynamic runs. Under observed COTS the coral model instead *overshoots* fast cover, so the dynamic model's trough grazing (2–3 adults/ha) and space competition interact.
4. **Abrupt 1999 crash.** Not reproduced. The juvenile pipeline (about 150 subadults/ha at the peak) refills adults. Tested and failed or not recommended:
   - ReefMod-style die-off (`COTS_OUTBREAK_COLLAPSE`);
   - post-grazing starvation;
   - shorter memory;
   - one-year starvation grace (evidence supports about 6 months).
5. **Manta detectability (untested).** Size- or density-dependent detectability, which could reconcile model trough adults with near-zero manta counts. Testing it changes the observation model and needs the lead's approval.
6. **Promotion gates.** No arm has passed the qualitative gate: 2 peaks on Lizard, MacGillivray and North Direction; 7/7 matched; worst trough ratio ≤ 0.10. Replicated seeds, `q` sensitivity and held-out testing come only after a pass.

**Suggested next steps:**
- (i) confirm the working coral configuration, then re-baseline the leading COTS candidate on it (with a reference replay);
- (ii) a composition-range sensitivity for first-wave height;
- (iii) revisit second-wave food supply once (i) is known.

## Experiment chronology (runs in `runs/`)

2026-10-12:

| Run | Result |
| --- | --- |
| 01–02 | switching retired |
| 03–04 | habitat |
| 05–06 | benthic `k` + cyclones |
| 07 | trough mechanisms: predator pit, maturation delay, gradual Allee |
| 08–10 | initiation: Ext×5 adopted; v1/v2 failed before simulation |
| 11 | `Q` tuning (`Q = 30`) |
| 12 | ReefMod die-off (counterproductive) × coral group |
| 13 | growth × starvation (memory 0 and post-grazing not recommended) |
| 14 | observed composition (troughs fixed, first peaks low) |
| 15 | composition range × starvation grace (grace not recommended; coral space limits fast coral) |
| 16 | coral calibration v1 (failed: DimVector loader bug) |

2026-10-13:

| Run | Result |
| --- | --- |
| 01 | coral calibration v2 (growth only; edges) |
| 02, 04, 05 | joint coral fit, split by the slowdown |
| 03 | Mellin background |
| 06–08 | runtime diagnostics |
| 09 | bleaching sensitivity |

## Known issues

- **Runtime slowdown: fixed 2026-10-13.** `set_initial_cover!` rebuilt DataCubes from `collect(cube.species)`, which nested DimensionalData lookups each arm. The fix captures plain label vectors once.
- **DimensionalData wrappers (general rule).** Always `collect`/`Array` lookups into plain vectors before reusing them. Three bugs have come from this.
- **ADRIA cannot run `k = 0` sites.** Runs floor `k` at 0.01. The upstream fix is not made.
- **ADRIA misindexes non-contiguous `CB_CALIB_GROUPS`.** The runner relabels them. The Lizard loader defaults every site to group 1; the RME layer gives groups 5/1/4.
- **Missing metadata.** Joint-fit parts 02 and 04 lack `metadata.toml` (stopped processes). This is disclosed in the protocol.
- **Tracked caches.** `sandbox/calibration/__pycache__/*.pyc` files are tracked or showing in git. Consider adding them to `.gitignore`.

## How to add the next screen

1. **Preregister.** Append to `crash_mechanism_protocol.md`: hypothesis, equations, evidence, arms including an exact-replay reference, gate, failures, rollback.
2. **Switch and tests.** Add an opt-in COTSMod switch whose default reproduces current behaviour, plus tests: `julia --project=. -e "using Pkg; Pkg.test()"` in `..\COTSMod.jl`.
3. **ADRIA wiring.** Add the switch to `_cots_runtime_params` and to the env parsing and validation in `scenario.jl`, plus an adapter test.
4. **Runner.** Add the experiment to the runner's `EXPERIMENT` list, `MIGRATION` list, treatment list and the per-arm flag blocks. Run the parse test.
5. **Run.** Use an unnumbered `YYYYMMDDT_<name>` run ID; it is numbered automatically.
6. **Audit.** Copy the nearest `audit_1991_*.py`. Check the reference replay, then inspect the PNG.
7. **Log.** Add the MODEL_LOG row (one line per table row: no line breaks inside cells), the protocol decision and the README section. Rerun the inventory, then update this page.

## Starting a new chat

Paste something like:

> Read `sandbox/calibration/CURRENT_STATE.md`, then `AGENTS.md`, and the end of `sandbox/calibration/crash_mechanism_protocol.md`. We are calibrating the COTS model in ADRIA (with sibling repo `..\COTSMod.jl`) on the Lizard v0.3 domain. Confirm the working coral configuration with me, then re-baseline the leading COTS candidate on it, following the preregister → opt-in → exact-replay → audit → log workflow.
