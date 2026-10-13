# COTS calibration: current state and handoff (updated 2026-10-13, after runs 11–18)

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

### Leading COTS candidate (COTS dynamics), re-baselined 2026-10-13

**Current baseline:** run `20261013T_11_cots_rebaseline_1991_v1` arm `wc_b150_q20` (replayed exactly as `wc_b150_q20_acrpoc` in run 12).
- **Change from run 15.** The same COTS model now runs on the confirmed working coral configuration (below), with predator pit **`Q = 20`** (re-checked; interior of 10–70).
- **Score.** Summed frozen loss **9.06** (run 15: 12.18). All three two-wave reefs have two peaks, but 6/7 match. **Gate fails** on matched peaks and the Lizard trough ratio (0.25 against ≤ 0.10).
- **β.** β = 1.5 was chosen over 1.0 by coral `J` (0.364 vs 0.369, a negligible margin). β 1.0 at `Q = 20` scores 9.12.

**Per reef** (COTS/tow at `q = 0.015`; observed in brackets):

| Reef | First peak | Second wave |
| --- | --- | --- |
| Lizard | 0.35 (1.92) | 0.16 (0.74) |
| MacGillivray | 0.35 (0.86) | 0.40 (0.67) |
| North Direction | 0.62 (0.56) | 0.83 (0.71) |
| Eyrie | 0.75 (about 0.13) | 0.75 (0.85) |

Troughs are 2.7 adults/ha.

**Detection-consistent candidate (labelled; not adopted):** run `20261013T_17_initial_state_rebaseline_1991_v1` arm `wc_b150_q20_eb_h10`, read at `A½ = 10` (run 18 re-score).
- **How it differs.** It is the same model with empirical-Bayes 1991 seeds, inverted and observed through `CPUE = q·A²/(A + 10)`.
- **Score.** It passes the qualitative gate under the raw peak count (7/7; worst trough 0.088), with a summed loss of 9.92 on its own observation model.
- **Caveats.** It depends on an unidentified `A½` and the lead's gate-count decision. First-wave heights remain far too low.

The previous candidate is described next for provenance.

**Previous candidate:** run `20261012T_15_food_grace_1991_v1` arm `acrpoc_g0`, identical to run 14's `obs_m100_q30`.

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

**Working coral configuration (confirmed by the project lead, 2026-10-13):**

| Component | Setting |
| --- | --- |
| Growth | CoralBlox untouched (`s = f = 1`) |
| Background | Mellin table + 3%/yr chronic loss |
| Cyclone strength | κ = 5 |
| Bleaching scalar | β = 1–1.5 |

In the forced runs this gives `J` about 0.32–0.34 and a 2017–19 trough of about 0.13–0.16 (observed 0.07).

**Residual.** Fast coral is about 2× too high in 2005–13 under observed COTS.

## Where the COTS work got to: open problems

Work order set by the lead (2026-10-13): (1) re-baseline, (2) first-wave height, (3) second-wave food, (4) the 1999 crash, (5) manta detectability. The lead then brought item 5 forward and chose the reef-differentiated starting state for item 2. Items 1, 5 (as a sensitivity) and 2a are done.

1. **Re-baseline: done (run 11).** Loss 12.18 → 9.06 at `Q = 20`, β 1.5. The Mellin + 3% background carries the gain; κ 5 alone is worse.
   - Slow overgrowth is reduced, but still 2–3× the photo transects.
   - Inter-wave fast coral on Lizard and MacGillivray did **not** rise (about 0.027), because adults stay at 2–12/ha and subadults at 11–130/ha between waves.
2. **First-wave height: composition closed (run 12).** The plausible composition range reaches MacGillivray's first peak, but not Lizard's (≤ 0.62 vs 1.92), and it raises troughs above 4 adults/ha.
   - **The between-reef ranking of first waves is inverted:** observed Lizard > MacGillivray > North Direction > Eyrie; model Eyrie > North Direction > Lizard ≈ MacGillivray.
   - **Causes:** the 1989–91 IDW initial state (North Direction 3.3 vs Lizard 0.5 adults/ha), and near-uniform Ext×5 boundary supply (70–95% of first-wave recruits; local retention ≤ 1.4/ha/yr).
   - **Starting state: closed (run 17).**
     - **Change tested.** An empirical-Bayes surface (`build_1991_eb_surface.py`) shrinks single-animal surveys; North Direction's seed falls from 3.3 to 0.8 adults/ha. Seeds were also inverted with the detection curve.
     - **Result.** First waves barely change, because they are set by supply and local production in 1992–96. The ranking stays inverted (rank correlation −0.8).
     - **Remaining levers (lead to choose):** local retention or self-recruitment at initiation-region reefs, or spatially structured external supply.
3. **Second waves.** North Direction and MacGillivray now build them. Lizard's stays small (0.16 vs 0.74).
   - The same coral model gives fast cover **about 2× too high** under manta-forced COTS and **too low** under COTS dynamics.
   - The difference is trough COTS density, so this is tied to item 5.
4. **Abrupt 1999 crash.** Not reproduced. The juvenile pipeline (80–130 subadults/ha in 1999–2002) refills adults.
   - **Observed.** The drop is synchronous: Lizard 1.92 → 0.05, MacGillivray 0.40 → 0.00 and North Direction 0.17 → 0.00, all between 1998 and 1999.
   - **Partial lead (run 17).** The detection-consistent pairing (`eb_h10` read at `A½ = 10`) makes the post-peak drop abrupt on MacGillivray (one-year ratio 0.09) and North Direction (0.19), but not Lizard (0.46). Modelled peaks fall in 1997, a year before the observed crash year.
   - **Tested and failed or not recommended:** ReefMod-style die-off (`COTS_OUTBREAK_COLLAPSE`), post-grazing starvation, shorter memory, and one-year starvation grace (evidence supports about 6 months).
5. **Manta detectability: tested as a labelled observation sensitivity (runs 13/14, 18).**
   - **Curve.** The lead chose `CPUE = q·A·A/(A + A½)`. The [explanatory plot](runs/20261013T_14_detectability_rescore_1991_v2/manta_detection_curve.png) and the re-scoring tool (`rescore_detectability.jl`) are checked in.
   - **At `A½ = 10`.** Troughs read like the near-zero manta counts (worst trough ratio 0.25 → 0.10), but peaks read about 22% lower, so the first-wave deficit grows. The frozen loss rises.
   - **`A½` is not identified.** It needs independent evidence, such as paired manta–SCUBA or culling data.
   - **The frozen observation model is unchanged.**
6. **Promotion gates.**
   - **Audit finding (run 14).** The gate's "7 observed peaks" exist only on the **raw** series. Audits counted smoothed matches, where at most 6 are possible, so "7/7" was unattainable as audited.
   - **Under the raw count.** The baseline already matches 7/7 and fails only the trough ratio. **The lead must decide the gate's peak count.**
   - **Arms that pass the raw-count gate:**
     - the detection-consistent `eb_h10` at `Q = 20` (worst trough 0.088), read at `A½ = 10`;
     - several re-scored arms at `Q ≥ 30`.
   - **Not promotable.** These passes rest on an unidentified `A½`, and first-wave heights remain far too low.
   - **After a pass.** Replicated seeds, `q` sensitivity and held-out testing come only after an accepted pass.

**Suggested next steps (lead to decide):**
- (i) **Two decisions:**
  - the gate's peak count: raw 7, or smoothed 6;
  - whether to seek evidence for `A½`, so that the detection-consistent configuration (`eb` seeds inverted at `A½`, read at `A½`) can become the working model.
- (ii) First-wave height and ranking: a supply mechanism (local retention or self-recruitment in the Lizard–Cairns initiation region, or spatially structured external supply), preregistered with evidence.
- (iii) Items 3–4 (Lizard's second wave; Lizard's abrupt crash) on the configuration (i) settles.

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
| 10 | COTS re-baseline: failed before any arm (runner lowering error; `FAILED.txt`) |
| 11 | COTS re-baseline on the working coral (new baseline `wc_b150_q20`, loss 9.06) |
| 12 | first-wave composition range (no level adopted; ranking of first waves inverted) |
| 13 | detectability re-score v1 (smoothed only; superseded) |
| 14 | detectability re-score v2 (raw + smoothed; gate-count audit finding) |
| 15 | empirical-Bayes 1991 COTS surface (derived input) |
| 16 | starting-state screen: failed before any arm (exact float check; `FAILED.txt`) |
| 17 | reef-differentiated starting state (first waves insensitive; `eb_h10` passes the raw-count gate at `A½ = 10`) |
| 18 | detectability re-score of run 17 |

## Known issues

- **Runtime slowdown: fixed 2026-10-13.** `set_initial_cover!` rebuilt DataCubes from `collect(cube.species)`, which nested DimensionalData lookups each arm. The fix captures plain label vectors once.
- **DimensionalData wrappers (general rule).** Always `collect`/`Array` lookups into plain vectors before reusing them. Three bugs have come from this.
- **ADRIA cannot run `k = 0` sites.** Runs floor `k` at 0.01. The upstream fix is not made.
- **ADRIA misindexes non-contiguous `CB_CALIB_GROUPS`.** The runner relabels them. The Lizard loader defaults every site to group 1; the RME layer gives groups 5/1/4.
- **Missing metadata.** Joint-fit parts 02 and 04 lack `metadata.toml` (stopped processes). This is disclosed in the protocol.
- **Parse test blind spot.** `test_calibration_scripts_parse.jl` only calls `Meta.parseall`, so lowering errors (for example, mixing `,` and `;` in an array literal) pass and then fail at run time (run 10). A `Meta.lower` check on each top-level expression catches them.
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
