# COTS calibration: current state and handoff (updated 2026-10-14)

This is the single entry point for resuming work, in this chat or a new one. It records:
- where the work stands;
- how to run things;
- what is decided and what is open;
- what to do next.

Nothing here is promoted to the ADRIA default model.

**Resume here:** the section **"Close-out plan (2026-10-14)"** below. Done: A1, A2, B1, B2, B3 (`candidates/candidates.toml`), B4 replicates, C1 and C2 (switch plus screen). Next: the lead's decisions, listed under "Open for the lead" in Progress.

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
| `build_observed_cots_forcing.py --year-basis survey` | survey-year COTS forcing (`observed_cots_forcing_1991_2024_survey_year/`) |
| `build_1991_eb_surface.py` | empirical-Bayes 1991 COTS surface (`runs/20261013T_15_1991_initial_surfaces_eb`) |
| `build_initial_share_range.py` | healthy-year 1991 fast-share levels (`derived/initial_fast_share_range/`) |
| `build_survey_timing.py` | LTMP survey timing vs bleaching (`derived/ltmp_survey_timing/`) |
| `analyse_manta_detectability.py` | paired diver–manta detectability (`derived/manta_detectability_pairs_lizard/`) |

**Scoring and evaluation tools (2026-10-13/14):**
- `rescore_detectability.jl [--q=<q>]`: re-scores stored trajectories under `CPUE = q·A·A/(A + A½)` for `A½` in {0, 1, 3, 10, 30}, both observation treatments. `A½ = 0` at `q` 0.015 replays exactly.
- `evaluate_gate_g2prime.py <out> <rescore dirs> [--trough-reefs …]`: gate G2′ (base miniconda Python).
- `connectivity_retention.jl` with `test_connectivity_retention.jl`: source-conserving larval retention.
- `test_bleaching_floor.jl`: tests the ADRIA switches `ADRIA_BLEACHING_DHW_FLOOR`, `ADRIA_BLEACHING_LOW_DHW_SCALE` and `ADRIA_BLEACHING_LOW_DHW_THRESHOLD` (opt-in; unset = default path). These switches are **uncommitted ADRIA core changes** in `growth.jl` and `scenario.jl`.

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

### Confirmed working configuration (project lead, 2026-10-14)

**Run** `20261014T_25_initial_share_1991_v2`, arm `wc_b150_q20_eb_lin_cv2_pick` (runner experiment `initial_share`, with `OWEN_INITIAL_SHARE_PICKS="Lizard Island Reef=min;MacGillivray Reef=q50;North Direction Reef=q25"`).

| Component | Setting |
| --- | --- |
| Coral, "cv2" | CoralBlox growth untouched; Mellin + **12%/yr** chronic background (identified in re-fit v2, run 20); cyclone κ 5; bleaching scalar β 1.5; bleaching mortality **×0.25 below 4 DHW** (`ADRIA_BLEACHING_LOW_DHW_SCALE=0.25`) |
| 1991 fast share (partially calibrated within healthy-year photo-transect ranges) | Lizard 0.147 (min), MacGillivray 0.129 (q50), North Direction 0.307 (q25), other reefs 0.255 (q50 default); totals preserved |
| COTS | working model (per-capita grazing, H1 starvation, H2 senescence, preferred condition); predator pit `Q = 20`; Ext×5; gradual Allee |
| 1991 COTS seeds | empirical-Bayes surface ÷ `q` (linear), 6:3:1 preload |
| Observation model | frozen `q = 0.015` (a Pratchett re-anchor to 0.030 is the next step) |
| Score | summed frozen loss **7.97**; coral `J` **0.260**; 7/7 raw peaks matched |
| Gate G2′ (trough gated on Lizard and MacGillivray) | fails at `A½ = 10` (troughs 0.056 / 0.070 against ≤ 0.054) |

**Known weaknesses:**
- **Late first waves.** In survey-date terms: North Direction 1999 against October 1994; MacGillivray 1997 against October 1995.
- **North Direction's first wave** (1.96 COTS/tow) far exceeds the observed 0.56.
- **Lizard's first wave** (29 adults/ha) is below Pratchett's diver-observed 40–55/ha.

### Earlier candidates (superseded by the confirmed working configuration; kept for provenance)

**Current baseline (2026-10-13, after the lead adopted the empirical-Bayes 1991 surface):** run `20261013T_17_initial_state_rebaseline_1991_v1` arm `wc_b150_q20_eb_lin` (replayed exactly as `wc_b150_q20_eb_lin_owen` in run 20).
- **Model.** The same model as `wc_b150_q20` below, but seeded from the EB surface.
- **Score.** Loss 9.10, 7/7 raw · 6/6 smoothed peaks matched.
- **Gate.** It fails G2′ only on the absolute trough test.

**Previous baseline:** run `20261013T_11_cots_rebaseline_1991_v1` arm `wc_b150_q20` (replayed exactly as `wc_b150_q20_acrpoc` in run 12).
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
   - **Under the raw count.** The baseline already matches 7/7 and fails only the trough ratio.
   - **Decided (2026-10-13).** The lead replaced the gate with **G2′**: both counts, distinct waves, and an absolute trough test (see the protocol). Under G2′ the remaining failure is North Direction's inter-wave level.
   - **Arms that pass the raw-count gate:**
     - the detection-consistent `eb_h10` at `Q = 20` (worst trough 0.088), read at `A½ = 10`;
     - several re-scored arms at `Q ≥ 30`.
   - **Not promotable.** These passes rest on an unidentified `A½`, and first-wave heights remain far too low.
   - **After a pass.** Replicated seeds, `q` sensitivity and held-out testing come only after an accepted pass.

**Lead's decisions (2026-10-13 and 2026-10-14):**
- **Gate.** G2′, with the trough test gated on Lizard and MacGillivray.
- **Starting state.** The empirical-Bayes 1991 surface is adopted.
- **Local retention.** Closed.
- **Time convention.** Keep the LTMP reporting year for scoring, but **track survey timing** (`build_survey_timing.py`; most surveys are made in Sep–Dec, before the scored year's bleaching summer).
- **Observation model.** Try a Pratchett-anchored `q`.
- **Coral.** Re-fit chronic loss under floors 3 and 4.
- **Lizard initialisation.** Test more hard coral (with RHIS for relative levels).

**Lead's decisions (2026-10-14, third round):**
- No crash mechanism.
- Identify the chronic background μ.
- A strong reduction, not a floor, for bleaching below 4 DHW: the 1993–94 transect loss was early-outbreak COTS (corrected).
- Partially calibrate the 1991 fast share within each reef's manta-limited photo-transect range.

**Where it stands (2026-10-14, after runs 20–27).**
- **Recommended working configuration (lead to confirm):**
  - **Coral "cv2":** low-DHW bleaching ×0.25 below 4 DHW (`ADRIA_BLEACHING_LOW_DHW_SCALE`), Mellin + 12%/yr (identified), κ 5, β 1.5, growth untouched.
  - **1991 fast share, partially calibrated:** Lizard 0.15, MacGillivray 0.13, North Direction 0.31, others 0.25.
  - **COTS:** `eb` seeds, `Q` 20, Ext×5, gradual Allee.
  - **Score:** summed loss **7.97** (from 9.10) and coral `J` **0.260** (from 0.363), at run `20261014T_25_initial_share_1991_v2` arm `wc_b150_q20_eb_lin_cv2_pick`.
- **Gate.** G2′ fails on Lizard and MacGillivray troughs at `A½ = 10` (0.056 and 0.070 against ≤ 0.054).
- **Caveats:**
  - North Direction's first wave (1.96 against 0.56) is barely penalised, because smoothing removes its observed 1995 peak.
  - Lizard's loss-picked share (minimum; 29 adults/ha) conflicts with Pratchett's 40–55. The median share gives 53/ha. This is coupled to the `q` choice.
- **Earlier baseline:** `wc_b150_q20_eb_lin` (loss 9.10, coral `J` 0.363), which passes G2′ at `A½ = 10`.

## Close-out plan (2026-10-14)

**Lead's decisions:**
- (1) The working configuration above is confirmed.
- (2) Re-score it under Pratchett's `q` (0.030).
- (3) Score North Direction's first peak against the **raw** series (its 1995 peak counts), and do **not** penalise Eyrie for a missing observed first peak (it was probably missed by sparse surveys).
- (4) **Goal:** close out with a short list of the most plausible candidate parameter sets for simulations, fix the late peak timing, promote the best options, then introduce COTS removal from the culling data.
- **Closed:** no crash mechanism.

**Method rules still apply.** Every step: preregister → opt-in switch with a default replay → exact-replay reference arm → audit → log. Never overwrite runs. The lead commits.

**Progress (2026-10-14).**
- **A1 done.** Score v2 is implemented (`score_v2.jl`, `test_score_v2.jl` 18/18; `loss_v2` and `score_v2_primary` columns in `rescore_detectability.jl`; `summed_loss_v2` in the evaluator).
  - Run 28 re-scores run 25 under v2: cv2 min 6.45, `ref` 6.85, `acrpoc` 7.33, working configuration 8.30.
  - Run 29 is G2′.
- **A2 done.** Score-only `q` 0.030 (runs 30–31), and consistent `q` 0.030 arms in the timing run. `q` 0.030 is worse under v2 and doubles troughs.
- **B1 done.** The rise is about a year late everywhere (external settlement plus maturation). The peak year is set by fast-coral exhaustion, so North Direction's q25 share delays and inflates its wave.
- **B2 done** (run 32; re-scores 33–34; G2′ 35–36; `timing_audit.json`).
  - **C4 = `wc_b150_q20_eb_lin_obs015_cv2_pnd`** (C1 with North Direction share at min, 0.191): v2 **6.44**, coral `J` 0.251, North Direction first peak 1997 at 0.93/tow.
  - **C2 = the same at `q` 0.030** (v2 8.57).
  - Lower `Q` is not selected. No arm passes the G2′ troughs (0.056 / 0.070 against 0.054).
- **C1 done.** `build_cull_removals.py` produces `sandbox/data/derived/cots_cull_removals_lizard_v0_3/` (git-ignored). 171,948 COTS culled on 72 domain reefs.
  - Calibration reefs: North Direction 2013 had 417 culled (about 20/ha, at its observed second peak); Eyrie 2021–24 had 1–3/ha a year; Lizard and MacGillivray had little.
  - The heaviest culling was on other reefs (Startle, Irene, Eye, U/N 14-133), mostly 2012–15.
- **Fifth-round decisions (lead):**
  - the working configuration is now **W** = North Direction share at min (run 32 `obs015_cv2_pnd`);
  - Lizard median share tested with both detectabilities;
  - cull design fixed (protocol);
  - replicates vary both seeds.
- **B3 / B4 replicates done** (run 37; re-scores 38–39; G2′ 40–41; `candidates_audit.json`).
  - Loss v2 over 5 seeds: W 6.79 ± 0.79, C3a (Lizard median) 6.72 ± 0.70, C5 (run-17 coral) 6.84 ± 0.49. They are indistinguishable.
  - C3a matches Pratchett's density but triples the Lizard trough.
  - C5 alone passes the troughs, but its coral fit is much worse.
  - Pratchett `q` is dominated.
- **C2 done** (COTSMod `apply_removals!`; ADRIA `COTS_CULL_REMOVALS_CSV`; run 42; re-score 43; G2′ 44; `culling_audit.json`).
  - The failure flag fired: in 2021–24 the model has no outbreak to cull (Eyrie 88% of culls unmet).
  - Where culls apply, the effects are small (North Direction 2013 peak 0.82 → 0.55).
- **Open for the lead** (see the latest protocol entries):
  - which candidates go forward;
  - how to treat the coral-fit vs trough trade-off;
  - whether culling goes into the simulations, and how to read the 2020s cull numbers (a missed 2020s outbreak?).
- **Remaining AGENTS gates:** held-out years and leave-one-reef-out (gate 4), ±25% sensitivity (gate 5).

**Phase A: lock the scoring (before any more model changes).**
- **A1. Frozen score v2.** Preregister and implement score v2 as a wrapper; `cots_cycle_metrics.jl` stays unchanged.
  - **North Direction:** use the raw observed series for peaks and loss (1995 and 2013 peaks).
  - **Eyrie:** drop the first-peak terms (peak count, timing and height for the missing observed first peak), keeping the 2013 peak and the series terms.
  - Report v1 and v2 side by side for every candidate.
  - **Natural place:** a `--score v2` option in `rescore_detectability.jl` (so stored trajectories are re-scored without re-simulation), plus a matching column in `evaluate_gate_g2prime.py`.
- **A2. Pratchett `q` re-score.**
  - **Re-score** run 25 (`cv2_pick`, `cv2_acrpoc`, `cv2_q50`) with `--q=0.030`, under score v2 and G2′.
  - **Consistent run.** Also run `cv2_pick` with seeds inverted at `q` 0.030. The `lizard_coral` runner pattern (`obs030`) shows how; it needs a small runner extension for the `initial_share` experiment.

**Phase B: candidate set and timing.**
- **B1. Diagnose late first waves on the working configuration.** Use stored fluxes (settlement, maturation, food survival, seeds) for North Direction and MacGillivray in 1992–98. Levers already in the code, with their earlier effects:
  - `Q` (run 11: `Q` 10 gives a mean first-peak year of 1997.0 against 1998.3 at `Q` 30);
  - seed level (EB vs IDW: run 17 shows little effect);
  - share levels (higher shares give later peaks: q50 / q75 push North Direction to 1999);
  - the Allee scale (held fixed by AGENTS; counterfactual only);
  - external supply timing.
- **B2. Bounded timing factorial on the working configuration**, scored with v2. Suggested: `Q` {10, 15, 20} × North Direction share {q25, min} (MacGillivray and Lizard fixed) × observation `q` {0.015, 0.030 consistent}. The reference arm replays `cv2_pick`.
- **B3. Candidate set** (target 3–6 named configurations), chosen by score v2, G2′, coral `J` and peak timing. Likely members:
  - **C1** working configuration;
  - **C2** C1 with Pratchett `q`;
  - **C3** C1 with Lizard at q50 share (matches Pratchett's density);
  - **C4** C1 with the best timing fix from B2;
  - **C5** the run-17 baseline (old coral) as a structural control.
  - Record each as a TOML/CSV parameter set under `sandbox/calibration/candidates/` (checked in), with run IDs, switches and hashes.
- **B4. Promotion gates** (AGENTS.md), for the candidate set only, labelled "calibration candidates for simulation", not the ADRIA default:
  1. **Tests:** legacy control replay, unit tests (COTSMod `Pkg.test`, `test_bleaching_floor.jl`, `test_connectivity_retention.jl`), parse test, data provenance (hashes).
  2. **Qualitative gate:** G2′ (with score v2).
  3. **Replicated seeds:** at least 5 recorded `ADRIA_RNG_SEED` and `COTS_CONNECTIVITY_SEED` combinations. Report peak timing, height, prominence and inter-peak period (mean ± sd).
  4. **Held-out tests:** temporal (fit 1992–2012, check 2013–24) and leave-one-reef-out; regional validation where observations exist.
  5. **Sensitivity:** ±25% on `Q`, μ, the low-DHW scale, Ext×5 and the shares.
  6. **Records:** run metadata complete (revisions, seeds, hashes, failures).
  7. **Documentation:** a README and protocol entry with the exact reproduction commands.

**Phase C: COTS removal (culling).**
- **C1. Data.**
  - **Source:** the CSIRO workbook (`sandbox/data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx`, sheet `Cull`: date, reef, site, bottom time, cohorts 1–4 by size). Cull dives in the Lizard region run 2012–26; e.g. Eyrie 2022 removed 1,592 COTS; Lizard Lagoon 2024 removed 566.
  - **Builder:** a checked-in `build_cull_removals.py` producing reef × year (and site where possible) removals per ha of reef area, split by cohort / size class. Record the cohort size definitions, units, hashes, the AMPTO/EOTR/legacy sources, and dive coverage. Output to the git-ignored `sandbox/data/derived/`.
- **C2. Mechanism.** An opt-in COTSMod switch applying observed removals each year:
  - absolute COTS/ha, bounded by the available adults and subadults, mapped by size class to model stages;
  - applied at a stated point in the annual step;
  - pure unit tests in COTSMod, and an ADRIA adapter env switch (e.g. `COTS_CULL_REMOVALS_CSV`) with an adapter test.
  - Preregister the timing, units and the size-to-stage mapping first.
- **C3. Experiment.** The candidate set × {no cull, observed cull}, 2012–2024: the effect on the 2013–15 second wave, troughs and coral. A no-cull reference replays each candidate exactly. Note that culling also biases the 2021–26 SALAD/manta detectability pairs.

**Open data notes for the next context.**
- **Survey timing.** LTMP reporting year = survey year + 1 for September–December surveys. Scoring keeps the reporting year (lead), and `derived/ltmp_survey_timing/survey_timing.csv` flags the bleaching-misaligned surveys. Read peak timing against survey dates.
- **Pratchett (2005) Lizard pairs:**

| Belt-transect survey | Density (COTS/ha) | LTMP manta (October survey) |
| --- | --- | --- |
| October 1996 | 40 | 0.90 COTS/tow |
| November 1997 | 43 | 1.92 |
| November 1998 | 13.5 | 0.05 |
| September 1999 | about 0.5 | 0.00 |

  Reef-wide coral was 32% (1996) and 23% (January 1999). The crash was disease-associated.

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
| 19 | G2′ evaluation of runs 11, 15 and 17 (re-scores 14 and 18) |
| 20 | Lizard local retention (negligible; Lizard food-limited) |
| 21 | detectability re-score of run 20 |
| 22 | G2′ evaluation of run 20 |

2026-10-14:

| Run | Result |
| --- | --- |
| 01 | detectability evidence v1 (pseudo-replicated pairing; superseded) |
| 02 | detectability evidence v2 (47 paired surveys; strong hyperdepletion; `A½` not identified) |
| 03 | bleaching floor × fast-share multiplier (floor 4 lifts Lizard; coral re-fit needed; baseline kept) |
| 04 | G2′ re-evaluation with the trough gated on Lizard and MacGillivray (renamed from a colliding 03) |
| 05 | detectability re-score of run 03 |
| 06 | G2′ evaluation of run 03 |
| 07 | survey-timing tracker |
| 08–11 | Pratchett-anchored observation re-scores (`q` 0.030; `q_max` 0.055 / `A½` 30) and their G2′ evaluations (08–09 renamed from unnumbered IDs) |
| 12 | coral re-fit: floor × μ × forcing `q` (floor 1 kept; μ at bound 0.08) |
| 13 | Lizard hard coral × `q` × μ (more Lizard coral lowers its wave; μ 0.08 improves both fits) |
| 14–17 | re-scores and G2′ of run 13 at `q` 0.015 and 0.030 |
| 18–19 | exploratory re-score (`q_max` 0.04, `A½` 10; post hoc, used for no decision) |
| 20 | coral re-fit v2 (μ identified at 0.10 / 0.12; cv2 selected) |
| 21 | 1991 fast-share ranges (healthy-year photo transects) |
| 22 | starting fast-share levels on cv2 (stage 2a) |
| 23–24 | re-score and G2′ of run 22 |
| 25 | starting fast-share confirmation with per-reef picks (stage 2b; loss 7.97) |
| 26–27 | re-score and G2′ of run 25 |
| 28–29 | score v2 re-score and G2′ of run 25 |
| 30–31 | score-only `q` 0.030 re-score and G2′ of run 25 |
| 32 | timing factorial (`Q` × North Direction share × consistent `q`) |
| 33–36 | re-scores (`q` 0.015 / 0.030) and G2′ of run 32 |
| 37 | candidate set × 5 seed pairs |
| 38–41 | re-scores and G2′ of run 37 |
| 42 | observed-culling screen |
| 43–44 | re-score and G2′ of run 42 |

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
