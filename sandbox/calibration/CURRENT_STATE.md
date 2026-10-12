# COTS calibration: current state (2026-10-12)

This is a short orientation page for resuming work. For full detail see:
- the [model log](MODEL_LOG.md) for the decision ledger;
- the [crash-mechanism protocol](crash_mechanism_protocol.md) for preregistrations and decisions;
- the [calibration README](README.md);
- the repository-root `AGENTS.md` for the rules.

Nothing described here is promoted to the default model.

## Where the code lives

| Piece | Location |
| --- | --- |
| COTS equations, switches and pure tests | **Separate repo** `..\COTSMod.jl` (`src/COTSMod.jl`, `test/runtests.jl`), origin `github.com/open-AIMS/COTSMod.jl`. Loaded by local path from `sandbox/Manifest.toml`, so uncommitted edits affect runs. |
| ADRIA adapter | `ADRIA/src/ecosystem/cots.jl` (`_cots_runtime_params`, initial-state reset, diagnostics) |
| Environment switches and logs | `ADRIA/src/scenario.jl`: `COTS_*` env vars. `cots_mechanism_log` has 7 channels: food memory, food survival, senescent adults, senescent deaths, preferred diet share, juvenile predation, maturation fraction. |
| 1991 experiment runner | `sandbox/calibration/run_1991_initialization_test.jl` (one `OWEN_1991_EXPERIMENT` per screen) |
| Audits | `sandbox/calibration/audit_1991_<experiment>.py` (reference-replay check, gate table, plot, JSON) |
| Derived habitat and cyclone inputs | `sandbox/domain_building/build_site_habitable_fraction.py`, `build_cyclone_hours_rrn.py`. Outputs go to the git-ignored `sandbox/data/Lizard_Historical_v0.2/derived/` because the sources are internal or licensed. |
| Run index | `python sandbox/calibration/build_model_run_inventory.py` |

Python environments:
- `rasterio` (GDAL 3.9, reads the GBR10 FileGDB rasters) is in the user Python.
- `openpyxl` (for the RRN workbook) is in base miniconda (`C:\Users\smatthew\AppData\Local\miniconda3\python.exe`).
- Conda GDAL 3.6 cannot read FileGDB *rasters*.

## Working model ("condition, no switching")

This is the project lead's preferred configuration. Its reference arm is `k05_nocyc` in [`20261012T_benthic_cyclone_1991_v2`](runs/20261012T_benthic_cyclone_1991_v2/metadata.toml), identical to `k05_bleaching` in the habitat run and `condition_ref` in the switching-repair run.

- **Grazing:** `COTS_PER_CAPITA_CONSUMPTION=true`, 10 m² per adult per year, fast prey first, ÷ site `k`.
- **Starvation:** `COTS_PREFERRED_PREY_STARVATION=true`, remembered fast cover, threshold `0.05`, memory `0.5` yr.
- **Senescence:** `COTS_ADULT_SENESCENCE=true`, ages 2–5 at archived `m3`, a 6+ class at `0.8`/yr.
- **Condition:** `COTS_PREFERRED_CONDITION=true`, condition from fast cover, saturation `0.10`, memory `1` yr.
- **Units:** `q = 0.015` COTS/tow per adult COTS/ha on every reef (approved convention). Archived parameters migrated by `k_unit = 19.61`. Allee scale `3` COTS/ha.
- **Reproducibility:** `ADRIA_RNG_SEED=90972102` must be pinned. ADRIA otherwise seeds the coral RNG from `floor(sum(params))`.
- **Forcing:** 1991 IDW CPUE/0.015 with an assumed 6:3:1 preload; Owen annual connectivity sampled with seed `20260930`; historical DHW; waves zero; domain `k = 0.5`; cyclones off unless stated.

Performance:
- **Lizard:** first peak 2001 at 1.40 COTS/tow; trough `0.103`; second wave 2014 at 0.57 (observed 0.74).
- **Across reefs:** 5/7 matched; two-wave reefs 2/3.
- **Fails on:** small second waves away from Lizard; late first peaks (about 2001 vs 1996–98); trough adults 4–10/ha vs observed about 0–3/ha.

## Site versus reef resolution

The model is **site-resolved** (2,895 V2 sub-reef sites on 113 reefs; COTS and coral state per site) and observations are compared after area-weighted aggregation to reefs. In practice it **behaves reef-level**:
- COTS larval supply uses reef×reef Owen matrices divided evenly across a reef's sites, and the external boundary is per reef;
- initial coral is per RME reef, DHW is per reef, and `k = 0.5` is uniform.

Within Eyrie, MacGillivray and North Direction, all sites follow identical trajectories. "Lizard Island Reef" differs only between its three GBRMPA reefs (14-116b/c/d). Site-level habitable fraction from the benthic map is the first input that varies within reefs.

## What has been established

1. **Legacy structure.** Starvation on total cover, linear grazing and a single adult class produced a plateau. Fast coral was about 0 from 1992.
2. **Units.** The fitted scales implied 1 model unit ≈ 20 COTS/ha; this is now fixed at `q=0.015`. The archived 1997 match relied on an inconsistent preload.
3. **RNG coupling.** Parameter edits changed the coral random stream; seeds are now pinned.
4. **Mechanisms.** Per-capita grazing + starvation + senescence gives a boom, starvation crash and coral recovery. Preferred-food condition cuts trough larval supply about 70%.
5. **Prey switching** is retired: it stabilises a plateau, even with diet-linked starvation or reduced feeding on massives.
6. **Habitat.** RME `k` (mean 0.715) and benthic site `k` (total mean 0.435; mapped 0.625; geomorphic-infilled 0.58) improve site-area coral against LTMP. But they enlarge first outbreaks and trough adults, so second waves vanish. No-bleaching does the same.
7. **Cyclones.** The RRN `cyc_maxHrs4mw` hindcast (Ita 2014 about 14 h, Nathan 2015 about 50–74 h), with mortality `1-exp(-0.00444·h)` from the COTS-Coral-OutcomesV2 LTMP dose-response, improves 2015–16 coral and barely affects COTS.
8. **Trough decomposition.** Trough adults are mostly new each year, from settlement of 15–60/ha/yr (70–87% internal immigration) through a coral-independent recruit stage. The 3/ha Allee scale is weak at 4–10/ha.
9. **First-wave timing.** The observed Lizard rise is about ×2.5/yr (0.23→1.92 COTS/tow over 1994–98); the model's is about ×1.5/yr. Early settlement is about 7 vs about 35+/ha/yr needed. Candidates: a weak 1990s external pulse, or low manta detectability of the 1991 initial stock.

## Latest result: collapse × fast-coral growth (2026-10-12)

**Leading candidate (unchanged):** working model + **External ×5** (assumption) + predator pit **`Q = 30`** (tuned) + gradual Allee, on `Lizard_Historical_v0.3`, with CoralBlox calibration group 1. Runner arm: `predation_capacity` / `ext5_PA_q30`, also replayed as `outbreak_collapse` / `g1_none`.

- **[Q tuning](runs/20261012T_predation_capacity_1991_v1/metadata.toml).** `Q = 30` gives 5/7 matched and second waves on Lizard, Eyrie and North Direction. Troughs stay at 5–10 adults/ha.
- **[Collapse × coral](runs/20261012T_outbreak_collapse_1991_v1/metadata.toml):**
  - **ReefMod-style die-off** (`COTS_OUTBREAK_COLLAPSE`): an abrupt drop, then a rebound. Second waves are lost because the crash no longer exhausts coral. Not recommended.
  - **Coral calibration group.** Lizard runs use group 1, but the RME layer puts all four calibration reefs in group 5. Group 5 raises first peaks to observed levels but overpredicts coral.
  - **Fast linear extension ×1.2** gives second waves on all four reefs (6/7 matched), including MacGillivray for the first time. They are too large and late (up to 2.5 COTS/tow, 2016–17).
- **Still failing:** the trough floor; first peaks around 2000–01; second-wave amplitude and timing when fast growth is raised.

**Next (project lead to choose):**
1. Fit fast-coral growth to observed coral recovery, ideally with group-level cover data, then rerun the COTS screen.
2. A faster starvation-coupled crash (shorter preferred-food memory or a sharper response).
3. Decide on the calibration group (group 1 vs RME group 5; a v0.4 rebuild).

## Known issues and open decisions

- **ADRIA cannot run `k = 0` sites.** `scenario.jl:850` uses `length(habitable_locs)`, and `constrain_recruitment!` mixes all-site and habitable-site vectors. Runs floor `k` at `0.01` and zero those sites' initial cover. The upstream fix is not made.
- **Labelled assumptions (project lead, 2026-10-12):** External ×5 (a magnitude on an unidentified production term; timing stays evidence-derived) and `Q = 30` (tuned; no field bound on the model's effective juvenile scale).
- Before any optimizer: a qualitative gate pass, then replicated seeds, `q` sensitivity and held-out reefs/years.
- **ADRIA: non-contiguous `CB_CALIB_GROUPS` are misindexed.** The sampler drops absent groups and `scenario.jl` indexes by group number. The `outbreak_collapse` runner relabels RME groups contiguously.
- **The Lizard loader defaults `CB_CALIB_GROUPS` to 1.** The RME layer assigns the domain reefs to groups 5 (58 reefs), 1 (44) and 4 (11).

## How to add the next screen

1. **Preregister.** Append to `crash_mechanism_protocol.md`: hypothesis, equations, evidence-backed settings, arms including an **exact-replay reference arm**, gate, failures, rollback.
2. **Implement and test.** Add the opt-in COTSMod switch, with a default that reproduces current behaviour, plus tests (`julia --project=. -e "using Pkg; Pkg.test()"` in `..\COTSMod.jl`).
3. **Wire it into ADRIA.** Add it to `_cots_runtime_params` and to the env parsing and validation in `scenario.jl`. Add an adapter test in `sandbox/calibration/test_1991_initial_state.jl`.
4. **Add the experiment.** Add it to the runner's `EXPERIMENT` list and treatment list, and to the `MIGRATION`-style blocks (`q`, `k_unit`, RNG pin, per-capita / H1 / H2 / condition flags). Avoid top-level variable names that shadow Base, such as `values`. Run `julia --project=sandbox sandbox/calibration/test_calibration_scripts_parse.jl`.
5. **Run.** Set `OWEN_1991_EXPERIMENT` and a new `OWEN_1991_RUN_ID` (runs never overwrite). About 2–4 minutes after domain load.
6. **Audit.** Copy the latest `audit_1991_*.py`, run it, and inspect the PNG.
7. **Log.** Add a MODEL_LOG row, a protocol decision and a README section, then rerun the inventory script. Update this page.
