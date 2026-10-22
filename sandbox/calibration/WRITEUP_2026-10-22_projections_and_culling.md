# Lizard COTS model: calibration, projections and culling scenarios

*Write-up of the work from 19 to 22 October 2026: rounds 29–33 of `crash_mechanism_protocol.md`. Nothing here is promoted to the ADRIA default model. All runs replay their controls exactly unless stated otherwise.*

## 1. The model used throughout

The working configuration is **W-wave-lf04**:
- the W-wave calibration of the COTS model coupled to ADRIA and CoralBlox on the Lizard domain (113 reefs, 2,895 sites);
- plus a 1991 starting increase of 0.04 (4 percentage points of habitable area) in fast-growing coral at Lizard Island Reef.

The increase was adopted on 19 October. Its basis:
- it brings Lizard's first outbreak (1997) to 36 adults/ha, closer to the 40–55 adults/ha that Pratchett (2005) counted on belt transects;
- it does not worsen the calibration score;
- shallow snorkel surveys (CSIRO RHIS) show about 1.6 times more branching and plating coral at Lizard than the 9 m photo transects.

**Main features of the configuration:**
- **Outside supply:** most larval supply comes from reefs outside the domain: 60–90% of settlement at the calibration reefs.
  - It is driven by an outbreak-probability map for every GBR reef. For 2018–2026 the map is corrected with diver and cull observations.
  - It is multiplied by 5.
- **Observation:** manta-tow counts are compared through a fixed conversion of 0.015 COTS per tow per adult/ha. A density-dependent alternative fitted to paired diver–manta data was tested and rejected, because it misreads the first outbreaks at North Direction and Eyrie by 5–30 times.
- **Coral:** calibrated CoralBlox growth, with background mortality and bleaching. Natural heat-tolerance adaptation (heritability 0.3) is part of the model.

## 2. How well the model reproduces the observations (1992–2026)

The figures are in `sandbox/data/derived/calibration_report_figures/2026-10-22_wwave_lf04_fit_statistics/`.
- Each panel shows the model (mean of 5 random seeds, with the seed range and a ±25% parameter envelope) against the observations.
- Every panel carries the number of matched years (n), the bias (model minus observed), the mean absolute error (MAE), and the Pearson (r) and Spearman (ρ) correlations.

**COTS, LTMP manta tows (COTS per tow), culling off:**

| Reef | n | Bias | MAE | Pearson r | Spearman ρ |
| --- | --- | --- | --- | --- | --- |
| Lizard Island | 28 | −0.06 | 0.22 | 0.70 | 0.53 |
| MacGillivray | 31 | +0.01 | 0.14 | 0.48 | 0.32 |
| North Direction | 28 | +0.14 | 0.17 | 0.58 | 0.63 |
| Eyrie | 10 | +0.13 | 0.26 | 0.04 | 0.58 |

**Coral, LTMP manta hard coral (fraction):**

| Reef | n | Bias | MAE | Pearson r | Spearman ρ |
| --- | --- | --- | --- | --- | --- |
| Lizard Island | 28 | +0.14 | 0.14 | 0.56 | 0.60 |
| MacGillivray | 31 | +0.02 | 0.05 | 0.73 | 0.44 |
| North Direction | 28 | −0.07 | 0.11 | 0.27 | 0.10 |
| Eyrie | 10 | +0.06 | 0.08 | 0.51 | 0.35 |

**Reading these statistics:**
- **Pearson r and MAE** are dominated by a few outbreak years; one year of timing offset lowers r sharply. That is why Eyrie has r 0.04 but ρ 0.58.
- **Spearman ρ** is more robust, but the many zero-count years tie.
- **Coral bias is partly structural.** The model reports coral as a fraction of habitable reef area; the LTMP data are from 9 m transects or manta tows. The +0.14 at Lizard is largely this difference.
- **Extra columns.** The file `fit_statistics.csv` adds relative bias, Nash–Sutcliffe efficiency (skill against simply predicting the observed mean), error on a log scale, and a true skill statistic for outbreak state (outbreak or not at 0.22 COTS per tow).
- **Recommendation:** for the report, quote the four statistics on the figures plus outbreak-state skill and peak timing (`peak_timing_summary.csv`). The latter two measure what the calibration targets.

**Known weaknesses:**
- The model's 1990s outbreak at Lizard is still lower than observed (0.54 against 1.92 COTS per tow).
- The model overshoots the 2026 densities: Lizard 23.6 adults/ha against 5.0 counted by divers; North Direction 30.1 against 7.2.
- The overshoot comes from the outside supply in 2024–26.

## 3. Projection set-up (1991–2100)

- **Heat stress:** annual maximum degree heating weeks (DHW) from 10 CMIP6 climate models (RME data set).
  - Each model is bias-corrected per reef against observed DHW for 2000–2026, using quantile delta mapping. This keeps each model's warming trend.
  - Uncorrected, the climate models were between 0.9 and 4 times too hot in 2000–2026.
- **Cyclones and the chlorophyll larval-survival signal:** whole observed years (1993–2024) are drawn at random for each future year, the same year for all reefs.
- **Outside larval supply:** observed 1991–2026 outbreak years are resampled in 15-year blocks.
  - Supply is scaled down as domain coral cover declines, because outside reefs would lose coral under the same climate.
- **Culling:** observed removals to 2026, then forward culling driven by diver effort.
  - Catch per diver-hour rises steeply with density (from the diver–cull relationship D = 7.874 × CPUE^0.549).
  - Recent effort is 3,733 diver-hours per year across 36 domain reefs (mean of 2023–2025).
- **Optional switches** (tested as scenarios):
  - a 2026 reset of COTS densities to diver counts;
  - extra thermal acclimation;
  - larvae of coral from outside the domain.

## 4. Projections to 2100: what drives the outcome

**Ensemble:** 5 climate models × SSP2-4.5 and SSP3-7.0 × 5 random seeds × 9 culling levels × with or without the 2026 reset (900 runs).

1. **Climate dominates.**
   - Under SSP2-4.5, mean coral cover in 2051–2075 ranges from 0.01 (UKESM) to 0.16 (EC-Earth3), depending on the climate model.
   - Under SSP3-7.0, all five models fall below 0.08 by 2051–2075 and to 0.02 or less by 2076–2100.
   - Natural adaptation is already included; additional acclimation of 0.12 DHW per year delays the collapse by about 5 years.
2. **Culling reduces outbreaks, but has little effect on coral after 2050.** Medians over climate models and seeds, mean of the four calibration reefs, SSP2-4.5:

   | | No forward culling | Recent effort ×1 | Recent effort ×3 |
   | --- | --- | --- | --- |
   | Years in outbreak, 2027–50 | 8.75 | 7.5 | 5.5 |
   | Peak adults/ha, 2027–31 | 63 | 51 | 36 |
   | Mean coral cover, 2027–50 | 0.160 | 0.169 | 0.174 |
   | Mean coral cover, 2051–75 | 0.100 | 0.101 | 0.102 |

3. **Targeted culling.** Sending the same hours to reefs currently in outbreak removes about 2.5 times more COTS across the domain, but helps the four calibration reefs less. Its domain-wide benefit is not yet measured, because the model only reports the calibration and a few extra reefs.
4. **Larvae of coral from outside the domain made no difference.** Recruits are killed by heat before they add cover.

The 2050 figures are in `.../2026-10-21_wwave_lf04/projections_2050/`.

## 5. The 2020s outbreak: how much culling would have prevented it?

**Question.** How much extra culling effort, above what actually happened, would have been needed:
- in 2020–2025, to suppress the build-up; or
- in 2026–2030, to stop COTS exceeding the outbreak threshold?

**Design.**
- **Reefs:** success means no year above the threshold on Lizard, MacGillivray, North Direction, Eyrie and Martin, in at least 4 of 5 seeds.
- **Thresholds:** 14.7 adults/ha (equivalent to 0.22 COTS per manta tow) and 3 adults/ha.
- **Outside supply settings:** three, because the current one overshoots in 2026:
  - **current** (diver/cull-corrected);
  - **graded 0.5** (a weaker correction that tracks Lizard to 2025);
  - **original** (no correction; no overshoot, but also no 2020s build-up).

**Results** (figures in `.../2026-10-22_buildup_2020s/`):
1. **The outbreak threshold is not crossed during 2020–2025** on these five reefs, even without extra culling. The highest density in those years is 6–12 adults/ha.
2. **The breach happens in 2026 and later.** Juveniles that settled from outside supply in 2024–2026 mature then.
3. **Extra culling of adults in 2020–25 cannot stop the 2026 surge.** At 32 times recent effort (about 120,000 extra diver-hours per year), 2020–25 densities fall to 3–5 adults/ha, but 2026 still reaches 11–15 adults/ha in the two corrected settings. Culling adults on these reefs does not reach the juveniles, or the outside reefs that produce them.
4. **Stopping breaches in 2027–2030 needs between 16 and 32 times recent effort**, about 60,000 to 120,000 extra diver-hours per year.
   - We tested 16 times and 32 times but nothing in between, so for some settings the true requirement lies somewhere between these two values.
   - In the corrected settings this requirement is inflated by the model's 2026 overshoot.
5. **The stricter threshold of 3 adults/ha was never reached at any tested effort.** Catch per hour collapses at low density, and outside supply keeps arriving.

## 6. Culling on reefs outside the domain

The earlier runs did **not** reduce outside supply for COTS culled on outside reefs. We added that discount:
- for each outside reef and year, the share of its adult population removed = COTS culled ÷ (cull-derived density × reef habitable area), capped at 90%;
- that share is subtracted from its larval supply to the domain.

The reconstruction of the supply reproduces the original exactly before the discount is applied.

**Four scenarios** (figures in `.../2026-10-22_outside_cull_2020s/`):
- observed outside culls to 2025, none after;
- recent outside effort carried on through 2026–2030;
- outside effort doubled in 2020–2030;
- outside effort doubled in 2026–2030 only.

**Result.** At the scale recorded in the cull data, outside culling changes the 2020s build-up by at most 4%.
- A typical culled outside reef loses about 1% of its estimated adult population per year.
- The few reefs with large removals (40–90%) are far south of Lizard and supply little of its larvae.
- Even doubled outside effort lowers the larval supply reaching Lizard, Eyrie and MacGillivray by at most 4% in any year.
- The extra inside effort needed is unchanged, except one case (graded setting, early targeted culling: 16 instead of 32 times recent effort).

**Caveat.** The share removed is probably underestimated, because it divides by the whole habitable reef area while culling happens on COTS-suitable habitat.

## 7. Limitations

1. **The 2026 overshoot** in the corrected supply settings inflates the projected 2026–2030 outbreak and the effort needed to stop it.
2. **COTS density and area.** Densities are calibrated to diver and manta data from COTS-suitable habitat, but the model counts animals over the whole site or reef area. If only a share s of the area is suitable, every cull (inside, forward and outside) removes too small a share of the population, by a factor 1/s. This is the lead's next proposal (section 8).
3. **No juvenile culling.** Juveniles are not culled, and outside reefs are not modelled (only their outbreak probabilities).
4. **Stationary outside supply.** Projected outside supply is resampled history, scaled by domain coral; it does not respond to culling or climate at the source reefs.
5. **Extrapolation.** The catch-per-hour relationship is extrapolated above about 20 adults/ha and capped at 20 COTS per hour.
6. **Licensed data.** The SALAD diver densities, cull records, RHIS records and the RRN cyclone data are licensed. Derived files stay in the git-ignored `sandbox/data/`; check permissions before publishing figures 2 and the culling figures.

## 8. Next steps

1. **COTS-suitable habitat basis** (lead's proposal):
   1. Compute the suitable share of each site from the lead's habitat layers, and how much larger the share culled would have been.
   2. Make culling and the outside discount use suitable area, as an opt-in switch.
   3. Assign cull removals to the site where the dive took place. Where the model holds too few COTS there, pass the remainder to the nearest sites with COTS, within the reef and a distance limit.
   4. Then decide whether larval export should move to the same area basis; that needs a calibration re-check.
2. **Agreed but not yet implemented:**
   - a downward-only 2026 reset;
   - thermal-acclimation scenarios in the ensemble;
   - domain-wide outputs, so that targeted culling can be judged across all reefs.
3. **Longer term:** an outside boundary from ReefMod-GBR, or a GBR-wide calibration, so that culling and climate act on the source reefs too.

## 9. Reproduction and files

- **Protocol and decisions:** `sandbox/calibration/crash_mechanism_protocol.md`, from "Lead's decision (2026-10-19, twenty-ninth round)".
- **Run ledger:** `MODEL_LOG.md`.
- **Commands:** the last three sections of `README.md`.
- **Key runs:**

  | Purpose | Runs |
  | --- | --- |
  | W-wave-lf04 control | 20261019T_51 |
  | Culling control | 20261021T_01 |
  | Scenario screen | 20261021T_03 |
  | Ensemble | 20261021T_04–08 (audit 09) |
  | Report ensembles | 20261021T_10/11 |
  | 2020s build-up | 20261022T_01–03 |
  | Outside culling | 20261022T_04–15 |

- **Builders:**
  - `sandbox/domain_building/build_projection_inputs.py`
  - `build_outside_cull_discount.py`
  - `sandbox/calibration/build_cull_effort.py`
  - `build_reset_targets.py`
- **Tests:**
  - `test_projection_inputs.jl` (179 checks)
  - `test_projection_mechanisms.jl` (35 checks)
- **Julia:** run with `--pkgimages=no`. The ADRIA compiled-package file in `~/.julia/compiled` became locked during this work.
