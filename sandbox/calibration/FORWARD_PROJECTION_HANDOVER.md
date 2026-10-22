# Forward projection to 2045 (extended to 2100): handover (2026-10-19)

**Update (round 29):** the lead extended the horizon to **2100**, to align with other forecasts. Read "2045" below as "2100". The projection design is preregistered in `crash_mechanism_protocol.md`, "Preregistered: projection design to 2100".

This note hands the Lizard COTS–coral calibration over to a new chat whose job is forward projection to 2045. Read it together with:
- repository-root `AGENTS.md`: the working rules (preregister, opt-in switches, exact replay, never overwrite runs);
- `CURRENT_STATE.md`;
- the end of `crash_mechanism_protocol.md`: all decisions to date.

The project lead makes the modelling decisions. Ask them in plain ecological language (current vs proposed, why, what is still wrong, recommendation), with numbered items.

## 1. What is calibrated (the starting point)

**Working configuration: W-wave-lf04** (lead, 2026-10-19, round 29): W-wave plus Lizard 1991 fast cover +0.04. It is not promoted to the ADRIA default, because gates 4 and 5 were only partly passed (on W-wave).
- Arm `wc_b100_q25_eb_h10_obs015_s<k>_t60_e50_lf04_aprmak_cv2_pnd_wv025`.
- Replay control: run 20261019T_51 (its +0.04 arms equal run 48, its 0 arms equal run 15).

| Component | Setting | Where |
| --- | --- | --- |
| Domain | Lizard v0.4: 2,895 sites, 113 reefs, 1985–2026; runs start in 1991 | `sandbox/data/Lizard_Historical_v0.4` |
| COTS | Separated recruitment; gradual Allee `exp(-(3/A)^0.5)`; juvenile predator pit `Q·J/(J+Q)` with Q 25 COTS/ha/yr; 1991 start from the 1989–91 surface inverted at A½ 10 (`h10`); fecundity 37.9 larvae/adult/yr pre-pelagic; observation q 0.015 COTS/tow per adult/ha | COTSMod.jl (sibling repo) and the ADRIA adapter `ADRIA/src/ecosystem/cots.jl` |
| Internal dispersal | Owen hydrodynamic COTS matrices 2018–22, sampled by seed | dataset `owen_global_2018_2023_corr_w3` |
| External supply | `B = Σ_j C[j,i]·A_j/A_i·s_j·P_j(y)` over non-domain reefs × production 2,798 × Ext 5. `P_j(y)` is the outbreak-probability layer (`gbrPredsAdj`), with the diver/cull correction counted 3 times (`corr_w3`). **This is 60–90% of settlement.** | `sandbox/domain_building/build_outbreak_layer_face_value.py`, `build_lizard_cots_connectivity_owen.py` |
| Chl larval survival | p90 anomaly index (`aprmak`, b 3, 1990s corrected), domain sources only | `COTS_LARVAL_SURVIVAL_CSV` |
| Wave exposure | Settlement × ramp on reef-mean wave period, 6.6 → 7.0 s, floor 0.25 | `COTS_SETTLEMENT_MODIFIER_CSV`, table `wave_settlement_modifier/tp66_70_f025` |
| Coral | CoralBlox, β 1, bleaching ×0.25 below 6 DHW, chronic μ 0.12 | runner arm `wc_b100_..._cv2_pnd` |
| Culling | Reef-proportional allocation of observed culls (`COTS_CULL_ALLOCATION=reef_proportional`). **Off in calibration scoring**; available for simulation. | `build_cull_removals.py`, `removals_sigma210.csv` |

**Example arm:** `wc_b100_q25_eb_h10_obs015_s<k>_t60_e50_aprmak_cv2_pnd_wv025`. Runner: `sandbox/calibration/run_1991_initialization_test.jl`, experiment `wave` or `sens` (`OWEN_1991_SENS_WAVE=true`), with:
- `OWEN_1991_DOMAIN_VERSION=v0.4`;
- `OWEN_1991_COTS_DATASET=owen_global_2018_2023_corr_w3`;
- `OWEN_1991_EXTRA_REEFS=true`.

Reference runs: 20261019T_15 (W-wave, 5 seeds) and 20261019T_46 (W-wave ±25% ensemble, report figures).

**Labelled alternatives for simulation:**
- graded correction weight 0.5 (`owen_global_2018_2023_corr_w05_graded`, culling off): Lizard better, Eyrie too low;
- graded weight 2 (`_corr_w2_graded`, culling on): holds up best with culling; Eyrie fits; Lizard, MacGillivray and North Direction overshoot.

**Known limitations, which carry into any projection:**
1. The 2020s build-up is about 2–3 years late. The reason is the external layer: source-reef outbreak probabilities are 10–40× lower in 2018–21 than in 2014.
2. 2026 overshoot: model 22–30 adults/ha at Lizard and North Direction against observed 5–7.
3. One regional supply cannot give Eyrie a large 2020s outbreak and MacGillivray a small one. That contrast is reef-level and not yet explained.
4. Martin's 2014 outbreak is under-predicted (0.15 against an observed 0.51 COTS/tow).
5. Gate 5 (±25%) fails 2/12: model behaviour depends on external supply at about ×5 and on chronic coral mortality of at least about 0.12.

## 1a. The most conspicuous gap: Lizard's outbreak peaks are too low, and its coral too high

**Symptom.**
- **COTS:** W-wave's Lizard 1990s peak is about 25 adults/ha in 1997 (seed mean; 22–30 across the ±25% ensemble), shown as 0.38 manta COTS/tow. Observed: 1.92 COTS/tow in reporting year 1998 (surveyed October 1997) and 0.74 in 2013.
- **Coral:** with too few COTS, model coral at Lizard stays at about 0.30 of habitable area through the first outbreak, while LTMP 9 m manta and photo transects show 0.08–0.15.
- This feeds into projections: low Lizard outbreaks mean too little coral loss on the main calibration reef.

**What the target really is.** Two numbers have been used, and they differ by about 2.5×.
- With the frozen manta conversion (q 0.015 COTS/tow per adult/ha), 1.92 COTS/tow means about 128 adults/ha.
- Pratchett (2005) belt transects at Lizard measured 40–55 COTS/ha in 1996–98 (12.5 in October 1994), and reef-wide hard coral of 32% (1996) falling to 23% (January 1999). The crash was disease-like, not starvation.
- Paired manta and diver data give an effective q of 0.022–0.045 at outbreak densities, but 0.002–0.006 at 2020s densities. Manta detection is strongly density-dependent.
- Against the diver density, the model's deficit is about **1.6–2×** (25 against 40–55), not 5×.
- The model's coral level (0.30) is close to Pratchett's reef-wide 0.32. The 9 m LTMP sites sit lower (0.11–0.13 in 1996), so part of the "inflated coral" is a difference in support. What the model misses is the outbreak-driven decline (Pratchett: 32% → 23%).

**Diagnosis** (W-wave, run 20261019T_15, 1992–2001, Lizard against North Direction, which reaches about 60 adults/ha):
1. **Food.** Model fast coral at Lizard is only 0.04–0.05 of habitable area before the outbreak and falls to 0.014 by 1998. Food survival drops to 0.37 (1998), and body condition never exceeds 0.68 (North Direction 0.90). Larval output scales with condition², so Lizard produces 150–200 larvae/ha/yr at its peak against North Direction's 570–610.
   - Lizard's starting fast-coral share is set at its plausible **minimum** (the W share set), because the frozen manta loss penalised the higher troughs that larger first waves bring.
2. **Supply.** External settlement at Lizard (62–110/ha/yr in 1994–2000) is about two-thirds of North Direction's (95–158). Local retention is tiny (Owen diagonal 0.06–0.33%). Lizard's own larvae add little: retention up to ρ 0.9 raised the peak by only 3%.
3. **The predator pit** removes 12–15 juveniles/ha/yr at Lizard. That is a larger share of Lizard's smaller cohorts.
4. **Diet.** Model COTS eat preferred (fast) coral. When it runs out they starve rather than switch to slow coral, so coral stays high. Prey switching was tested (2026-10-12) and retired because the refuge stabilised the boom–bust cycle.
5. **Not yet resolved: the "food-scale" question** (protocol, local-retention decision, 2026-10-13, option (a)). Observed Lizard reached outbreak densities on 9 m fast cover of only about 0.02–0.04. The likely reason is that the prey pool is broader than *Acropora* + Pocilloporidae at 9 m: *Montipora*, *Isopora*, and crest or shallow *Acropora*. The starvation threshold θ 0.05 and condition saturation 0.10 are defined on the narrow 9 m pool. This was never tested; it needs evidence (e.g. RHIS benthic records in the CSIRO workbook, or Pratchett's reef-wide composition), not a fit.

**What has been tried** (from `crash_mechanism_protocol.md`; peak in adults/ha unless stated):

| Lever | Effect on Lizard's first wave | Why not adopted |
| --- | --- | --- |
| Local retention ρ 0.25–0.9 (2026-10-13) | +3% | Lizard's own production is too small to matter |
| Starting state (empirical-Bayes surface; `h10` start) | ≈ 0 | External supply dominates the build-up |
| Starting fast share: Lizard median instead of minimum (2026-10-14, cv2) | **29 → 53** (inside Pratchett's 40–55) | The frozen q 0.015 loss penalised the higher troughs that came with it |
| Starting share ×1.25 (W-wave ±25%, run 46) | 25.5 → 29.7 (+16%) | v3 +0.6 |
| Bleaching floor 4 DHW (old coral) | 29 → 39 | Coral fit 13–16% worse; replaced by ×0.25 below 6 DHW (adopted) |
| β 1.0 instead of 1.5 (adopted) | 0.35 → 0.41 COTS/tow | Small |
| Uniform food multiplier m 2–3 | 107–148 (the manta-scale target) | North Direction and Eyrie overshoot heavily; troughs explode |
| Larger habitable area (H-k) | 1.93 COTS/tow, but in 2002 | Wrong timing; doubled residual adults |
| Local larval output ×2 / ×3 (run 20261015T_45) | +14–23% (27–29) | Not evidence-based |
| Chl survival on external supply (run 27) | 25.5 → 29.3 (+15%) | Diver error worse (2026 overshoot) |
| External ×0.75 / ×1.25; pit Q ×0.75 / ×1.25; chl slope; low-DHW scale; μ (run 46) | 23–27 | Insensitive |
| Observation q 0.030 (Pratchett pairs, re-score only) | Reads 2× higher (Lizard 0.69 smoothed) | Troughs read 2× too; fails G2′. A linear q cannot fit both peaks and troughs |
| Low-density detection `A½` 10–30 | Peaks read lower | Fixes troughs, worsens peak reading |

**Latest test (2026-10-19, run 20261019T_48; lead's hypothesis that 1991 fast coral was underestimated).** Raising Lizard's 1991 fast cover by an absolute +0.02 / +0.04 (slow coral unchanged) gives:
- Lizard's first wave: 25.5 → 30.9 / 36.3 adults/ha;
- v3: 11.011 → 10.915 / 11.010 (G2′ 5/5);
- Lizard coral in 1999: 0.31 → 0.30 / 0.28 (Pratchett 0.23);
- other reefs: unchanged.

Both meet the fit rule. **Adopted: +0.04** (lead, round 29). It changes the 1991 start of every projection run.

**Round-29 follow-up** (run 20261019T_51; protocol "Round-29 tests: results"):
- The density-dependent (power-law) manta observation selects no larger share. Fitted to LTMP pairs, it misreads North Direction's and Eyrie's first waves by 5–30×, so it cannot serve as a regional observation model.
- +0.06 reaches Pratchett's density (42.6 adults/ha) but fails the G2′ trough test.
- RHIS shallow fast cover at Lizard is 1.6× the 9 m photo cover, which supports +0.04. Other reefs show 2–4×, so the under-sampling is general.
- Lizard's coral now falls 0.305 → 0.280 (1996–99), against Pratchett's 0.32 → 0.23.

**Settings that do raise the peak.** Two are evidence-compatible:
- Lizard's starting fast share at the median of its plausible range (53 adults/ha);
- less low-DHW bleaching before the outbreak.

Both were held back by the frozen loss, built on linear q 0.015, which reads the higher troughs as failures.

**Recommended next tests** (preregister; lead to decide):
1. **Score with a density-dependent manta observation model** (power or saturating, from the paired SALAD and Pratchett data) instead of linear q 0.015. Then re-select Lizard's starting share. This is the cheapest test and it targets the coupling found on 2026-10-14.
2. **Food scale:** define preferred prey on a broader, evidence-based pool (reef-wide or shallow fast cover, *Montipora*, *Isopora*) for θ and condition, using RHIS or Pratchett composition data.
3. **An outbreak-driven coral decline needs more COTS, not prey switching.** If tests 1–2 raise Lizard to 40–55/ha, check coral against Pratchett's 32% → 23%.

## 2. What a 2045 projection needs

| Need | Status now | What to do |
| --- | --- | --- |
| **Climate (DHW) 2027–45** | The Lizard v0.4 domain stops in 2026; its `dhw_RCP45.nc` is historical only, despite the name | GBR-wide CMIP6 annual max DHW exist locally: `sandbox/data/rme_ml_2025_06_05/data_files/dhw_csv` (3,806 reefs × 2000–2100; several GCMs × SSP1-1.9 … SSP5-8.5). Map them to Lizard sites (reef id, then sites; check against NOAA 2016–26 for bias) and build a domain v0.5 with a scenario dimension. The Moore domain (`Moore_2025-08-22_v080/DHWs/dhwRCP*.nc`) shows ADRIA's projected-DHW layout. |
| **Cyclones 2027–45** | Historical RRN mortality to 2024; 2025–26 zero | `rme_ml_2025_06_05/data_files/cyc_csv` has 100 GBR-wide cyclone simulations. Map them to reefs and convert to CoralBlox mortality as for the historical series. |
| **Chl larval survival 2027–45** | Index ends in 2026 | No observations. Options: block-resample historical years (seeded), climatology (multiplier 1), or water-quality scenarios (`water_csv/COTS_LARVAL_REDUCTION_*`, RME). |
| **Internal connectivity** | 2018–22 matrices sampled per seed | Keep resampling (record the sampled year); no new data needed. |
| **External larval supply 2027–45** | `P_j(y)` ends in 2026 | **The key decision**; see section 3. |
| **Starting state in 2026** | Calibrated runs end in 2026 with an overshoot | Option A: run continuously from 1991 to 2045 (simplest, keeps state consistency). Option B: re-initialise in 2026 from observations (adjust 2026 adults to the diver/cull densities). B is a lead decision, and it is a data-assimilation step that should be logged as such. |
| **Culling 2027–45** | Only observed removals (to 2026) | A forward cull scenario needs effort (dives or hours per reef per year), not removals. That means a new adapter: removals = f(effort, model density), using the SALAD–cull anchor `D = 7.874 × CPUE^0.549` inverted. Scenarios: none / recent effort / targeted. |
| **Coral under future heat** | CoralBlox with the calibrated settings | Decide on thermal adaptation and acclimation settings for projections (ADRIA has `a_adapt` / tolerance options); keep the calibrated values as the default. |
| **Runner** | `run_1991_initialization_test.jl` is a calibration runner | Write a projection runner (`run_projection.jl`): reuse the W-wave set-up, extend the time axis, add GCM / SSP / supply / cull scenario factors, and use the same output tables (`trajectories`, `flows`, `mechanisms`, extra reefs). |
| **Uncertainty design** | Seeds and a ±25% one-at-a-time ensemble | Projection ensemble: GCM × SSP × external-supply scenario × cull scenario × seeds (≥ 10) × a small parameter set (the ±25% ensemble members, or a Latin hypercube within those bounds). |
| **Metrics** | Calibration scores | Outbreak frequency, duration and peak density per reef; years ≥ 0.22 COTS/tow; coral cover trajectories and time below thresholds; culling benefit (coral saved per COTS removed). |

## 3. External larval supply after 2026: options

Settlement at the calibration reefs is mostly external supply, and its timing comes from `P_j(y)`. Projected outbreak waves will mostly come from how this boundary is projected.

| # | Option | What it gives | Effort | Main risks |
| --- | --- | --- | --- | --- |
| 1 | **Calibrate a GBR-wide reef-level version** (ADRIA reef-level GBR domain, the same COTSMod) so supply is endogenous; Lizard is then a nested sub-domain or a subset | Supply responds to climate, coral and culling everywhere; consistent physics | Large (months): GBR-wide calibration on LTMP manta (about 100+ reefs), outbreak initiation in the northern GBR has to emerge, compute | Same structural issues at larger scale (initiation, timing, heterogeneity between reefs); the wave and chl mechanisms need GBR-wide inputs |
| 2 | **Nested coarse boundary model**: a simpler reef-level COTS model over the source region only (e.g. Cooktown–Lizard and Cairns sectors), calibrated to regional outbreak timing, providing `P_j(y)` | Climate-responsive supply at lower cost; one-way nesting | Medium | The source region must contain the initiation area; calibration targets are coarser |
| 3 | **Use ReefMod-GBR (RME) projections as the boundary.** RME is on disk with GBR-wide CMIP6 DHW, cyclone simulations, COTS connectivity and a COTS control module (`rme_ml_2025_06_05`, MPL 2.0 code; check the data licences). Run 2026–45 under the same GCM / SSP, convert RME COTS density per source reef to `P_j(y)`, and calibrate the conversion so the RME hindcast matches `gbrPredsAdj` for 1991–2026 | An existing calibrated GBR-wide COTS model; climate-consistent; can include GBR-wide cull scenarios | Medium (RME runs, conversion, checks); ADRIA already has `RMEDomain` / `ReefModDomain` loaders | Structural mismatch between RME COTS and COTSMod; the conversion adds a layer of assumption |
| 4 | **Statistical boundary.** Either project the outbreak-probability SDM (`gbrPredsAdj`) forward if its covariates can be projected, or generate stochastic supply timelines by block-resampling the 1991–2026 `P_j(y)` series by regional wave (about 15-year cycles; seeded) | Fast; keeps observed wave periodicity, amplitude and spatial pattern | Small (days) | Not climate-responsive at the sources; resampling assumes the past wave regime continues; cannot represent culling outside the domain |
| 5 | **Bounding scenarios**: closed domain (no external supply), climatological mean `P`, repeat of the last wave | Brackets the importance of supply | Small | Not forecasts; for sensitivity and communication only |

**Recommendation.**
- **Now:** use option 4 (resampled waves) with option 5 as bounds to get a first projection quickly.
- **Next:** use option 3 (RME boundary) if the RME runs and licences are workable; it is the quickest climate-consistent boundary.
- **Long term:** keep option 1 (GBR-wide calibration) as the goal; it is the only route to endogenous, cull-responsive supply. It is a project in its own right, and the Lizard results (waves, pit, juvenile pipeline, external-layer timing) are its preregistered starting hypotheses.
- Option 2 is the fallback if option 1 proves too slow and option 3 too mismatched.

## 4. Suggested task order for the projection chat

1. Read the references above; reproduce run 20261019T_15 exactly. It is the replay control for everything that follows.
2. **Preregister** the projection design: scenarios, ensemble, metrics, the external-supply option, the 2026 start-state decision and the cull scenarios. Put the choices to the lead.
3. Build domain v0.5 (checked-in builder): CMIP6 DHW for 2027–45 mapped to sites, cyclone simulations, and a chl continuation. Before 2026 it must equal v0.4 exactly (replay test).
4. Implement the external-supply continuation (option 4 first) as a builder producing `P_j(y)` for 2027–45 per seed. Keep it opt-in; years up to 2026 unchanged.
5. Write `run_projection.jl`. A projection with every scenario set to "historical" must replay the calibration years exactly.
6. Forward culling adapter (effort → removals), opt-in and tested.
7. Run a small pilot (1 GCM, 2 SSPs, 3 supply scenarios, 5 seeds), audit it, log it, then scale up.
8. Separately, scope option 3 (RME boundary): a 1991–2026 RME hindcast against `gbrPredsAdj`.

### Progress (2026-10-21/22, rounds 30-33)

Details are in `WRITEUP_2026-10-22_projections_and_culling.md` and in the protocol from "Lead's decisions (2026-10-21, thirtieth round)".

- **Culling is on in projections.** Observed removals run to 2026, then forward effort culling: none, recent effort or targeted, each at ×1 to ×3.
  - Also added: a 2026 reset factor, coral-coupled outside supply, thermal-acclimation and outside-coral-larvae switches (all opt-in ADRIA switches).
- **Ensemble done:** 5 climate models × SSP2-4.5 / SSP3-7.0 × 5 seeds × reset × 9 cull levels (runs 20261021T_04–08).
- **Report figures:**
  - calibration figures with fit statistics;
  - culling projections to 2050;
  - 2020s build-up study;
  - outside-reef culling scenarios (`sandbox/data/derived/calibration_report_figures/2026-10-2*`).
- **Task 6** (forward culling adapter) and the **source-reef cull discount** (`build_outside_cull_discount.py`) are done.
- **Open:**
  - COTS-suitable habitat basis for culling (lead's proposal);
  - downward-only reset;
  - thermal-acclimation scenarios in the ensemble;
  - domain-level outputs;
  - RME boundary scoping (task 8).

### Progress (2026-10-19, round 29)

1. Done: the replay control is reproduced (run 20261019T_51: 0 = run 15, +0.04 = run 48, exactly).
2. Done: the design is preregistered (protocol, "Preregistered: projection design to 2100") and the lead's choices are recorded.
3. Done, as an **input set** rather than a copied domain: `build_projection_inputs.py`, with `test_projection_inputs.jl` 179/179.
4. Done: supply option 4 (wave replicates) and option 5 bounds (`mean`, `closed`, `repeat`) are in the input set.
5. Done inside the calibration runner (`OWEN_1991_PROJECTION=true`) rather than a separate `run_projection.jl`, so the W-wave-lf04 set-up is shared exactly. 1991–2026 replays exactly.
6. **Open:** the forward-cull adapter needs effort scenarios from the lead.
7. Done: pilot run 20261019T_58 is audited and logged. See the protocol for findings and the decisions needed before scale-up.
8. **Open:** RME boundary scoping. The pilot shows that stationary outside supply holds late-century coral at zero, which raises the priority of option 3 or a coral-coupled supply.

## 5. Data and working constraints

- **Internal or licensed data.** `SpatialData/GBR-HabitatMapping_INTERNALONLY`, `GBRMPA_RRN_2025_AIMS.xlsx`, the SALAD master workbook and the CSIRO Manta-Cull-RHIS workbook stay unpublished. Derived outputs stay in the git-ignored `sandbox/data/`; only summary statistics go in the protocol. The report figures (`sandbox/data/derived/calibration_report_figures/`) include SALAD/cull-derived densities (Figure 2), so check permissions before publication.
- **Commits.** The lead commits; do not commit unless asked.
- **No new coral parameters.** Do not add parameters to the ADRIA/CoralBlox coral codebase.
- **Runs.** Never overwrite runs; run `renumber_runs.py --apply` after creating re-score or gate directories.
- **G2′ checks** need `--trough-reefs "Lizard Island Reef" "MacGillivray Reef"` and the miniconda Python (`C:/Users/smatthew/AppData/Local/miniconda3/python.exe`).
- **Uncommitted work at handover (2026-10-19):** ADRIA (wave hook, cull allocation, external-survival hook, runner experiments, audits) and COTSMod (`starvation_bank`, parked). Ask the lead to commit before projection work starts.
