# COTSMod Calibration Workflow

The curated [model log](MODEL_LOG.md) tracks what has been tried, which runs are
valid evidence, and the current decision. The [crash-mechanism protocol](crash_mechanism_protocol.md)
freezes the next 1991 cohort-state factorial before interpretation; neither
document changes the supported model defaults. The generated
[run inventory](MODEL_RUN_INVENTORY.csv) indexes every non-hidden run directory,
including incomplete/unreviewed runs; its metadata statuses are not biological
promotion decisions.

This folder is the working seed for the standalone COTSMod calibration study repository. The goal is to calibrate COTS outbreak dynamics against Lizard Island reef observations while keeping the ecological model in `COTSMod.jl` and the ecosystem orchestration in `ADRIA.jl`.

> [!IMPORTANT]
> This file is the authoritative status and execution guide for the active COTS
> calibration workflow. Parameter values and commands in `AGENT_HANDOFF.md`,
> `sandbox/README.md`, and archived scripts describe earlier experiments unless
> they are explicitly restated here.
> Core-model changes are also governed by the repository-root `AGENTS.md`.

## Current Status

- COTSMod is integrated into ADRIA through the compatibility adapter in
  `ADRIA/src/ecosystem/cots.jl`.
- Package, adapter, and synthetic metric tests pass.
- The focused, expanded, pulse-enabled, and COTS-connectivity bounded pilots
  all failed the biological promotion criteria. None is a calibrated result.
- Every expanded candidate produced at most one detected peak per reef. The
  post-2005 temporal holdout matched zero of the four observed later peaks.
- A dedicated Lizard COTS matrix is built reproducibly from six ReefMod spawning
  seasons. The arithmetic mean remains the static control. Annual reef-scale
  matrices, source-specific larval survival, an external-GBR boundary flux, and
  juvenile storage have now been tested without changing the Allee threshold.
  None generated a second outbreak peak, so none is promoted.
- Peak timing/count and peak amplitude are the primary calibration targets.
  Pointwise error and correlation are supporting diagnostics.
- The calibration study must not be described as complete until it passes the
  pilot, held-out validation, and regional validation gates below.

### Amplitude convention

ADRIA adult COTS density is not currently known to be numerically equivalent to
AIMS COTS-per-tow observations. The hardened objective therefore fits one
non-negative observation-scale factor per reef on training years, then compares
peak heights and peak prominences in observed COTS-per-tow units. This preserves
relative outbreak amplitudes without silently treating incompatible units as
identical. The fitted scale is a nuisance observation-model parameter and must
be logged. A fixed empirical density-to-CPUE conversion should replace it when
one is available.

## Repository Boundary

The intended split is:

- `COTSMod.jl`: COTS state transitions, predation, larval dispersal, external supply hooks, and pure model tests.
- `ADRIA.jl`: domain loading, CoralBlox coupling, scenario sampling, result logging, and reef/site aggregation.
- calibration study repo: optimisation scripts, validation data, objective functions, plots, run metadata, and reproducible outputs.

During development these scripts still live under `sandbox/calibration`. Once stable, copy this folder, the plotting scripts, and shareable validation inputs into the calibration repo.

## Current Objective

The calibration target is not just pointwise agreement with observed COTS-per-tow values. We want simulations that reproduce outbreak cycles: presence of peaks, approximate timing, inter-peak spacing, amplitude, and non-flat dynamics.

The main reusable metric is implemented in:

```text
sandbox/calibration/cots_cycle_metrics.jl
```

The primary function is:

```julia
cots_cycle_score(sim_years, sim_values, obs_years, obs_values)
```

It returns a `CotsCycleScore` where `total_loss` is lower for better fits.

## Legacy Objective Components

The objective below describes the July-September 2026 implementation. It is
retained for traceability while the calendar-aware peak and amplitude objective
is being hardened. Do not use these weights for a large production sweep until
the metric validation gate is complete.

The cycle-aware loss combines:

- `matched_rmse`: RMSE at years where simulation and observations overlap.
- `percent_bias_abs`: absolute percent bias over matched years.
- `peak_count_penalty`: penalises missing or extra outbreak peaks.
- `peak_timing_penalty`: penalises simulated peaks that are too early or late.
- `period_penalty`: penalises inter-peak spacing far from the expected COTS cycle period, currently 15 years.
- `amplitude_penalty`: penalises trajectories with too much or too little relative variability.
- `flatline_penalty`: strongly penalises trajectories that do not produce meaningful cycles.
- `lag_correlation_penalty`: rewards shape similarity after allowing a bounded lag, while peak timing remains penalised separately.

The default total loss is:

```text
total_loss =
    0.35 * matched_rmse +
    0.002 * percent_bias_abs +
    1.00 * peak_count_penalty +
    1.50 * peak_timing_penalty +
    0.75 * period_penalty +
    0.75 * amplitude_penalty +
    2.00 * flatline_penalty +
    0.60 * lag_correlation_penalty
```

These weights are deliberately explicit so they can be reviewed and changed as calibration behaviour becomes clearer.

## Hardened Peak and Amplitude Objective

The active metric implementation now:

- collapses duplicate survey years and smooths by elapsed calendar years;
- restricts simulation and observations to their common time window;
- requires observations on both sides before declaring a peak;
- matches simulated and observed peaks one-to-one within a timing tolerance;
- scores peak height and prominence after fitting the documented per-reef
  observation scale;
- compares simulated and observed matched-peak periods rather than imposing a
  period when the observations contain fewer than two detected peaks; and
- gives peak timing and amplitude substantially more weight than generic RMSE,
  bias, or lagged correlation.

The legacy loss is retained in candidate logs as a diagnostic but is no longer
added to the optimisation target.

## Lagged Correlation

Supervisor feedback suggested using lags in correlation optimisation. The implementation now evaluates Pearson and Spearman correlation over a bounded lag window.

```julia
lagged_correlation(sim_years, sim_values, obs_years, obs_values; max_lag=6)
```

Interpretation:

- `best_lag_years == 0`: simulated cycle shape aligns without shifting.
- `best_lag_years > 0`: shifting the simulation later improves correlation, so the raw simulation is probably too early.
- `best_lag_years < 0`: shifting the simulation earlier improves correlation, so the raw simulation is probably too late.

Lagged correlation is used as a shape diagnostic and a modest objective component. It does not replace peak timing penalties. This is important because otherwise an optimiser could find a cycle with the right shape but wrong timing and still score too well.

## Cycle Metric Integration

The BlackBoxOptim driver (`calibrate_cots_blackbox.jl`) evaluates candidates using the cycle-aware objective function in `cots_cycle_metrics.jl`.

Each run directory records full cycle metrics in:
- `evaluated_candidates.csv`
- `best_summary.csv`
- `evaluated_by_reef.csv`

The reusable score object includes the following by-reef fields:
- `cycle_loss`
- `cycle_rmse`
- `cycle_abs_percent_bias`
- `cycle_peak_count_penalty`
- `cycle_peak_timing_penalty`
- `cycle_period_penalty`
- `cycle_amplitude_penalty`
- `cycle_flatline_penalty`
- `cycle_lag_correlation_penalty`
- `cycle_best_lag_pearson`
- `cycle_best_lag_spearman`
- `cycle_best_lag_years`
- detected simulated and observed peak years

The candidate CSV records reef-mean summaries, while `evaluated_by_reef.csv`
records each reef's scale, component penalties, matched-peak count, and detected
simulated and observed peak years. Successful and failed evaluations are both
logged.

The overall objective sorts candidates by:
```text
loss = mean_cycle_loss
```
The old pointwise `legacy_loss` is logged for diagnosis only and does not enter
the active optimization target.


## Validation Tests

Run the metric tests with:

```powershell
julia --project=sandbox sandbox\calibration\runtests.jl
```

The synthetic tests check that:

- correctly timed two-peak trajectories score best;
- shifted trajectories are recognised by lagged correlation but still penalised for timing;
- one-peak trajectories are penalised for peak count;
- flat trajectories are strongly penalised.

## BlackBoxOptim Calibration Driver

The BlackBoxOptim driver is implemented in:

```text
sandbox/calibration/calibrate_cots_blackbox.jl
```

### Search Modes & Parameter Sets

The driver supports three configurable search modes via the `BBO_MODE` environment variable:

1. **`FOCUSED` (Default)**: 4 core parameters for fast tuning and smoke testing:
   `a_F`, `a_S`, `IMM`, `seed_mult`.
2. **`EXPANDED`**: 14 demographic, functional response, mortality, and dispersal parameters:
   `a_F`, `a_S`, `IMM`, `seed_mult`, `a_ricker`, `b_ricker`, `m1`, `m2`, `m3`, `p_tilde`, `C_max`, `tau_condition`, `imm_threshold`, `eta_imm`. The Allee threshold is fixed separately and is not optimized.
3. **`EXPANDED_PULSE`**: Includes all expanded parameters plus pulse timing and magnitude:
   `pulse_start`, `pulse_duration`, `pulse_relative_magnitude`.

### Environment Controls

```powershell
$env:BBO_MODE = 'FOCUSED'     # Options: FOCUSED, EXPANDED, EXPANDED_PULSE
$env:BBO_MAX_STEPS = '20'     # Maximum evaluation steps
$env:BBO_MAX_EVALS = '0'      # Function-evaluation request (0 = package default)
$env:BBO_MAX_TIME = '0.0'     # Time limit in seconds (0.0 = unlimited)
$env:BBO_METHOD = 'adaptive_de_rand_1_bin_radiuslimited'
$env:BBO_SEED = '20260928'    # Optimizer and scenario-template seed
$env:BBO_RUN_ID = 'pilot_001' # Reusable run identifier
$env:BBO_RESUME = 'false'     # Reuse logged candidate evaluations when true
$env:ADRIA_COTS_CONNECTIVITY_MODE = 'cots' # auto|cots|coral
$env:COTS_ALLEE_THRESHOLD = '3.0' # COTS/ha; fixed mechanism contract
$env:COTS_CONNECTIVITY_TEMPORAL_MODE = 'mean' # mean|cycle|sample
$env:COTS_CONNECTIVITY_SEED = '20260929' # used by sample mode
$env:COTS_APPLY_LARVAL_SURVIVAL = 'false'
$env:COTS_EXTERNAL_SOURCE_DENSITY = '0.0' # legacy name: external pre-dispersal recruits ha^-1 yr^-1
$env:COTS_JUVENILE_STORAGE = 'false'
$env:COTS_HABITAT_MEDIATION = 'false'
$env:COTS_SEPARATED_RECRUITMENT = 'false'
$env:COTS_LARVAL_FECUNDITY = '6.0' # pre-pelagic COTS ha^-1 yr^-1 scale
$env:COTS_SETTLEMENT_PROBABILITY = '1.0'
$env:COTS_EXTERNAL_OUTBREAK_PRODUCTION = '0.0' # probability-1 source production
$env:COTS_CALENDAR_START_YEAR = '1985' # required when domain time is relative
$env:BBO_EXCLUDE_REEFS = ''    # Semicolon-delimited Lizard holdouts
julia --project=sandbox sandbox\calibration\calibrate_cots_blackbox.jl
```

Each run writes to `sandbox/calibration/runs/<run-id>/`. The directory contains
run metadata, successful and failed candidate rows, by-reef peak diagnostics,
the sorted summary, simulations, and figures. Run metadata records git
revisions, input hashes, search bounds, seeds, environment controls, and Julia
version. `BBO_RESUME=true` restarts the optimiser but reuses exact candidates
already present in the evaluation cache; BlackBoxOptim state itself is not
checkpointed.

### Candidate Evaluation Logging

Evaluations are saved under the run directory:

```text
sandbox/calibration/runs/<run-id>/evaluated_candidates.csv
sandbox/calibration/runs/<run-id>/evaluated_by_reef.csv
sandbox/calibration/runs/<run-id>/best_summary.csv
```

Each successful row contains the candidate vector and reef-mean metric
breakdown. Failed rows include status and error text. Exact previously evaluated
candidate vectors are cached when `BBO_RESUME=true`; optimizer population state
is not checkpointed.

## Unified Calibration & Visualization Pipeline

The calibration workflow consists of three clean, non-duplicated steps:

```text
1. Optimization Sweep:   calibrate_cots_blackbox.jl
                           │ Outputs: bbo_best_summary.csv
                           ▼
2. Unified Simulation:  simulate_best_calibration.jl
                           │ Outputs: best_calibrated_trajectories.csv
                           ▼
3. Figure Rendering:     render_calibration_figures.py
                           │ Outputs: plots/*.png
```

### Execution Example

```powershell
# Step 1: Run BlackBoxOptim driver (Focused or Expanded mode)
$env:BBO_MODE = 'FOCUSED'
$env:BBO_MAX_STEPS = '50'
julia --project=sandbox sandbox\calibration\calibrate_cots_blackbox.jl

# Step 2: Run simulation on top candidate (Single-run or Stochastic Ensemble)
$env:BBO_RUN_ID = 'pilot_001'
$env:COTS_N_STOCHASTIC_SCENS = '1'
$env:COTS_STOCHASTIC_MODE = 'demographic' # environmental|demographic|combined
julia --project=sandbox sandbox\calibration\simulate_best_calibration.jl

# Step 3: Render publication-quality plots
python sandbox\calibration\render_calibration_figures.py
```

## Hardening and Development Gates

The active development sequence is:

1. **Authoritative status:** keep this README current and mark older calibration
   instructions as historical.
2. **Metric validation:** use calendar-aware observations, a common scoring
   window, explicit peak matching, and observation-scaled peak height and
   prominence errors. Synthetic tests must cover irregular surveys, boundary
   years, missing peaks, incorrect second-peak amplitude, and flat trajectories.
3. **Reproducible runs:** record configuration, seeds, git revisions, package
   versions, data fingerprints, environment controls, failures, and by-reef
   scores in a unique run directory. Runs must be resumable.
4. **Pipeline modes:** smoke-test focused, expanded, pulse, deterministic, and
   stochastic modes. Ensemble figures must show median and uncertainty bands.
5. **Pilot calibration:** run a bounded search with replicated candidates and
   convergence checks. Freeze the objective before a production-scale sweep.
6. **Sensitivity and mechanisms:** apply Sobol/PAWN analysis around the validated
   region and compare coral-proxy versus COTS-specific connectivity and
   pulse/no-pulse mechanisms.
7. **Validation and extraction:** hold out reefs or years, then validate on
   Moore/Cairns domains. Extract the study repository and publish/version
   COTSMod only after these gates pass.

Gates 1-5 now have working implementations and bounded results. Gate 6 has a
validated COTS-connectivity build and paired six-candidate comparison, but a
formal Sobol/PAWN analysis is intentionally deferred because no plausible
two-peak region exists. Gate 7's temporal holdout was executed and failed on all
four reefs. Regional validation remains pending because the local Moore package
covers 2025-2099 rather than the 1985-2024 observation period.

The next model-development gate is therefore mechanistic: identify why the
current internal dynamics decay after the first outbreak and cannot regenerate
a second peak. The exploratory screen prioritizes `a_ricker`, followed by
`seed_mult` and `IMM`, for targeted recurrence experiments, but its six-point
rank correlations are not a formal sensitivity result. Do not launch a
production or multi-seed convergence sweep until at least one candidate produces
the observed second peak.

### Recurrence diagnosis

The best COTS-connectivity trajectory now exports recruits, juveniles, adults,
body condition, and coral cover. By 2024 the four calibration reefs have
recovered body condition of 0.91-0.96 and coral cover of 0.41-0.46, but adult
density remains only 0.094-0.114 and recruits only 0.044-0.047. The candidate's
Allee threshold is 3.20, so the final adult-to-threshold ratio is 0.029-0.036
and the fecundity multiplier
`N_adult^2 / (allee_threshold^2 + N_adult^2)` is only 0.00085-0.00127.

This identifies a low-density/Allee trap rather than persistent food limitation.
It does **not** establish that the independently supported Allee threshold is
wrong. The calibration contract now fixes `COTS_ALLEE_THRESHOLD=3.0` COTS/ha
during recurrence experiments. Threshold variation is retained only as a
labelled counterfactual.

### Seven-step recurrence mechanism experiment

The next bounded experiment follows these ordered steps:

1. **Units and fixed threshold.** Treat COTS state, immigration, and the Allee
   threshold as COTS/ha. Hold the Allee threshold at 3.0 COTS/ha and keep the
   CPUE conversion in the observation model.
2. **Open boundary.** Partition recruitment into local production, internal
   Lizard transport, and an absolute external-GBR boundary flux. A sink may
   receive settlers while its adults remain below the spawning threshold.
3. **Annual connectivity.** Use the six individual ReefMod COTS matrices in a
   recorded deterministic cycle or recorded seeded draw. The arithmetic mean is
   the control, not the treatment.
4. **Source survival.** Multiply transported larvae by the reef- and year-
   specific COTS larval-survival data before dispersal. Record the water-quality
   scenario and hydrodynamic year.
5. **Juvenile storage.** Test an opt-in retained-juvenile state in which
   maturation increases with coral food and settlement responds to an open-
   habitat proxy. The legacy immediate-maturation path remains the default until
   promotion. A true rubble/recent-disturbance state is deferred until an ADRIA
   rubble input is available and validated.
6. **Bounded factorial.** Compare static/closed, annual/internal, annual plus
   larval survival, annual plus survival plus the external boundary, and the
   same treatment plus juvenile storage. Test habitat mediation separately so
   it is not confounded with storage. Use identical demographic parameters and
   seeds.
7. **Mechanism gate.** Log local fecundity, background immigration, internal
   immigration, external supply, settlement habitat, juvenile retention, and
   maturation. Reject a treatment unless it produces a biologically plausible
   second peak in timing and amplitude before optimizer expansion.

Do not copy CoCoNet's explicit approximately 16-year spawning oscillator into
the primary experiment. It is an exogenous synchronizer and would make it hard
to determine whether connectivity and cohort dynamics generate recurrence.

### Recurrence-factorial result

The corrected six-treatment run
`recurrence_mechanism_pilot_v3_seed20260929` uses the best prior expanded
candidate, seed `20260929`, a fixed Allee threshold of 3.0 COTS/ha, cyclic annual
connectivity, and q3baseline source survival. The table reports four-reef means;
detected and matched peaks are totals across reefs.

| Treatment | Mean loss | Detected peaks | Matched peaks | Timing penalty | Amplitude penalty |
|---|---:|---:|---:|---:|---:|
| Static mean, closed | 4.941 | 4 | 3 | 0.792 | 0.892 |
| Annual internal | 4.945 | 4 | 3 | 0.792 | 0.893 |
| Annual + source survival | 5.175 | 3 | 3 | 0.750 | 0.897 |
| Annual + survival + external | 5.175 | 3 | 3 | 0.750 | 0.897 |
| Same + juvenile storage | 6.425 | 1 | 1 | 0.896 | 0.965 |
| Same + storage + habitat proxy | 7.125 | 0 | 0 | 1.000 | 1.000 |

No treatment generated a second simulated peak. The static control peaks in
1993-95. After fitting the documented reef-specific observation scale, its peak
heights are 0.081, 0.161, 0.570, and 0.311 COTS/tow for Eyrie, North Direction,
Lizard, and MacGillivray respectively. Observed peaks are 0.85; 0.56 and 0.71;
1.92 and 0.74; and 0.86 and 0.67 COTS/tow. The main failure is therefore both
missing recurrence and insufficient amplitude.

Annual connectivity alone is effectively neutral. Applying q3baseline source
survival reduces mean internal immigration from 0.216 to 0.000234 COTS/ha per
year because the six annual mean survival factors range from 0.000106 to
0.003538. This is the supplied ReefMod factor, but the result shows that its
relationship to COTSMod's fitted Ricker recruitment scale must be audited before
the factor enters calibration.

At the screening value of 3 pre-dispersal recruits/ha/year, mean external
immigration is `1.55e-5` COTS/ha/year: about 6.6% of survival-weighted internal
immigration and 0.1% of background immigration. It is too small to test
connectivity-mediated escape from the Allee trap. The environment variable name
`COTS_EXTERNAL_SOURCE_DENSITY` is retained for compatibility, but its quantity
is pre-dispersal recruit production, not external adult COTS density. A future
boundary series must be derived from observed or simulated upstream outbreak
production; it must not be tuned as an undocumented pulse.

Storage alone lowers mean maturation from 0.166 to 0.112 COTS/ha/year and leaves
only the MacGillivray peak. Adding the open-habitat proxy lowers the mean
settlement gate to 0.340 and removes all detected peaks. The proxy is rejected
until a true rubble/recent-disturbance state is available.

The decision is **no promotion and no production optimization**. This result
motivated the separated production model and evidence-derived boundary below.
The Allee threshold remains fixed at 3.0 COTS/ha.

### Separated production and evidence-derived boundary

The experimental path now separates three stages that were previously combined
inside `a_ricker`:

1. `larval_fecundity` generates potential larvae from adult density, body
   condition, density dependence, and the fixed Allee effect;
2. reef- and hydrodynamic-year-specific q3baseline survival acts on production
   at the source; and
3. `settlement_probability` and the optional habitat gate act once at the
   sink.

The legacy effective-recruitment equation remains the default. The separated
path is enabled only with `COTS_SEPARATED_RECRUITMENT=true`. Local retention,
internal immigration, external potential production, external pelagic
survivors, and settled recruits are logged separately.

In the experimental separated path, potential local production is
`larval_fecundity × body_condition² × D × exp(-b_ricker × D) ×
D²/(D²+3²)`, where `D` is adult COTS/ha. At low density it is approximately
proportional to `D³`, so the existing fertilisation proxy is strongly
density-dependent. This is a phenomenological Allee term, not an explicit
simulation of sex ratios or fertilisation encounters. The prescribed upstream
boundary does **not** automatically inherit this local adult-density term.

The upstream boundary is derived from
`sandbox/data/gbrPredsAdj_20262408.csv`. It contains 35 annual predictions
(1991-2025) for the probability that a reef exceeds the outbreak definition of
0.22 COTS/tow. The probability is used continuously as expected source
activity; it is not thresholded at 0.22 or converted to a binary outbreak.
Predictions map to 3,805 of 3,806 ReefMod reefs after the six aliases documented
in the ReefMod ID file. Reef `99-407` is the sole missing source and is
assigned zero probability. Fifty-six prediction-only reef IDs are ignored.

For sink reef `j`, calendar year `t`, and hydrodynamic season `h`,
the potential boundary multiplier is:

```text
B_potential[j,t,h] =
    sum(i outside Lizard,
        outbreak_probability[i,t] *
        connectivity[i,j,h] *
        habitable_area[i] / habitable_area[j])
```

The pelagic multiplier additionally includes source survival. Multiplying these
dimensionless coefficients by `COTS_EXTERNAL_OUTBREAK_PRODUCTION` produces
potential and surviving sink-density fluxes in COTS/ha/year. Settlement is then
applied once in COTSMod.

### Separated-production pilot result

The bounded run
`separated_production_evidence_pilot_v4_seed20260929` used a
scale-preserving initial decomposition:

- mean Lizard q3baseline survival: 0.00108268;
- pre-pelagic fecundity: `6.72559 / 0.00108268 = 6211.96`;
- probability-1 production at the 3.0 COTS/ha Allee density: 7810.96
  COTS/ha/year; and
- probability-1 production at the Ricker maximum (17.93 adults/ha):
  37749.74 COTS/ha/year.

Settlement was held at its upper bound of 1.0, with juvenile storage and habitat
mediation disabled.

| Treatment | Mean loss | Detected peaks | Matched peaks |
|---|---:|---:|---:|
| Legacy annual internal | 4.945 | 4 | 3 |
| Separated, hydrodynamic mean, closed | 7.124 | 0 | 0 |
| Mean + evidence boundary at Allee-source production | 7.121 | 0 | 0 |
| Mean + evidence boundary at Ricker-peak production | 6.639 | 2 | 1 |
| Annual cycle + evidence boundary at Ricker-peak production | 5.535 | 4 | 3 |

The high-production evidence boundary creates a real recurrence response, but
only at Eyrie. Under hydrodynamic-mean forcing, Eyrie peaks in 2001 and 2017 at
0.107 and 0.076 COTS/tow after observation scaling; its observed peak is 0.85
COTS/tow in 2013. Cycling the six hydrodynamic matrices shifts the second Eyrie
peak to 2019 and amplifies supply whenever the strong 2010 matrix repeats.
Lizard, North Direction, and MacGillivray still do not produce their observed
second outbreaks.

Mean surviving external supply under the high-production/mean treatment is
0.218 COTS/ha/year at Eyrie, 0.042 at Lizard, 0.039 at North Direction, and
0.012 at MacGillivray. This spatial imbalance means that increasing one global
production scalar would overdrive Eyrie before repairing the other reefs.

The separated formulation is retained as an experimental mechanism because it
is unit-consistent and can generate recurrence. It is not promoted: timing,
amplitude, and cross-reef behavior fail the peak contract. The next gate is to
validate source/sink area scaling and the temporal relationship between
upstream predicted adult outbreaks, spawning, settlement, and observed adult
peaks. Juvenile maturation remains deferred.

### Lizard V2 domain migration

The final Far North GBR subreef polygons and interim distance-based site matrix
have been integrated as the additive `Lizard_Historical_v0.2` package. The
calibration control remains `Lizard_Historical_v0.1`; existing losses and peak
metrics must not be compared directly with v0.2 because v0.2 changes the site
tessellation, polygon areas, COTS raster sampling, and historical DHW.

V2 uses `indexV2` as the canonical order and stable site identifiers
`FNG_V2_0001:FNG_V2_2895`. The supplied `site_id` is retained only as
`source_site_id` because it is not unique. Six invalid polygons are repaired,
area is recomputed in EPSG:32755 in square metres, and the unlabelled 2,895 by
2,895 workbook is labelled in strict `indexV2` order after finite,
non-negative, dimension, and symmetry checks. Complete source and output hashes
are recorded in `sandbox/data/Lizard_Historical_v0.2/provenance.json`.

The six ReefMod COTS matrices remain the default `legacy` dataset on v0.2. The
five Owen annual matrices (2018–2022) are available only through
`ADRIA_COTS_CONNECTIVITY_DATASET=owen_global_2018_2023`. Their current artifact
is experimental: the first four matrices have wholly non-finite source rows, 56
new sources lack ReefMod survival/area data, and the supplied survival series
ends in 2017. The provisional build zeros whole non-finite rows, excludes the 56
unmapped external sources, and applies the latest available 2017 survival
modifier. Every policy is explicit in provenance and none is promoted.

Simulation and diagnostic scripts accept `LIZARD_DOMAIN_VERSION` while retaining
v0.1 as the default. Run metadata records the domain and connectivity dataset.
The BlackBoxOptim driver explicitly rejects v0.2: its current unweighted reef
mean is not the tested area-weighted comparison mapping, and no V2 treatment
has passed the qualitative peak gate. The area-weighted mapping is now frozen
for bounded comparisons in `LIZARD_OBSERVATION_CONTRACT.md`, not promoted into
the optimizer. The Allee threshold remains fixed; survival,
missing-source, and area contracts must be resolved before juvenile maturation
or settlement is revisited. Full build and validation commands are in
`sandbox/domain_building/README.md`.

A bounded, same-seed V2 pilot compared legacy versus Owen connectivity, each
with a closed or evidence-driven upstream boundary. It used the same frozen
biological candidate, seed `20260930`, an Allee threshold of 3 COTS/ha,
polygon-area-weighted reef means, and observation scales estimated once from
the legacy/closed control. The boundary input was `38827.55` potential
recruits/ha/year, derived as a provisional production ceiling and multiplied
by the source outbreak probabilities and transport; it is **not** measured
recruitment. Exact per-reef scores, CPUE peak heights, trajectories, forcing
years, and hashes are in
`sandbox/calibration/runs/20261001T103104_lizard_v2_connectivity_seed20260930`.

| Treatment | Mean loss | Detected peaks across four reefs | Matched peaks |
| --- | ---: | ---: | ---: |
| Legacy, closed | 5.799 | 2 | 2 |
| Legacy, evidence boundary | 5.515 | 4 | 3 |
| Owen, closed | 7.148 | 0 | 0 |
| Owen, evidence boundary | 4.511 | 5 | 4 |

The Owen/boundary treatment produced one matched peak per reef, but its main
2017–18 peaks were late and low: Eyrie `0.0232`, North Direction `0.0771`,
Lizard `0.1169`, and MacGillivray `0.0810` COTS/tow after the *frozen* observation
conversion. North Direction also had a 2004 peak (`0.0659` COTS/tow). None of
the four 2017–18 heights reaches the `0.22` COTS/tow outbreak threshold. Owen
alone removed the peaks in this candidate. Thus the external boundary can
reseed the open subset, but neither this boundary magnitude nor the provisional
Owen forcing passes the peak timing/amplitude gate. Do not promote it or tune
juvenile maturation against it. This exploratory run records revisions and
input/output hashes, but not package versions or a dirty-worktree diff; those
are required before a promotion-grade validation. Reproduce with:

```powershell
julia --project=sandbox sandbox\calibration\compare_lizard_v2_connectivity.jl
```

#### Paired domain gate, Owen coverage, and observation smoothing

The checked-in `reef_observation_mapping.jl` now computes the within-reef
polygon-area-weighted COTS density from a parent-reef-ID crosswalk, with unit
and alignment assertions. The V1-control CPUE conversion is fitted once on
raw observations at seed `20260930` and held fixed across V1/V2, connectivity,
boundary, seed, and smoothing treatments. The named three-year observation
sensitivity averages available CPUE in the survey year and adjacent calendar
years without filling gaps. Raw scores are retained. The five-year alternative
was not adopted: it shifted peaks further and introduced a low Eyrie peak.

The Owen coverage audit is
`20261001T115318_owen_coverage_audit`. Of 56 unmatched external source reefs,
only ten source-year rows have nonzero raw transport into the 113 Lizard reefs
under the provisional zero-nonfinite-row policy. The excluded share of *raw*
external connectivity averages 0.14% across sink-years, but reaches 10.6% at
one sink in 2022. This is not a recruitment fraction: the missing sources still
lack habitable area and survival. For mapped sources, replacing reused 2017
survival with each source's 2010–17 mean changes the external coefficient by a
median factor of 1.84 across sink-years. Neither historical mean nor extrema
are evidence for actual 2018–22 survival. These assumptions remain unresolved.

The paired V1/V2 `cycle` control is
`20261001T112847_lizard_domain_gate`. It was repeated with three seeds, but
all results were identical because `cycle` forcing is deterministic. The
separate sampled-hydrodynamic-year sensitivity
`20261001T113401_lizard_domain_gate` used three recorded seeds and produced
variation. Each run compared V1 and V2 legacy connectivity, plus V2 Owen,
with closed/evidence-boundary treatments and unchanged biological parameters.
V1/V2 is a *whole-domain* comparison: the target reefs' site counts and areas
also change substantially (for example Lizard Island has 65 sites/27.77 km²
in V1 versus 63 sites/6.13 km² in V2). It does not isolate connectivity.

| Sampled-year treatment | Mean loss, raw | Mean loss, 3-year observed | Matched / 21 observed peaks |
| --- | ---: | ---: | ---: |
| V1 legacy, closed | 5.970 | 5.779 | 6 |
| V2 legacy, closed | 5.963 | 5.770 | 6 |
| V1 legacy, boundary | 5.279 | 5.083 | 10 |
| V2 legacy, boundary | 4.826 | 4.648 | 12 |
| V2 Owen, closed | 7.144 | 7.122 | 0 |
| V2 Owen, boundary | 4.398 | 4.489 | 12 |

Owen/boundary again matches one peak per reef and seed, but its late peaks
(mostly 2017–18) range only `0.026–0.190` COTS/tow across reefs/seeds; none
reaches the `0.22` COTS/tow outbreak definition. Smoothing improves the legacy
mean losses slightly but does not repair this timing/amplitude failure (and
slightly worsens Owen/boundary mean loss). The [updated peak plot](runs/20261001T113401_lizard_domain_gate/plots/cots_peak_comparison_gapaware.png)
shows raw survey points, smoothed observations, model medians and seed bands,
detected peaks, and the outbreak threshold. The fixed-scale and smoothing
contract is in `LIZARD_OBSERVATION_CONTRACT.md`; all per-reef scores and
trajectory outputs are in the two run directories. `sample` is a forcing-year
stress test, not a historical hindcast, and no treatment is promoted.

```powershell
python sandbox\domain_building\audit_owen_forcing_gaps.py
julia --project=sandbox sandbox\calibration\test_reef_observation_mapping.jl
julia --project=sandbox sandbox\calibration\test_cots_cycle_metrics.jl
$env:LIZARD_GATE_TEMPORAL_MODE = 'cycle'
$env:LIZARD_GATE_RUN_ID = 'choose_an_unused_cycle_run_id'
julia --project=sandbox sandbox\calibration\compare_lizard_domain_gate.jl
$env:LIZARD_GATE_TEMPORAL_MODE = 'sample'
$env:LIZARD_GATE_RUN_ID = 'choose_an_unused_sample_run_id'
julia --project=sandbox sandbox\calibration\compare_lizard_domain_gate.jl
python sandbox\calibration\plot_lizard_domain_gate.py --run-id choose_an_unused_sample_run_id
```

#### Coral-cover diagnostic and old-wave replay (2026-10-01)

**Observation correction:** the first coral plot used an exact model-reef-name
lookup and therefore found only Eyrie. That was wrong: the LTMP coral files use
the same aliases already frozen for the COTS tow mapping (`Lizard Isles`,
`Macgillivray Reef`, `North Direction Island`, and `Eyrie Reef`). The corrected
[V1/V2 COTS and coral plot](runs/20261001T113401_lizard_domain_gate/plots/cots_coral_gate_ltmp_manta_photo.png)
and [old-wave replay plot](runs/20261001T132602_lizard_peak_mechanism_replay/plots/cots_coral_replay_ltmp_manta_photo.png)
use `reef_manta.csv` for LTMP `HC`/`MANTA` at all four reefs and
`reef_photo_transect.csv` for reef-level `HARD CORAL`/`GROUP_LEVEL` at Lizard,
MacGillivray, and North Direction. Both are 9 m methods; Eyrie has no matching
photo-transect series in the supplied file. Photo composition rows are **not**
added to the group-level hard-coral estimate. Within the 1985-2024 model
window, the four manta series contain 32, 33, 27, and 11 observations; the
three photo series contain 23, 24, and 23. The checked-in
`plot_lizard_ltmp_coral.py` writes the exact mapped rows, source hashes, and
method filters beside each plot. The earlier Eyrie-only plots remain archived
as superseded diagnostics, not as the current evidence.

The corrected panels also show the archived old model cover as an explicitly
**unweighted** site mean. LTMP observations are reef-level 9 m estimates,
whereas the current model line is polygon-area-weighted whole-reef cover.
Source-reported intervals are shown; neither method is silently converted to
a whole-reef observation or inserted into the frozen COTS objective.

The closed runs have very low adult COTS while coral cover rises to roughly
0.35-0.50 in later years. In this comparison, low modeled coral food is not
the immediate cause of the missing COTS peak. The provisional Owen/boundary
treatment increases adults and depresses coral, but still misses peak height.
The corrected observations show a substantial **multi-reef** coral-response
discrepancy in 2015, not only an Eyrie discrepancy:

| Reef | LTMP manta HC | LTMP photo hard coral | V2 Owen/boundary model, seed 20260930 |
| --- | ---: | ---: | ---: |
| Lizard Island | 0.079 | 0.085 | 0.280 |
| MacGillivray | 0.084 | 0.118 | 0.283 |
| North Direction | 0.068 | 0.088 | 0.235 |
| Eyrie | 0.105 | unavailable | 0.308 |

These are diagnostic comparisons across different spatial measurement
supports, not direct coral-cover calibration residuals. Coral cover is not
consistently high in the model throughout the record: the corrected plots
also show early and mid-period method/model differences that vary by reef.
The low observed 2015 cover does not by itself prove whether missing COTS
predation, coral disturbance/recovery, or spatial observation support causes
the mismatch. It does rule out the previous interpretation that Eyrie was the
only reef with usable coral observations.

The one-seed, static-transport V1 replay
`20261001T132602_lizard_peak_mechanism_replay` uses a frozen per-reef CPUE
scale and polygon-area-weighted site states. Its
[paired COTS/coral plot](runs/20261001T132602_lizard_peak_mechanism_replay/plots/cots_coral_peak_replay.png)
compares the old parameter candidate with coral-proxy and COTS connectivity
at the evidence-based `3 COTS/ha` Allee threshold. The `1 COTS/ha` replay is
**only a labelled historical counterfactual**, not a proposed threshold.

| V1 replay treatment | Mean 3-year-observation loss | Matched / 7 observed peaks | 2015 Eyrie model coral cover |
| --- | ---: | ---: | ---: |
| Old candidate, coral proxy, Allee 1 (counterfactual) | 2.785 | 7 | 0.053 |
| Old candidate, coral proxy, Allee 3 | 5.456 | 2 | 0.260 |
| Old candidate, COTS matrix, Allee 3 | 4.891 | 3 | 0.260 |
| Current candidate, COTS matrix, legacy recruitment, Allee 3 | 5.311 | 2 | 0.384 |
| Current candidate, separated production, no source survival, Allee 3 | 15.701 | 3 | 0.000 |
| Current candidate, separated production plus source survival, Allee 3 | 6.352 | 1 | 0.446 |

The counterfactual nearly recreates the archived *model* coral trajectory on
several reefs and supplies second peaks on three of four reefs. It does not
validate coral dynamics against LTMP: in 2015 its area-weighted cover is
0.084, 0.175, 0.141, and 0.053 across the four reefs, compared with manta
medians 0.079, 0.084, 0.068, and 0.105. Its late COTS peaks also remain
several years later than observed. Raising only the Allee threshold to 3 removes
most recurrence; changing the network alone does not restore it. Thus the old
figure is not evidence that the V2 site migration caused the missing waves.
The archived environment was not fully recorded, so this is a replay of its
parameter path, not a certified exact reproduction.

The replay's 2008-18 mean fluxes expose a second bottleneck: with separated
production and source survival, mean local fecundity across the four reefs is
`0.783` pre-pelagic recruits/ha/year, but pelagic survivors and settled
recruits are only about `0.0005`/ha/year. Background immigration is
`0.0303`/ha/year and dominates recruitment; maturation is `0.0126`/ha/year.
Switching source survival off creates gross overproduction and zero modeled
coral cover by 2015; it is a diagnostic control, not a rescue setting. This
large bracketing response makes the provenance and units of source survival,
transport, and settlement the next priority before tuning fecundity or
juvenile parameters.

```powershell
$env:LIZARD_REPLAY_RUN_ID = 'choose_an_unused_replay_run_id'
julia --project=sandbox sandbox\calibration\replay_lizard_peak_mechanisms.jl
python sandbox\calibration\plot_lizard_ltmp_coral.py --run-id choose_an_unused_replay_run_id --kind replay
# On a fresh paired gate run, use --kind gate with its unused run ID.
```

A subsequent one-seed V2 Owen/boundary stage-survival screen is
`20261001T133241_lizard_stage_survival_screen`. It retains `3 COTS/ha`, the
same sampled hydrodynamic-year seed, separated production, source survival,
and the evidence-derived boundary within the screen. It restores the older
*default* juvenile (`m2=0.2`) and/or adult (`m3=0.1`) mortality only as
labelled sensitivities; these defaults are not independently validated bounds.

| Stage-survival treatment | Mean 3-year-observation loss | Matched / 7 peaks | 2010-20 Lizard maximum COTS/tow | 2015 Eyrie coral cover |
| --- | ---: | ---: | ---: | ---: |
| Current `m2=0.317`, `m3=0.195` | 4.643 | 4 | 0.138 | 0.336 |
| Legacy `m2` only | 4.618 | 4 | 0.161 | 0.320 |
| Legacy `m3` only | 4.865 | 3 | 0.222 | 0.266 |
| Legacy `m2` and `m3` | 4.864 | 3 | 0.254 | 0.247 |

Lower adult mortality raises the late Lizard trajectory past the `0.22`
outbreak line, but does not rescue Eyrie or MacGillivray peaks and worsens the
frozen aggregate score. With both older mortality defaults, 2015 modeled
cover is 0.198, 0.206, 0.152, and 0.247 across the four reefs: still above
all four 9 m manta medians, though these are not like-for-like spatial
estimates. No treatment passes the qualitative multi-reef peak/coral gate.
This screen recalculates pre-pelagic fecundity from the **Owen** source-survival
mean (`4532` recruits/adult/year), whereas the earlier paired V1/V2 gate used
the V2 legacy mean (`6389`); its control should therefore be compared with
the other treatments *within this screen*, not treated as an exact repeat of
the earlier Owen/boundary control. The external-source omissions and reused
survival remain provisional. Do not promote mortality settings from this
single-seed exploratory check.

```powershell
$env:LIZARD_SURVIVAL_RUN_ID = 'choose_an_unused_survival_run_id'
julia --project=sandbox sandbox\calibration\screen_lizard_stage_survival.jl
```

#### Owen survival/flux and LTMP coral-support gate (2026-10-01)

The checked-in `audit_owen_flux_contract.jl` replays the frozen V2 Owen plus
evidence-boundary treatment at seed `20260930`, with separated production,
`3 COTS/ha` Allee threshold, source survival, unit settlement probability, and
the gate's seed-controlled sampled hydrodynamic years. Its source audit
confirms that **all 113 local reef survival factors** in each 2018–22 forcing
year equal that reef's RME
`q3baseline` **2017** value (mean `0.001484`; the 2010–17 per-reef mean,
averaged over reefs, is `0.001121`). This establishes what the adapter uses,
not what survival actually was in 2018–22. The RME connectivity notes define
the 0–1 factor as a multiplier on larvae produced at the source. The supplied
Owen package contains five matrix CSVs and reef geometries but no method note
that establishes whether those transport probabilities already include the
same mortality. That double-counting question, the 56 unmatched external
sources, and the post-2017 survival gap remain unresolved; the treatment is
still experimental.

For the implemented separated path, the audit reconstructs each year's
source-to-sink chain from the logged site production and selected matrix:
source production is the **equal-site mean** within each Owen reef, multiplied
once by source survival; transport then multiplies by `source reef area / sink
reef area` and the source-row/sink-column matrix; site settlement applies the
site gate once. The external boundary is separately generated from the
time-varying outbreak-probability coefficient, with source survival already
inside its pelagic coefficient. Background immigration remains a separate
flux. Lizard's survey reef grouping covers **three** Owen nodes (`14-116b`,
`14-116c`, `14-116d`), so the audit reconstructs at site level before applying
the frozen polygon-area observation weights. In the passing
[`20261001T170000_owen_flux_audit`](runs/20261001T170000_owen_flux_audit/metadata.toml),
all reconstructed internal, external, settlement, and recruit-stage balances
agree with the logs to below `1e-7 COTS/ha` (maximum absolute difference about
`3.4e-13`). All 156 target reef-years match the archived gate's selected
hydrodynamic years and adult/coral trajectories to floating-point precision.
This verifies code arithmetic and timing, **not** biological
validity of the source factor or boundary units.

In 2008–18 the four survey-reef means are `18.8–129.8` pre-pelagic production,
`0.0020–0.0085` local-network pelagic supply, and `0.153–0.309` external
pelagic supply, all in `COTS/ha/year` equivalents. This explains why the
provisional boundary drives the outcome more strongly than local reef
transport. The current equal-site source mean differs from a polygon-area
mean; for North Direction the area-weighted/equal-site production ratio has a
2008–18 median of `1.15` (other three survey groupings are near `1.00`). This
is a **new candidate mechanism/implementation hypothesis**, not an authorized
silent correction: a source-area-weighted treatment needs a mass-balance test
and same-seed control comparison before promotion.

The separate [`LTMP coral-support audit`](runs/20261001T172000_ltmp_coral_mapping_audit/metadata.json)
keeps raw manta and photo-transect medians distinct and writes paired
reef/report-year records with source IDs and survey dates. Paired surveys
have the same recorded date. Across Lizard, MacGillivray, and North Direction
there are `23`, `24`, and `23` paired years;
Eyrie has no matching photo series. The median photo-minus-manta cover
differences are `+0.017`, `-0.070`, and `-0.073`, with within-reef Pearson
correlations `0.51`, `0.59`, and `0.26`. Thus a pooled or universal
photo-to-manta conversion is not defensible. Manta is the common four-reef
**survey diagnostic**, and photo is an independent three-reef method check;
neither enters the frozen COTS score. The model output is area-weighted
whole-reef cover from all 2,895 V2 sites, each confirmed in the domain
GeoPackage to have placeholder median depth `5 m`;
the selected LTMP observations are `9 m`. Their year-matched differences are
diagnostic, not a calibrated cover likelihood or direct 9 m prediction.

Next promotion gate: obtain Owen matrix methods confirming whether source
water-quality mortality is absent or already embedded, and source-specific
2018–22 survival/area evidence (or run labelled bracketing scenarios). Then
test the area-weighted production hypothesis as an opt-in control/treatment
with flux and mass-balance checks. For coral, define a depth/habitat-aware
site-to-transect observation operator from actual survey locations or matched
depth strata, validate it against both LTMP methods without pooling them,
and only then consider a separately frozen coral score. Do not tune fecundity,
settlement, or maturation to the current support-mismatched cover residuals.

```powershell
$env:OWEN_FLUX_AUDIT_RUN_ID = 'choose_an_unused_flux_audit_id'
julia --project=sandbox sandbox\calibration\audit_owen_flux_contract.jl
python sandbox\calibration\test_owen_flux_audit.py --run-id choose_an_unused_flux_audit_id
$env:PYTHONDONTWRITEBYTECODE = '1'
python sandbox\calibration\audit_ltmp_coral_mapping.py --run-id choose_an_unused_coral_audit_id
```

#### Assumption-first Owen peak-capacity screen (2026-10-02)

We also ran the simpler, explicitly **exploratory** assumption that Owen's
COTS connectivity already includes larval mortality. The separate RME
multiplier is off in this treatment only; the default model and fixed
`3 COTS/ha` Allee threshold are unchanged. Merely turning survival off while
retaining the old `6389` pre-pelagic fecundity and `38828` boundary production
would increase effective supply by roughly three orders of magnitude. The
checked-in [`embedded-mortality screen`](runs/20261001T162000_owen_embedded_mortality_screen/metadata.toml)
therefore uses **one global, diagnostic mean-flux match**, not per-reef or
per-year fitting: effective fecundity `9.48` pre-pelagic recruits per adult
per year and upstream production `35.67` pre-pelagic recruits/ha/year. (The
exploratory metadata key ending `fecundity_pre_pelagic_ha_year` is a unit-label
misnomer; it stores the per-adult coefficient.) It then changes local
production and boundary supply separately and jointly, with the paired
gate's sampled hydrodynamic years, seed `20260930`, frozen CPUE scales, and
unchanged stage parameters. These numbers are not biological estimates.

| Treatment | Model peaks / 7 observed | Matched observed peaks | Mean frozen loss | What it shows |
| --- | ---: | ---: | ---: | --- |
| RME-on paired control | 5 | 4 | 4.591 | Original one-seed reference |
| Embedded, mean-flux matched | 5 | 4 | 4.591 | Removing the separate multiplier alone changes little |
| Embedded, local production 2× | 4 | 4 | 4.737 | Local production alone does not restore recurrence |
| Embedded, upstream boundary 2× | 8 | 6 | 4.139 | Two formal maxima on all four reefs |
| Embedded, local and boundary 4× | 8 | 7 | 2.923 | Seven matched peaks, but not the observed wave shape |

The two informative treatments were repeated with **the same parameters**
under three recorded hydrodynamic seeds in
[`20261001T163000_owen_embedded_replicates`](runs/20261001T163000_owen_embedded_replicates/metadata.toml).
Boundary 2× yielded eight model peaks and six matches in each seed; joint 4×
yielded eight model peaks and seven matches in each seed. This demonstrates
that the existing equations can make repeated maxima under the assumption.
It does **not** yet demonstrate distinct observed outbreaks. The
[peak-capacity plot](runs/20261001T162000_owen_embedded_mortality_screen/plots/cots_peak_capacity.png)
and [coral guardrail](runs/20261001T162000_owen_embedded_mortality_screen/plots/coral_guardrail.png)
show broad, elevated inter-peak populations. At seed `20260930`, joint 4×
peaks at Lizard are `0.52` and `0.66` COTS/tow versus observed smoothed
`1.42` and `0.74`; Eyrie peaks only reach about `0.11` versus its observed
`0.85`. A peak matcher can pair years despite those amplitude misses and an
extra Eyrie peak. Joint-4× 2015 whole-reef model coral-cover medians over
three seeds are `0.079`, `0.079`, `0.040`, `0.111` (Lizard, MacGillivray,
North Direction, Eyrie); the 9 m manta medians are `0.079`, `0.084`, `0.068`,
`0.105`. The methods' different spatial/depth support still prevents a
like-for-like coral calibration claim.

The post-hoc, **diagnostic only**
[inter-peak trough audit](runs/20261001T163000_owen_embedded_replicates/shape_audit/interpeak_troughs.csv)
defines trough depth as the minimum 3-year-smoothed COTS/tow between two
detected peaks divided by the smaller peak. The observed ratios are
`0–0.005` at Lizard, MacGillivray, and North Direction; boundary-2× model
ratios are `0.48–0.61` and joint-4× ratios `0.41–0.53` across seeds and
reefs. The current frozen COTS score does not sufficiently penalize this
persistent mid-cycle population; do not promote a treatment on peak count
alone or retroactively alter its reported score.

Two narrow diagnostics localize the remaining gap. An **oracle
counterfactual** setting only the external boundary to zero in 2000–10
([`20261002T090000_owen_boundary_shutdown`](runs/20261002T090000_owen_boundary_shutdown/metadata.toml))
reduces the four modeled trough ratios to `0.17–0.21`, but local recruitment
and maturation continue and no near-zero collapse appears. At the
configured (not independently evidence-validated) upper bounds, constant
adult mortality `m3=0.3` and Ricker density dependence `b_ricker=0.5`,
separately or jointly, give trough ratios about `0.30–0.47` in the
[`crash-parameter screen`](runs/20261002T094000_owen_crash_parameter_screen/metadata.toml).
Increasing mortality worsens peak matching; increasing density dependence
keeps seven matches but still leaves broad mid-cycle abundance. See the
[fixed-seed shape comparison](runs/20261002T094000_owen_crash_parameter_screen/plots/cots_peak_shape_diagnostics.png).

Conclusion: it was reasonable to proceed under a clearly labelled Owen
mortality assumption. Existing parameters can generate **recurrence**, and
upstream production matters far more than simply doubling local fecundity.
The tested constant parameter bounds and an upstream-off interval do **not**
reproduce the observed outbreak/crash contrast. This points to a missing or
mis-specified *decline/persistence* process, but does not prove a new core
mechanism is necessary: the time-varying boundary and other existing
parameters still need a discriminating test. Before another optimizer, freeze
an inter-peak trough/depletion diagnostic from LTMP (with observer uncertainty)
alongside timing, height, and coral guardrails. A small next factorial can
then separate sustained internal recruitment from adult survival and a
temporally intermittent upstream boundary. Keep all such treatments opt-in.

```powershell
$env:OWEN_EMBEDDED_RUN_ID = 'choose_an_unused_embedded_screen_id'
julia --project=sandbox sandbox\calibration\screen_owen_embedded_mortality.jl
python sandbox\calibration\plot_owen_embedded_mortality.py --run-id choose_an_unused_embedded_screen_id
python sandbox\calibration\audit_owen_peak_shape.py --run-id choose_an_unused_embedded_screen_id
$env:OWEN_EMBEDDED_REPLICATE_RUN_ID = 'choose_an_unused_replicate_id'
julia --project=sandbox sandbox\calibration\replicate_owen_embedded_mortality.jl
python sandbox\calibration\audit_owen_replicate_shape.py --run-id choose_an_unused_replicate_id
$env:OWEN_SHUTDOWN_RUN_ID = 'choose_an_unused_shutdown_id'
julia --project=sandbox sandbox\calibration\test_owen_boundary_shutdown.jl
$env:OWEN_CRASH_RUN_ID = 'choose_an_unused_crash_screen_id'
julia --project=sandbox sandbox\calibration\screen_owen_crash_parameters.jl
python sandbox\calibration\audit_owen_peak_shape.py --run-id choose_an_unused_crash_screen_id
python sandbox\calibration\plot_owen_peak_shape_diagnostics.py --run-id choose_an_unused_crash_screen_id
```

#### Inter-wave recruitment-blackout diagnostic (2026-10-02)

To test the missing crash directly, the protocol froze a deliberately
generous *screening* threshold before this experiment: on Lizard,
MacGillivray, and North Direction, the minimum three-year-smoothed modeled
COTS/tow strictly between two detected peaks must be at most **10%** of the
smaller peak, with seven observed-peak matches retained across all four reefs.
The observed ratios are `0-0.005`. This is not a formal observation-error
interval, does not replace the frozen score, and cannot by itself certify
realistic peak heights. Eyrie has only one observed peak and is excluded from
the trough ratio, not from timing, amplitude, or coral checks.

The opt-in adapter diagnostic `COTS_INTERNAL_SETTLEMENT_BLACKOUT_START_YEAR`
and `...END_YEAR` sets the reef-network settlement scalar to zero during the
inclusive **calendar** interval. It does not change the connectivity matrix,
fecundity, mortality, or external boundary supply. With both variables unset,
the archived joint-fourfold seed `20260930` replays adult densities, coral
cover, and selected hydrodynamic years to floating-point tolerance. The
[`20261002T120000_owen_recruitment_blackout` run](runs/20261002T120000_owen_recruitment_blackout/metadata.toml)
compares 2000-10 internal-only and external-only source blackouts, their
combination, and the combination at configured-upper adult mortality
`m3=0.3`. The 11-year blackout is an **oracle counterfactual**, not inferred
historical transport or production. Local and internal settlement and/or
external boundary flux were verified zero in the designated treatment years;
pre-existing juveniles still matured. The fixed Allee threshold remained
`3 COTS/ha`, Owen was assumed to embed larval mortality, and all treatments
used the same reef mapping and CPUE scales.

| Treatment | Inter-peak trough / smaller peak on the three two-wave reefs | Matched observed peaks / 7 | Mean frozen loss | Trough screen |
| --- | ---: | ---: | ---: | --- |
| Joint-fourfold reference | 0.405-0.518 | 7 | 2.923 | No |
| External supply off, 2000-10 | 0.169-0.198 | 7 | 3.261 | No |
| Internal settlement off, 2000-10 | 0.270-0.366 | 7 | 3.120 | No |
| Both sources off, 2000-10 | 0.051-0.071 | 6 | 3.522 | No: lost one match |
| Both off plus `m3=0.3` | 0.037-0.048 | 7 | 3.494 | Yes, screening only |

The [four-reef crash plot](runs/20261002T120000_owen_recruitment_blackout/crash_audit/cots_recruitment_blackout.png)
shows the observed raw tow points, three-year-smoothed observations, and
three-year-smoothed model trajectories on the frozen COTS/tow scale. The
[shape table](runs/20261002T120000_owen_recruitment_blackout/crash_audit/peak_shape.csv)
and [gap-flux table](runs/20261002T120000_owen_recruitment_blackout/crash_audit/gap_flows.csv)
give reef-level values. For example, Lizard's 2000-10 mean modeled internal
settled immigration was `0.381 recruits/ha/year` and external pelagic supply
was `0.442 pre-settlement recruits/ha/year` in the reference (settlement
probability and gate are both 1 in this screen). Either single-source blackout left
enough recruitment for an elevated trough; both together drove those fluxes
to zero, though maturation continued. The double-blackout's later Lizard
second peak (2020 rather than 2017) explains the lost match.

This demonstrates *crash capacity* under a long imposed recruitment gap, not
historical wave fidelity. The only treatment passing the new trough screen
uses a configured-upper, not evidence-validated mortality rate. Its mean
frozen loss worsens from `2.923` to `3.494`; Lizard's first modeled peak shifts
from 1997 to 2000 and falls from `0.52` to `0.41 COTS/tow`, compared with an
observed smoothed peak near `1.42`. Eyrie's modeled peak stays below `0.09`
against an observed `0.85`. Lizard 2015 modeled whole-reef coral cover rises
from `0.097` in the reference to `0.198`, whereas the 9 m manta median is
`0.079`; spatial/depth support remains mismatched. The model also declines
later than the observed first crash. Thus the treatment is **not promoted**
and does not justify a new core decline mechanism or fitted 11-year cutoff.
Next test evidence-based, time-varying local production/settlement and
upstream supply separately, while checking adult survival, peak amplitude,
decline timing, coral, and the frozen score across seeds.

The first run ID `20261002T113000_owen_recruitment_blackout` completed all
simulations but Julia's sandboxed Git subprocess failed during metadata
finalization. It was retained, then superseded by the complete run above;
their four output CSVs are byte-identical. To reproduce in a normal Julia
environment (choose an unused run ID):

```powershell
$env:OWEN_RECRUITMENT_BLACKOUT_RUN_ID = 'choose_an_unused_blackout_id'
julia --project=sandbox sandbox\calibration\screen_owen_recruitment_blackout.jl
$env:PYTHONDONTWRITEBYTECODE = '1'
python sandbox\calibration\audit_owen_recruitment_blackout.py --run-id choose_an_unused_blackout_id
```

If Julia cannot start `git` in a restricted runner, set
`OWEN_ADRIA_REVISION` and `OWEN_COTSMOD_REVISION` from each repository's
`git rev-parse HEAD` before the run; output and script hashes remain in
`metadata.toml`.

#### Upstream-source probability audit (2026-10-02)

The checked-in [source-response audit](audit_owen_source_probability_response.py)
reconstructs the currently used Owen upstream boundary from the global
source-to-sink annual matrices, reef areas, and
[`gbrPredsAdj_20262408.csv`](../data/gbrPredsAdj_20262408.csv). Its
[immutable run](runs/20261002T140000_owen_source_probability_response/metadata.json)
records source and output hashes, selected hydrodynamic years, source
contributions, and units. Reconstructed area-weighted boundary coefficients
match the archived forcing to `1.1e-16`; multiplying by the archived
`142.6624 pre-settlement recruits/ha/year` production matches logged external
pelagic flux to `1e-9`. The prediction file **does** mediate upstream supply,
but it is used *linearly* as source outbreak probability, not as an
outbreak-only switch.

As an offline test, the audit replaces each source probability `p` with
`p^2` or `p^4`; one global multiplier per alternative preserves the 2012-15
four-reef mean upstream supply. These are counterfactual response shapes,
not inferred fecundity or biological estimates. Mean 2004-09 incoming
pre-settlement flux (recruits/ha/year) is:

| Reef | Linear `p` | `p^2` | `p^4` |
| --- | ---: | ---: | ---: |
| Eyrie | 0.166 | 0.057 | 0.043 |
| Lizard | 0.229 | 0.076 | 0.061 |
| MacGillivray | 0.179 | 0.055 | 0.042 |
| North Direction | 0.379 | 0.178 | 0.172 |

This also reduces *first-wave* 1994-97 Lizard supply from `0.913` to
`0.667` (`p^2`) or `0.435` (`p^4`), and North Direction from `1.246` to
`0.921` or `0.616`. Since Lizard's first modeled peak is already far below
the smoothed observation, a steeper source response could worsen peak
amplitude. Clack Reef (14-017) alone contributes about `9-17%` of the
2004-09 boundary across the four sinks, with mean predicted outbreak
probability `0.457` in that window. The top five sources contribute about
`30-41%`; a steep response does not eliminate their supply. More decisively,
the full-model *external-zero oracle* still left `0.169-0.198` trough ratios,
above the frozen `0.10` screen. Thus a probability response alone is
unlikely to deliver the required trough; this is an inference from an
offline screen plus a stronger archived intervention, **not** a simulation
of the `p^2` or `p^4` treatments. No production map or default changed.
The reference still received mean 2004-09 internal settled immigration of
`0.060`, `0.204`, `0.255`, and `0.294 recruits/ha/year` at Eyrie, Lizard,
MacGillivray, and North Direction respectively; that pathway must be
examined independently of the upstream-probability map.

The existing model's separate food-survival pathway warrants a bounded test
before a new core mechanism. At roughly `0.055-0.063` modeled whole-reef
coral cover in 2004-09, the active `C_max` and `p_tilde` parameters can
change starvation survival, and hence juvenile maturation and adult
persistence. Test each configured upper bound separately and jointly under
the same forcing, seed, observation operator, and fixed `3 COTS/ha` Allee
threshold; report trough ratios alongside first/second peak heights, timing,
frozen loss, coral, and stage fluxes. Bounds are not evidence-validated and
passing a one-seed shape screen would not promote them. Do **not** tune the
exposed `eta_starve` or `fecundity_gate` factors: although passed through
ADRIA, neither appears in COTSMod's transition equations; the starvation
power is hard-coded as 3 and fecundity uses body condition squared. Any
attempt to make those factors active is a separate opt-in core-model change
requiring compatibility tests and the promotion gates.

```powershell
$env:PYTHONDONTWRITEBYTECODE = '1'
python sandbox\calibration\audit_owen_source_probability_response.py --run-id choose_an_unused_probability_audit_id
```

#### Existing food-survival parameter screen (2026-10-02)

The [bounded one-seed run](runs/20261002T150000_owen_food_survival_screen/metadata.toml)
replayed the archived joint-fourfold reference exactly (adult state and coral
to `1e-12`, hydrodynamic years exactly), then changed only `C_max` and/or
`p_tilde` to their configured upper bounds of `1.0`. Production, Owen
transport, upstream boundary, CPUE mapping, seed, frozen score, and the
`3 COTS/ha` Allee threshold were fixed. These bounds have **not** been
independently validated. The [shape and coral audit](runs/20261002T150000_owen_food_survival_screen/food_audit/peak_shape_and_coral.csv)
and [COTS/coral plot](runs/20261002T150000_owen_food_survival_screen/food_audit/cots_coral_food_survival.png)
show:

| Treatment | Lizard trough / peak | MacGillivray | North Direction | Observed peaks matched / 7 | Mean frozen loss |
| --- | ---: | ---: | ---: | ---: | ---: |
| Joint-fourfold reference | 0.481 | 0.405 | 0.518 | 7 | 2.923 |
| `C_max=1.0` | 0.502 | No two peaks | 0.407 | 6 | 3.678 |
| `p_tilde=1.0` | 0.486 | 0.409 | 0.519 | 7 | 2.931 |
| Both upper | 0.506 | No two peaks | 0.413 | 6 | 3.680 |

Increasing `C_max` did lower *absolute* Lizard trough COTS/tow from
`0.250` to `0.172`, but it also lowered the first peak from `0.519` to
`0.343`, delayed that peak from 1997 to 2000, and removed MacGillivray's
first formal peak. Lizard 2015 whole-reef coral cover rose from `0.097` to
`0.166`, against a 9 m manta median of `0.079`; these supports are not
directly equivalent, but the direction is unfavorable. Mean 2000-10
internal settled immigration fell from `0.381` to `0.081` at Lizard and
`0.489` to `0.104 recruits/ha/year` at MacGillivray, whereas external
pelagic supply remained `0.442` and `0.381`, respectively. Thus stronger
food-dependent survival suppresses the waves and internal production broadly,
not selectively enough to create a deep inter-wave crash. `p_tilde` alone
has little leverage near its fitted value `0.951`. None passes the frozen
qualitative trough screen; no setting or core default is promoted.

The next discriminating hypothesis is a *time-varying recruitment gap* that
is independently tied to source outbreak probabilities, larval survival,
settlement, or spatial connectivity, and suppresses both boundary and
internal supply during quiet years without weakening the first and second
wave supplies. The earlier 2000-10 double-blackout shows crash capacity but
is only an oracle bound; it is not a defensible historical schedule. A
small, prospectively defined forcing/settlement factorial should be checked
against peak height and timing, trough depth, coral, and the frozen score
before any optimizer or core mechanism change.

```powershell
$env:OWEN_FOOD_SURVIVAL_RUN_ID = 'choose_an_unused_food_screen_id'
julia --project=sandbox sandbox\calibration\screen_owen_food_survival.jl
$env:PYTHONDONTWRITEBYTECODE = '1'
python sandbox\calibration\audit_owen_food_survival.py --run-id choose_an_unused_food_screen_id
```

If a restricted Julia runner cannot launch `git`, set
`OWEN_ADRIA_REVISION` and `OWEN_COTSMOD_REVISION` from the two repositories'
`git rev-parse HEAD` values before the run; the metadata also hashes the
scripts, core/adapter source, forcing provenance, inputs, and outputs.

#### External-source maternal-condition proxy test (2026-10-02)

The user-suggested mechanism is biologically distinct from an arbitrary
recruitment blackout: low adult density reduces fertilisation, and depleted
coral prey may reduce the quality of larvae produced by surviving adults.
The **local** COTS source equation already contains the fixed Allee factor
and body-condition-squared fecundity factor. The external boundary instead
uses upstream outbreak probability multiplied by a constant production
scale. We do not have annual adult density and coral/maternal condition for
all 3,692 external sources, so the full common-source equation cannot yet be
estimated or tested for the entire GBR boundary.

The [opt-in source-proxy builder](build_owen_source_condition_proxy.py) tested
one explicitly provisional bridge: for each external reef `u`, set
`q[u,1991]=0.8` and update
`q[u,t]=(1-1/tau)q[u,t-1]+(1/tau)*(1-p[u,t-1])`, using the frozen local
candidate's `tau=8.333153725198475` years. Replace the source probability
weight `p[u,t]` with `p[u,t]*q[u,t]^2`. Here `1-p` is **not measured coral**;
it is a lagged source-condition proxy. A single global multiplier `1.58115`
preserves mean four-reef 2012-15 external flux. Owen matrices, source-to-sink
orientation, habitable-area conversion, and the embedded-pelagic-mortality
assumption are unchanged. The [builder run](runs/20261002T160100_owen_source_condition_proxy/metadata.json)
replays the archived linear boundary to maximum error `4.5e-16` and records
input/output hashes. External adult density and fertilisation remain
untested; probability was **not** substituted for adult density in the Allee
equation.

After normalisation, the proxy barely changes 2004-09 incoming flux
(pre-settlement recruits/ha/year): Lizard `0.229` to `0.236`,
MacGillivray `0.179` to `0.179`, and North Direction `0.379` to `0.359`.
First-wave 1994-97 flux falls slightly at all three. The paired full-model
[run](runs/20261002T170000_owen_source_condition_screen/metadata.toml)
uses the same seed `20260930`, local/internal COTS production equations,
fixed `3 COTS/ha` Allee threshold, observation mapping, and frozen score.
Its reference replays archived adult/coral trajectories to `1e-12` and
hydrodynamic years exactly. The
[shape/coral table](runs/20261002T170000_owen_source_condition_screen/source_condition_audit/peak_shape_and_coral.csv)
and [comparison plot](runs/20261002T170000_owen_source_condition_screen/source_condition_audit/cots_coral_source_condition_proxy.png)
show:

| Treatment | Lizard trough / peak | MacGillivray | North Direction | Peaks matched / 7 | Mean frozen loss |
| --- | ---: | ---: | ---: | ---: | ---: |
| Linear external boundary | 0.481 | 0.405 | 0.518 | 7 | 2.923 |
| Lagged external condition proxy | 0.505 | 0.433 | 0.488 | 7 | 2.838 |

Lizard's modeled first peak changes only `0.519` to `0.526 COTS/tow`, but
its trough rises `0.250` to `0.266`; MacGillivray's second peak falls
`0.382` to `0.355`. The apparent improvement in frozen mean loss does not
pass the separately frozen trough gate (`<=0.10` on all three reefs).
The proxy is **not promoted** and its weak effect must not be interpreted as
evidence against maternal nutrition generally. It only rejects this
particular `1-p` condition proxy, fixed lag, and late-wave normalisation as
an explanation for deep troughs. A stronger fitted penalty would be an
unvalidated forcing fit. The next substantive input is upstream
reef-specific adult density and coral/condition history (or a defensible
full-GBR source model), which would allow the same fertilisation and
maternal-viability equation to be evaluated on both sides of the Lizard
boundary.

```powershell
$env:PYTHONDONTWRITEBYTECODE = '1'
python sandbox\calibration\build_owen_source_condition_proxy.py --run-id choose_an_unused_proxy_builder_id
# Compare that builder output hash with the archived run. The paired runner
# deliberately reads the immutable archived proxy input recorded above.
$env:OWEN_SOURCE_CONDITION_RUN_ID = 'choose_an_unused_proxy_screen_id'
julia --project=sandbox sandbox\calibration\screen_owen_source_condition_proxy.jl
python sandbox\calibration\audit_owen_source_condition_proxy.py --run-id choose_an_unused_proxy_screen_id
```

If Julia cannot launch `git` in a restricted runner, set the two
`OWEN_*_REVISION` variables as above. Never overwrite a prior run directory.

### Local fecundity × coral consumption screen (2026-10-02)

The one-seed, four-treatment [screen](runs/20261002T180000_owen_fecundity_consumption_screen/metadata.toml)
tested whether larger local larval production and faster coral feeding can
produce sharper COTS waves through prey depletion, lower body condition, and
starvation. It kept the joint-fourfold Owen upstream boundary, connectivity,
fixed area-weighted observation conversion, frozen score, and `3 COTS/ha`
Allee threshold. Each treatment used the archived value or `1.5×` for local
fecundity, and independently the archived value or `1.5×` for both `a_F` and
`a_S`. The control replayed the archived adults, coral, and hydrodynamic years;
there were no failed runs. The `1.5×` perturbation is diagnostic, not a
biological bound. Moreover, the historic optimiser allowed `a_F<=2.0` but
the current factor table allows only `a_F<=1.0`, which the archived candidate
already exceeds. That contract must be resolved before any promotion.

The [peak/coral table](runs/20261002T180000_owen_fecundity_consumption_screen/fecundity_consumption_audit/peak_shape_and_coral.csv),
[flux table](runs/20261002T180000_owen_fecundity_consumption_screen/fecundity_consumption_audit/gap_flows.csv),
and [comparison plot](runs/20261002T180000_owen_fecundity_consumption_screen/fecundity_consumption_audit/cots_coral_fecundity_consumption.png)
show the result:

| Treatment | Detected peaks: Lizard / MacGillivray / North Direction | Matched observed peaks | Mean frozen loss | Trough gate |
| --- | --- | ---: | ---: | --- |
| Joint-fourfold reference | 2 / 2 / 2 | 7/7 | 2.923 | Fail: ratios 0.481 / 0.405 / 0.518 |
| Local fecundity 1.5× | 2 / 3 / 2 | 7/7 | 3.416 | Fail: extra peak; remaining ratios 0.406 / 0.422 |
| Consumption 1.5× | 2 / 3 / 2 | 7/7 | 3.194 | Fail: extra peak; remaining ratios 0.540 / 0.489 |
| Both 1.5× | 3 / 3 / 3 | 7/7 | 3.556 | Fail: extra peaks on all three reefs |

The `2004-09` mean modeled COTS/tow across the three two-wave reefs falls from
`0.224` (reference) to `0.148` (both), but this is not a deep trough and the
combined treatment shifts/adds early peaks. In `2000-10`, mean internal
immigration falls from `0.473` to `0.105` pre-settlement recruits/ha/year,
whereas external pelagic supply remains `0.477` under every treatment. This
helps explain why stronger local prey feedback alone does not create the
observed crash. Modeled coral also falls on three of four reefs in 2015; the
LTMP 9 m manta/photo measurements are shown as a guardrail, not a directly
equivalent whole-reef observation mapping. The tested parameter-only
explanation is **not promoted**. The next discriminating input remains an
independently supported upstream adult-density/condition history (or a
full-GBR source model), not a higher fecundity/feeding setting fitted to this
one seed.

The [BlackBoxOptim driver](calibrate_cots_blackbox.jl) was used for the earlier
v0.1 candidate that these screens replay; the current v0.2 Owen factorial is
**not** a black-box optimisation. The driver explicitly rejects v0.2 while
its qualitative mechanism and objective gates remain unresolved. Optimising
the existing loss now could reward matched peaks despite extra peaks, broad
troughs, and coral degradation.

```powershell
$env:OWEN_FECUNDITY_CONSUMPTION_RUN_ID = 'choose_an_unused_screen_id'
$env:OWEN_ADRIA_REVISION = (git rev-parse HEAD)
$env:OWEN_COTSMOD_REVISION = (git -C ..\COTSMod.jl rev-parse HEAD)
julia --project=sandbox sandbox\calibration\screen_owen_fecundity_consumption.jl
$env:PYTHONDONTWRITEBYTECODE = '1'
python sandbox\calibration\audit_owen_fecundity_consumption.py --run-id choose_an_unused_screen_id
```

### Low-cover COTS mortality trigger test (2026-10-02)

The [site-level diagnostic](runs/20261002T201000_owen_food_trigger_diagnostic/metadata.toml)
replayed the archived joint-fourfold Owen treatment and reconstructed the
existing food-survival multiplier from adult states and maturation, and from
body-condition/pre-feeding prey cover, to a maximum difference of `4.4e-15`.
The [trigger table](runs/20261002T201000_owen_food_trigger_diagnostic/trigger_audit/trigger_window_summary.csv)
and [timing plot](runs/20261002T201000_owen_food_trigger_diagnostic/trigger_audit/food_trigger_timing.png)
tested three fixed, uncalibrated triggers. Low condition and low prey cover
*per adult* were not selective for the trough; they also stressed the second
or first outbreak wave. A direct low-prey-cover trigger at `0.10` cover
fraction had trough stress `0.200/0.217/0.240` at Lizard, MacGillivray,
and North Direction—at least `2.84` times the larger wave stress. This
offline result warranted one coupled test, not a parameter estimate.

The opt-in core treatment multiplies the legacy food-survival factor for
juveniles and adults by `1 - lambda*max(0,1-(F+S)/0.10)`, where `F+S` is
pre-feeding fast-plus-slow coral cover. Its default switch is **off**.
Unlike raising `C_max`, it does not change the condition signal or the
immigration gate. The one-seed [dynamic run](runs/20261002T210000_owen_low_cover_mortality_screen/metadata.toml)
kept the external boundary, Owen hydrodynamic selection, local fecundity,
fixed `3 COTS/ha` Allee threshold, area-weighted observation mapping, and
frozen score identical across control, `lambda=0.5`, and `lambda=1.0`.
Control adults/coral replayed the archive to `1e-12`; no runs failed.
The [peak/coral table](runs/20261002T210000_owen_low_cover_mortality_screen/low_cover_audit/peak_shape_and_coral.csv),
[food/flux table](runs/20261002T210000_owen_low_cover_mortality_screen/low_cover_audit/gap_flux_and_food.csv),
and [comparison plot](runs/20261002T210000_owen_low_cover_mortality_screen/low_cover_audit/cots_coral_low_cover_mortality.png)
show:

| Treatment | Lizard trough / peak | MacGillivray | North Direction | Peaks matched / 7 | Mean frozen loss |
| --- | ---: | ---: | ---: | ---: | ---: |
| Control | 0.481 | 0.405 | 0.518 | 7 | 2.923 |
| Extra mortality `lambda=0.5` | 0.472 | 0.407 | 0.500 | 7 | 3.041 |
| Extra mortality `lambda=1.0` | No two peaks | 0.391 | 0.477 | 6 | 3.680 |

The stronger treatment removes Lizard's first formal peak, while neither
strength approaches the frozen `<=0.10` trough criterion. Mean effective
adult food survival across the three two-wave reefs in 2000-10 only falls
from `0.850` to `0.837` or `0.833`; local/internal recruitment declines,
but external pelagic supply remains `0.477` pre-settlement recruits/ha/year.
Higher mortality leaves more coral, which weakens the low-cover trigger in
the coupled model. This specific extra-mortality mechanism is **not
promoted**. The tested threshold and strengths are diagnostic rather than
evidence-derived, and whole-reef modeled coral is not directly equivalent
to LTMP 9 m manta/photo observations. The core switch remains off by
default; setting `COTS_LOW_COVER_MORTALITY=false` is the rollback.

```powershell
$env:OWEN_FOOD_DIAGNOSTIC_RUN_ID = 'choose_an_unused_food_diagnostic_id'
$env:OWEN_ADRIA_REVISION = (git rev-parse HEAD)
$env:OWEN_COTSMOD_REVISION = (git -C ..\COTSMod.jl rev-parse HEAD)
julia --project=sandbox sandbox\calibration\diagnose_owen_food_trigger.jl
python sandbox\calibration\audit_owen_food_trigger.py --run-id choose_an_unused_food_diagnostic_id

# The dynamic runner reads the immutable 20261002T201000 diagnostic gate.
$env:OWEN_LOW_COVER_RUN_ID = 'choose_an_unused_low_cover_screen_id'
julia --project=sandbox sandbox\calibration\screen_owen_low_cover_mortality.jl
python sandbox\calibration\audit_owen_low_cover_mortality.py --run-id choose_an_unused_low_cover_screen_id
```

### Owen upstream, initial-state, and low-boundary audit (2026-10-07)

The active experimental Owen dataset is `owen_global_2018_2023`, derived from
the five newer GBR-wide `conMatCots2018-19.csv` through
`conMatCots2022-23.csv` files under
`sandbox/data/OwenHydro/connectivityMatsYearlyGlobal`. Its 3,861-reef matrices
have **source rows and sink columns**; 113 reefs are internal to Lizard,
3,692 mapped reefs supply the outside boundary, and 56 reefs lacking the
required RME source match are omitted. The 2018-22 hydrodynamic matrices are
sampled over historical simulation years; they are not historical matrices
for 1985-2024. The outside boundary uses the supplied prediction probability
of `>0.22 COTS/tow` from 1991 onward, multiplied by a fixed absolute
pre-pelagic production scalar. Those outside reefs are **not** dynamically
simulated, so their density-dependent fertilisation and prey feedback do not
arise from the local model. Internal Lizard-to-Lizard production does use
modeled densities. Separate source survival is off for this Owen treatment,
on the stated assumption that the matrix embeds larval loss.

The read-only [initialization and source audit](runs/20261007T101500_owen_initialization_audit/metadata.json)
verified all 23 provenance-recorded source and derived-file hashes. The
independent upstream reconstruction in
`20261002T140000_owen_source_probability_response` matched the boundary to
`1.1e-16`; dominant trough suppliers include Clack, Combe, Conical Rock,
Pickersgill, and Corbett (the order varies by sink). The remaining policy
uncertainties are the 56 excluded source reefs and whether Owen connectivity
already embeds mortality; neither was varied in the screen below.

The frozen reef-specific observation factor converts modeled **adult COTS/ha**
to expected COTS/tow. It is a fitted observation mapping, not a universal tow
area. Thus 0.1-0.4 COTS/tow corresponds to approximately 0.34-1.36 adults/ha
at Lizard, 0.61-2.44 at MacGillivray, 0.79-3.18 at North Direction, and
1.54-6.18 at Eyrie. The supported 3 adult COTS/ha Allee threshold is
equivalent to 0.883, 0.492, 0.378, and 0.194 COTS/tow respectively under
those *observation* scales. The actual archived spatial seed was below this
threshold at every target site:

| Reef | Initial adults/ha | Subadults/ha | Juveniles/ha | Model 1989 COTS/tow |
| --- | ---: | ---: | ---: | ---: |
| Lizard | 0.337 | 1.010 | 2.019 | 0.341 |
| MacGillivray | 0.389 | 1.168 | 2.335 | 0.226 |
| North Direction | 0.100 | 0.299 | 0.599 | 0.073 |
| Eyrie | 0.357 | 1.072 | 2.145 | 0.074 |

The seed is determined by the domain's spatial `cots_density` field and the
archived multiplier `3.915`. Crucially, that field was sampled from
`COTS_prob_0.02_cpue_year2025_clean.tif` during the V2 build, although this
simulation begins in **1985**. Its source hash still matches the V2 domain
provenance. This is a 2025 spatial probability proxy, not a demonstrated
1985 abundance distribution: `2,881/2,895` sites and all 113 reefs exceed
the initializer's `p>0.1` seed cutoff. The spatial vector bypasses the
`COTS_SEED_FIRST_N`/`ADRIA_DEBUG_INIT_DENSITY` fallback settings. The
initialization formula preloads juvenile,
subadult, and adult classes in a **6:3:1** ratio. At the four target reefs,
the initial adult Allee multiplier is only `0.0035-0.0166`; early 1986-90
adult abundance largely comes from the preloaded subadult/juvenile pipeline
and within-domain dynamics. The 1985 COTS log row is currently an unfilled
zero row because the scenario loop starts at timestep 2. Treating it as the
biological initial state or smoothing through it is misleading. External
prediction supply is zero before 1991, so it cannot explain the premature
1989-90 modeled CPUE at Lizard and MacGillivray.

The preregistered [2 x 3 diagnostic run](runs/20261007T103000_owen_seed_boundary_screen/metadata.json)
changed only initial seed multiplier (reference or one-quarter) and outside
production (100%, 10%, or zero of `142.662` pre-pelagic recruits/ha/year).
Local fecundity stayed `37.927` pre-pelagic recruits/adult/year, the adult
Allee threshold stayed `3/ha`, and hydrodynamics, fixed observation mapping,
scoring, and coral model stayed unchanged. The full/full replay matched the
archived adult and coral trajectories exactly.

| Seed / outside supply | Observed peaks matched (of 7) | Two-wave trough ratio: Lizard / Mac / North Direction | Interpretation |
| --- | ---: | --- | --- |
| Reference / 100% | 7 | 0.481 / 0.405 / 0.518 | Peaks retained; troughs far too shallow |
| Reference / 10% | 3 | no second peaks | Initial cohort declines; no convincing resurgence |
| Reference / zero | 3 | no second peaks | No external rescue under current parameters |
| Quarter / 100% | 6 | 0.387 / 0.366 / 0.454 | First Lizard/Mac peaks shift 1997 to 2002 |
| Quarter / 10% | 5 | 0.693 / 0.693 / 0.713 | Two very low, late waves; troughs worse |
| Quarter / zero | 0 | no formal peaks | Near-extinction |

See the [full COTS/coral comparison](runs/20261007T103000_owen_seed_boundary_screen/seed_boundary_audit/cots_coral_seed_boundary.png),
[first-wave zoom](runs/20261007T103000_owen_seed_boundary_screen/seed_boundary_audit/first_wave_zoom.png),
and [per-reef peak and flux table](runs/20261007T103000_owen_seed_boundary_screen/seed_boundary_audit/peak_shape_flux_and_coral.csv).
At Lizard in 2000-10, mean local fecundity/internal immigration fell from
`1.437/0.381` with the full boundary to `0.197/0.022` with 10% boundary and
`0.081/0.004` with zero boundary (pre-settlement recruits/ha/year), while
external pelagic supply fell `0.442 -> 0.044 -> 0`. Model coral cover in
2015 rose from `0.097` to `0.333` and `0.422` respectively, reflecting much
less COTS predation, although whole-reef model cover is not directly equivalent
to LTMP 9 m manta/photo observations.

No factorial cell passed the qualitative peak/trough/coral gate. This is a
one-seed screen of the **current** parameterization, not proof that all
possible low-boundary configurations fail. Simply setting the current adult
seed near 3/ha by multiplying it about ninefold would also preload roughly
18 juvenile and 9 subadult COTS/ha because of the hard-coded stage ratio;
that would not test the user's adult-threshold hypothesis cleanly. The next
bounded test first needs a historically defensible **1985 seed footprint**,
then an opt-in, independently specified initial stage profile with adults
around the supported threshold but a small immature cohort, followed by a
local-production/transport sensitivity under zero and low external boundary.
Keep default initialization and the 3/ha Allee
threshold unchanged until that treatment passes replicated peak and coral
gates. The exact design and failure criteria are in
`CALIBRATION_PROTOCOL.md`.

```powershell
$env:OWEN_ADRIA_REVISION = (git rev-parse HEAD)
$env:OWEN_COTSMOD_REVISION = (git -C ..\COTSMod.jl rev-parse HEAD)
$env:OWEN_SEED_BOUNDARY_RUN_ID = 'choose_an_unused_seed_boundary_id'
julia --project=sandbox sandbox\calibration\screen_owen_seed_boundary.jl
python sandbox\calibration\audit_owen_seed_boundary.py --run-id choose_an_unused_seed_boundary_id

python sandbox\calibration\audit_owen_initialization.py --run-id choose_an_unused_initial_audit_id
```

The 2026-10-07 screen used a metadata recovery script after all six
simulations had saved successfully; the only runner correction was to flatten
the treatment list for TOML serialization. Its recovered manifest includes
exact control-replay checks, source hashes, and both repository revisions.

### Connectivity and validation utilities

```powershell
# Build the static control plus annual reef matrices, source survival, and the
# external-GBR boundary coefficients. Select q3baseline, q3A, q3B, or q3R.
$env:COTS_WATER_QUALITY_SCENARIO = 'q3baseline'
julia --project=sandbox sandbox\domain_building\build_lizard_cots_connectivity.jl
julia --project=sandbox sandbox\domain_building\validate_lizard_cots_connectivity.jl

# Compare same-seed coral and COTS-connectivity pilots
julia --project=sandbox sandbox\calibration\analyze_connectivity_factorial.jl

# Freeze the pre-2006 observation scale and score 2006 onward
$env:BBO_RUN_ID = 'expanded_cotsconn_pilot_seed20260930'
julia --project=sandbox sandbox\calibration\validate_temporal_holdout.jl

# Diagnose the post-outbreak age structure and Allee multiplier
julia --project=sandbox sandbox\calibration\diagnose_recurrence.jl

# Run the bounded recurrence-mechanism factorial after building annual forcing
$env:BBO_RUN_ID = 'expanded_cotsconn_pilot_seed20260930'
$env:COTS_FACTORIAL_SEED = '20260929'
$env:COTS_FACTORIAL_EXTERNAL_SOURCE_DENSITY = '3.0' # pre-dispersal recruits ha^-1 yr^-1
$env:COTS_FACTORIAL_RUN_ID = 'recurrence_mechanism_repeat_seed20260929' # choose an unused ID
julia --project=sandbox sandbox\calibration\run_recurrence_mechanism_factorial.jl

# Reproduce the separated-production/evidence-boundary factorial
$env:BBO_RUN_ID = 'expanded_cotsconn_pilot_seed20260930'
$env:COTS_FACTORIAL_SEED = '20260929'
$env:COTS_FACTORIAL_RUN_ID = 'separated_production_evidence_repeat_seed20260929' # choose an unused ID
julia --project=sandbox sandbox/calibration/run_separated_production_factorial.jl
```

The factorial external source is an explicit screening value in pre-dispersal
recruits/ha/year, not an adult-density estimate. Repeat the bounded treatment
only with independently justified upstream production before promotion. The
`sample` connectivity mode is deterministic for a recorded
`COTS_CONNECTIVITY_SEED`; `cycle` is the primary
same-seed comparison.

Promotion criteria, the locked observed-peak contract, expanded-search order,
sensitivity outputs, connectivity experiment, and validation design are defined
in `CALIBRATION_PROTOCOL.md`.

### 1991 LTMP-initialization and coral-surface screen (2026-10-07)

This **experimental** screen tests a cleaner historical start; it does not
replace the 1985-2024 control or promote a COTS/tow-to-density conversion.
[`build_1991_initial_surfaces.py`](build_1991_initial_surfaces.py) pools
1989-1991 LTMP COTS counts/tows and interpolates them to all 2,895 V2 sites
with tow-effort-weighted IDW. It independently interpolates 9 m LTMP manta
hard-coral cover from the same years; no matching photo-transect observations
exist in that window. The builder uses a 111 km radius, at least three distinct
donor reefs, a 1 km distance floor, and explicit reef-ID/name crosswalks. The
[surface provenance](runs/20261007T199104_1991_initial_surfaces/metadata.json)
records 85 matched donors for each variable, zero uncovered sites, source and
output hashes, unmatched names, and leave-one-reef-out checks. The surface's
median is `0.00545 COTS/tow` and `0.225` 9 m manta coral-cover fraction.
The 1991 probability of `>0.22 COTS/tow` is joined by parent reef ID and
recorded separately; it is not interpreted as expected COTS/tow or multiplied
into the seed. The existing Owen outside-boundary treatment continues to use
the probability time series.

The [bounded 1991 run](runs/20261007T199105_owen_1991_initialization/metadata.toml)
uses the actual 1991-2024 historical DHW slice and the same Owen sample-mode
hydrodynamics (seed `20260930`), `3 adult COTS/ha` Allee threshold, fecundity,
outside boundary, fixed reef observation scales, and COTS score. The first
1991 COTS log row is a placeholder, so it and the 1989-1991 seed observations
are excluded from scoring. The opt-in `COTS_INITIAL_STATE_CSV` adapter takes
site-keyed juvenile, subadult, and adult **COTS/ha**, leaving the default
probability initializer unchanged. The IDW treatments seed adults only rather
than inheriting the unsupported 6:3:1 immature preload. Since a physical
COTS/tow-to-adult-density conversion remains unavailable, `0.015` and `0.150`
COTS/tow per adult COTS/ha are clearly labelled **sensitivity assumptions**,
not calibrated values. The coral treatment rescales each site's existing
group/size composition to the interpolated manta-cover fraction. It is a
whole-reef initial-cover *proxy* from 9 m surveys, not an established
observation model.

| Treatment | First Lizard / MacGillivray modeled peak | Max trough / smaller peak across three two-wave reefs | Matched observed peaks | Sum of frozen losses, 1992-2024 |
| --- | --- | ---: | ---: | ---: |
| Archived 1985 reference, re-scored on common window | 1997 / 1997 | not comparable after window truncation | 6 | 11.540 |
| 1991 inherited seed and coral | 1999 / 1999 | 0.313 | 6 | 13.795 |
| 1991 IDW adults, `0.015`, inherited coral | 2002 / 2003 | 0.453 | 6 | 13.928 |
| 1991 IDW adults, `0.015`, IDW coral | 2002 / 2003 | 0.590 | 6 | 13.605 |
| 1991 IDW adults, `0.150`, inherited coral | 2003 / 2003 | 0.515 | 6 | 13.953 |
| 1991 IDW adults, `0.150`, IDW coral | 2003 / 2003 | 0.492 | 6 | 13.943 |

See the [COTS/coral trajectories](runs/20261007T199105_owen_1991_initialization/cots_coral_1991_initialization.png),
[per-reef peak/trough audit](runs/20261007T199105_owen_1991_initialization/peak_trough_audit.csv),
and [common-window archived-reference score](runs/20261007T199105_owen_1991_initialization/historical_reference_common_window_scores.csv).
The observed first smoothed peaks are 1997 at Lizard and MacGillivray and
1994 at North Direction. The 1991 IDW treatments therefore delay rather than
recover the first outbreak. All three two-wave reefs remain above the
preregistered trough ceiling `0.10`. The coral initialization moves starting
cover toward local manta estimates, but simulated coral falls to roughly
`0.01-0.05` at the target reefs, with no convincing recovery of the observed
coral trajectories. Manta 9 m versus modeled whole-reef support remains a
material limitation. No 1991 treatment passes the qualitative gate; do not
optimize or promote it on this evidence. The result suggests that the first
wave's timing depends on pre-1991 history/immature cohorts and/or recruitment
structure, rather than just the choice of 1991 adult/coral surface. This is an
inference from the bounded screen, not proof of a unique missing mechanism.
The [post-change 1985 control replay](runs/20261007T199106_legacy_replay/metadata.toml)
passed its archived adult/coral equality checks to `1e-12` and exact
hydrodynamic-year alignment, with no failed treatments. The opt-in initializer
therefore has not changed the legacy default path in this test.

```powershell
python -m unittest discover -s sandbox/calibration -p test_1991_initial_surfaces.py -v
python sandbox/calibration/build_1991_initial_surfaces.py --output sandbox/calibration/runs/<new_surface_id>
$env:OWEN_1991_SURFACE_DIR = (Resolve-Path sandbox/calibration/runs/<new_surface_id>).Path
$env:OWEN_1991_RUN_ID = '<new_1991_run_id>'
julia --project=sandbox sandbox/calibration/test_1991_initial_state.jl
julia --project=sandbox sandbox/calibration/run_1991_initialization_test.jl
$env:OWEN_1991_RUN_DIR = (Resolve-Path sandbox/calibration/runs/<new_1991_run_id>).Path
julia --project=sandbox sandbox/calibration/score_1985_reference_on_1991_window.jl
python sandbox/calibration/audit_1991_initialization.py sandbox/calibration/runs/<new_1991_run_id>
python sandbox/calibration/plot_1991_initialization_test.py sandbox/calibration/runs/<new_1991_run_id>
```

### 1991 cohort history, remnant grazing, and conversion-data gate (2026-10-07)

The [bounded cohort run](runs/20261007T199111_owen_1991_cohort_history_final/metadata.toml)
holds the 1991 inherited coral, Owen connectivity/boundary, ecological
parameters, seed `20260930`, area-weighted observation mapping and frozen
score fixed. It varies only the opt-in 1991 initial COTS stage profile. The
`6:3:1` juvenile:subadult:adult ratio is inherited from the **model
initializer**, not observed historical stage densities. `0.015` and `0.150`
COTS/tow per adult COTS/ha are provisional observation sensitivities, not
validated conversions. The prior adult-only `0.015` run is the zero-immature
reference. See [cohort and coral plot](runs/20261007T199111_owen_1991_cohort_history_final/cots_coral_cohort_history.png),
[peak audit](runs/20261007T199111_owen_1991_cohort_history_final/peak_trough_audit.csv),
and [gap/coral audit](runs/20261007T199111_owen_1991_cohort_history_final/cohort_history_audit.csv).

| 1991 state | Lizard first peak | MacGillivray first peak | Worst two-wave trough / smaller peak | Observed peaks matched / 7 |
| --- | ---: | ---: | ---: | ---: |
| Inherited control | 1999 | 1999 | 0.313 | 6 |
| IDW adults only, `0.015` | 2002 | 2003 | 0.453 | 6 |
| IDW adults + half 6:3:1 cohort, `0.015` | 2001 | 2002 | 0.501 | 6 |
| IDW adults + full 6:3:1 cohort, `0.015` | **1997** | 1999 | 0.594 | 6 |
| IDW adults + full 6:3:1 cohort, `0.150` | 2004 | 2003 | 0.485 | 5 |

Restoring the full immature cohort at `0.015` recovers Lizard's observed 1997
first-peak year, but makes the trough **shallower**, leaves North Direction
late, and can drive its modeled coral to zero. This supports sensitivity to
unobserved pre-1991 cohorts, not this profile's promotion. The full 1985
model warm-up still outperforms a 1991 adult-only restart on the common scoring
window; a historical juvenile/subadult state has not yet been reconstructed.

The corrected *zero-COTS* counterfactual starts from inherited 1991 coral,
sets initial COTS and both internal/external production to zero, and verifies
zero COTS at every timestep. By 2012, Lizard coral is `0.637` versus `0.130`
in the inherited COTS run, showing a large aggregate COTS effect but not
identifying the inter-wave contribution by itself. A first attempt
([invalid diagnostic run](runs/20261007T199107_owen_1991_cohort_history/metadata.toml))
set `ADRIA_COTS_ENABLED=false` alone; the experimental external boundary
still injected COTS. It is retained as a failed diagnostic and **must not** be
used as a coral-only comparison.

The [post-peak grazing diagnostic](runs/20261007T199109_owen_1991_gap_grazing/metadata.toml)
keeps COTS demographic updates and recruitment intact but applies predation
to a copy of coral cover during 2004-2012. The opt-in window is default-off;
the inherited control replays the prior adult/coral trajectories exactly, and
the two paths are identical through 2003. See the [paired plot](runs/20261007T199109_owen_1991_gap_grazing/cots_coral_gap_grazing.png)
and [audit](runs/20261007T199109_owen_1991_gap_grazing/gap_grazing_audit.csv).
Lizard coral in 2012 rises from `0.130` to `0.482` when gap grazing is bypassed,
so remnant grazing materially constrains modeled recovery. Yet Lizard's
trough/smaller-peak ratio worsens from `0.286` to `0.451`, and its second peak
falls from `0.610` to `0.546 COTS/tow`. More coral alone does not generate a
deeper COTS crash or a larger next wave; COTS remain present because the
recruitment/demographic pathways remain active. This is an intentionally
unphysical *diagnostic*, not a candidate mechanism. The 9 m manta versus
whole-reef coral support difference also remains unresolved.

For a defensible conversion, [RosettaCOTS](https://github.com/Scott-Foster/RosettaCOTS)
cannot directly supply the LTMP count-to-density factor: its manta-tow
calibration excluded observed COTS counts. The supplied SALAD and EOTR manta/
cull workbooks are therefore being treated as independent calibration data,
not a ready-made factor. The [checked-in overlap audit](audit_cots_conversion_overlap.py)
and [provenance](runs/20261007T199110_conversion_overlap/metadata.json)
found **14** named-reef/year SALAD–EOTR manta overlaps in 2021-2026, **9**
with surveys within 30 days, and **13** with cull records. Six of the nine
near-date groups have nonzero manta COTS counts. The audit preserves the raw
workbooks, checks SALAD's count/5 m swath-area density calculation, and makes
the reef-ID crosswalk explicit. [Reef-year counts and densities](runs/20261007T199110_conversion_overlap/reef_year_overlap.csv)
are an *overlap inventory*, not conversion estimates. SALAD `Density(Ha)`
includes all detected sizes; EOTR and LTMP manta are treated as comparable
for this analysis, but culling is targeted,
the Lizard/South Direction groups each pool several reef IDs, and matching by
reef/year does not establish the same tow path or detection probability. The
original [SALAD size audit](runs/20261007T199112_salad_stage_sizes/metadata.json)
finds 991 individually recorded COTS in these 22 reef-years, all with sizes.
Its >=250/260 mm counts are archival sensitivities, not the agreed adult
definition. The [revised audit](runs/20261008Tsalad_adult150_sizes_checked/metadata.json)
uses the agreed categories: 985 are adults >=150 mm, of which 56 are small
adults 150-250 mm inclusive and 929 are larger adults >250 mm. Track-level and individual
counts reconcile exactly in only 13 of 22 reef-years, and in only 6 of the 9
near-date manta matches. Only 3 of those 6 have nonzero manta counts. The
[revised size audit rows](runs/20261008Tsalad_adult150_sizes_checked/salad_stage_sizes.csv)
therefore report >=150 mm adult COTS/ha only for exact reconciliations; they
do not justify a fitted tow-to-density coefficient yet. The 250/260 mm
cutoffs remain labelled archival comparators, not an adult definition.
The next observation-model gate is to resolve the track/individual count
mismatches, use SALAD per-animal sizes to define an adult-density outcome,
pair surveys as closely as spatially/temporally
possible, explicitly model zero-heavy tow counts and tow distance, and
propagate uncertainty in detection and tow-distance despite treating EOTR
and LTMP manta counts as comparable for this calibration. The
RosettaCOTS cull-to-SALAD relationship can provide an additional noisy bridge
only if its outcome definitions (COTS **plus scars**) and effort/selection
match; do not chain point estimates as if they were direct measurements.

The [coordinate proximity screen](runs/20261007T199115_salad_manta_spatial_overlap_final/metadata.json)
uses a 1 km SALAD-track/EOTR-tow midpoint bound and at most 30 days. It finds
40 SALAD tracks with 106 recorded COTS near 153 **distinct** manta tows, but
none of those tows saw a COTS. A [stricter 0.5 km/7-day screen](runs/20261007T199116_salad_manta_spatial_overlap_strict_final/metadata.json)
still has 14 SALAD tracks with 51 COTS near 30 distinct zero-count tows.
The [track-level audit](runs/20261007T199115_salad_manta_spatial_overlap_final/track_overlap.csv)
is a midpoint/date *screen*: it does not prove path intersection, identical
habitat or temporal population stability, and a tow can match multiple SALAD
tracks. It therefore constrains the observation-model design but is not a
direct estimate of detectability or `COTS/tow per COTS/ha`. Do not calibrate
the provisional `0.015`/`0.150` factors to this inventory.

No treatment passes the qualitative trough gate or changes the model default.
Before further optimization we need a supported adult COTS/ha observation
mapping and either observed pre-1991 cohorts or a labelled model warm-up
state; then repeat the bounded peak/trough test with held-out reefs/years.

```powershell
$env:OWEN_1991_EXPERIMENT = 'cohort_history'
$env:OWEN_1991_RUN_ID = '<unused_cohort_run_id>'
julia --project=sandbox sandbox/calibration/test_1991_initial_state.jl
julia --project=sandbox sandbox/calibration/run_1991_initialization_test.jl
python sandbox/calibration/audit_1991_initialization.py sandbox/calibration/runs/<unused_cohort_run_id>
python sandbox/calibration/audit_1991_cohort_history.py sandbox/calibration/runs/<unused_cohort_run_id>
$env:OWEN_1991_EXPERIMENT = 'gap_grazing'
$env:OWEN_1991_RUN_ID = '<unused_gap_run_id>'
julia --project=sandbox sandbox/calibration/run_1991_initialization_test.jl
python sandbox/calibration/audit_1991_initialization.py sandbox/calibration/runs/<unused_gap_run_id>
python sandbox/calibration/audit_1991_gap_grazing.py sandbox/calibration/runs/<unused_gap_run_id>
python sandbox/calibration/audit_cots_conversion_overlap.py --output sandbox/calibration/runs/<unused_overlap_run_id>
python sandbox/calibration/audit_salad_stage_sizes.py --overlap sandbox/calibration/runs/<unused_overlap_run_id>/reef_year_overlap.csv --output sandbox/calibration/runs/<unused_size_run_id>
python sandbox/calibration/match_salad_manta_tracks.py --output sandbox/calibration/runs/<unused_spatial_run_id>
python sandbox/calibration/match_salad_manta_tracks.py --output sandbox/calibration/runs/<unused_strict_run_id> --max-days 7 --max-km 0.5
```

### Coherent-conversion cycle-capacity screen and adult size (2026-10-08)

The user's observation contract is now **adult >=15 cm**, comprising small
adults 15-25 cm inclusive and larger adults >25 cm. EOTR and LTMP manta-tow
counts are treated as comparable for this calibration. There are no observed
historical immature cohorts, so a 1991 6:3:1 juvenile:subadult:adult preload
is a sensitivity scenario, not a reconstruction. The revised SALAD inventory
above uses these categories but does not resolve tow detectability.

[Babcock, Milton & Pratchett (2016)](https://doi.org/10.1007/s00227-016-3009-5)
report a fitted female gonad-mass relationship of approximately
`G_f(D_mm) = 3.384 exp(0.0115 D_mm)` grams (`R^2 = 0.55`, 103 females), and
a comparable male fit. On this curve a 40 cm animal has approximately 10 times
the female gonad mass of a 20 cm animal. Their oocyte estimates use 90,000
oocytes per gram of female gonad, with substantial individual scatter and
uncertainty about mature egg release. The paper also gives a power-law fit of
oocyte production to *body mass*; exponential-in-diameter and power-in-mass
are not interchangeable parameterizations. The current COTS adult state has
no within-adult size distribution: its `body_condition^2` fecundity factor
does not reproduce this size effect. A future opt-in mechanism should track
small and large adults or a size distribution, apply size-specific gamete
output **before** the fertilisation term, and test growth/food dependence.
With no cohort observations, assumed maturation and growth histories must be
bounded sensitivities. Using `G(E[D])` in place of `E[G(D)]` would also
understate reproductive output when adult sizes vary.

Before adding that mechanism, the [paired capacity runner](screen_1991_cycle_capacity.jl)
tests whether existing parameters can produce the missing two-peak/deep-trough
shape. This is not a production optimiser. It uses 24 deterministic
space-filling design points plus an archived-state control, each at the
supported `theta=3 COTS/ha` and a separately labelled `theta=1` historical
counterfactual. Design axes are one provisional global `q=0.005-0.30`
`COTS/tow per adult COTS/ha` factor; 1991 immature preload 0, half or full;
inherited or IDW coral; fecundity 0.5-2x; external boundary 0-1.5x; adult
mortality 0.08-0.30/year; and feeding 0.6-1.5x. These are **diagnostic ranges**,
not validated biological bounds. The same `q` maps IDW CPUE to 1991 adult
density and modeled adult density back to CPUE, fixing the earlier
initialization/observation-scale inconsistency *within this screen*. A global
`q` is intentionally restrictive; it is not a fitted SALAD-manta conversion.
The frozen 3-year-smoothed observed peaks and scoring function are used, but
losses from this coherent-`q` screen are not directly comparable with old
reef-specific-scale losses. All stochastic draws use recorded seeds and Owen
hydrodynamic years; the 3 COTS/ha anchor checks the archived adult, coral and
forcing-year trajectory to numerical tolerance.

The [independent audit](audit_1991_cycle_capacity.py) requires exact 2/2/2/1
peak counts, all seven observed peaks matched, trough/minimum-adjacent-peak
ratio <=0.10 at Lizard, MacGillivray and North Direction, all peak years
within four years and all heights within 0.5-2 times the observations. Coral
minimum and 2015 cover and local/internal/external recruitment during 2004-09
are reported separately; whole-reef modeled coral and 9 m manta coral remain
different measurement supports. Passing this qualitative screen would permit
replication and stricter coral/hold-out validation, **not** promotion. If only
`theta=1` passes, test the size-weighted breeder hypothesis before changing
the supported threshold. If neither passes, expand the diagnostic design or
test an independently motivated time-varying recruitment/size structure;
do not infer impossibility from 24 space-filling points. The legacy/default
model is untouched, and all experimental files can be ignored without a
rollback migration.

```powershell
$env:OWEN_CAPACITY_PAIRS = '24'
$env:OWEN_CAPACITY_RUN_ID = '<unused_capacity_run_id>'
julia --project=sandbox sandbox/calibration/screen_1991_cycle_capacity.jl
python sandbox/calibration/audit_1991_cycle_capacity.py sandbox/calibration/runs/<unused_capacity_run_id>
```

The [2026-10-08 attempted 24-pair run](runs/20261008Tcycle_capacity_24pairs/incomplete_provenance.json)
completed the control plus 23 design pairs (48 treatments). The final pair
had `q=0.00849`, half the inherited immature preload, and low adult mortality;
its `theta=3` run exceeded ten minutes of sustained CPU and was stopped before
either threshold produced a result. Both treatments for pair 24 are marked
**untested**, not failures of the ecological hypothesis. The completed
treatments have no logged exception and passed the independent finite-value,
unique-year, and four-reef completeness checks. The control's `theta=3`
adult, coral and selected hydrodynamic-year trajectories replayed the archived
1991 cohort run to `1e-12`. See the [design](runs/20261008Tcycle_capacity_24pairs/design.csv),
[per-reef audit](runs/20261008Tcycle_capacity_24pairs/capacity_audit.csv),
[treatment summary](runs/20261008Tcycle_capacity_24pairs/capacity_treatment_summary.csv),
and [audit metadata](runs/20261008Tcycle_capacity_24pairs/capacity_audit.json).
The [paired COTS/coral trajectory figure](runs/20261008Tcycle_capacity_24pairs/cots_coral_capacity_pair_020.png)
shows the lowest-loss supported-threshold completed treatment (`pair 20`),
its `theta=1` counterfactual, and the archived-state control against raw and
3-year-smoothed COTS/tow and LTMP manta coral. Its first modeled peaks are
late on Lizard, MacGillivray and North Direction, modeled mid-cycle COTS
persist, and modeled coral falls much lower than the 9 m manta observations
around the first wave. Coral supports differ, but this is not an acceptable
peak-and-coral calibration.

| Allee branch | Completed treatments | Maximum matched / 7 peaks | Best single-reef trough ratios: Lizard / MacGillivray / North Direction | Full qualitative passes |
| --- | ---: | ---: | --- | ---: |
| Supported `3 COTS/ha` | 24 | 6 | 0.293 / 0.229 / 0.254 | 0 |
| `1 COTS/ha` counterfactual | 24 | 6 | 0.287 / 0.277 / 0.266 | 0 |

The three *best single-reef* ratios in each row can come from different
parameter settings; they do not describe one successful treatment. No
completed treatment reached the <=0.10 trough criterion even on **one** of
the three two-wave reefs. Only one treatment in each branch had exactly the
observed 2/2/2/1 peak counts, and neither matched all seven peaks. Across
paired reef comparisons that had two detectable peaks at both thresholds,
reducing `theta` changed the median trough ratio by only `-0.005`. In the
2004-09 gap the median internal settled immigration was 0.587 recruits/ha/year
at `theta=3` and 1.720 at `theta=1` (paired median ratio `2.61`); median
external settlement was 0.149 recruits/ha/year in both. Thus lowering the
local fertilisation half-saturation mainly amplified within-domain larval
supply without solving the adult-density floor. Across the 96 reef-treatment
comparisons in each branch, the median modeled 2015 coral fraction fell from
0.077 at `theta=3` to 0.055 at `theta=1`; this is a whole-reef diagnostic,
not a direct LTMP 9 m coral fit.

This screen **does not prove** that the supported-threshold model has no
solution: 23 non-control points sparsely cover several uncertain dimensions,
the conversion is not validated or reef-specific, and one extreme point is
unfinished. It does, however, reject the narrow idea that simply lowering
the Allee half-saturation will create the missing crash in this Owen/1991
configuration. The next bounded mechanism test should explicitly target the
joint gap-period **internal and external larval floor**, while retaining
the density-dependent fertilisation term and checking adult survival. A
size-structured breeder mechanism is biologically motivated by Babcock et
al., but should include growth/maturation and adult mortality or grazing
consequences; size-dependent fecundity alone cannot remove surviving adults
from the manta count. Fit a spatial/zero-aware SALAD-to-manta observation
operator as a separate evidence task, and only then widen optimization or
replicate/promote a mechanism. The model default remains unchanged.

### Opt-in size-weighted breeder candidate (2026-10-08)

The [mechanism protocol](size_weighted_candidate.md) records the hypothesis,
equations, units, bounds, control, rollback and qualitative failure criteria
before promotion. The sibling `COTSMod.jl` core now has an **off-by-default**
two-adult-size production switch. The `N[3]` state and the 1991 observation
mapping still mean **all adults >=15 cm in COTS/ha**. Within that total, a
tracked >25 cm class and the residual 15–25 cm class contribute to potential
larvae with weights 1 and `exp(beta*(D_small-D_large))`, respectively. The
default representative 200/300 mm classes and published female gonad slope
`beta=0.0115 mm^-1` imply a small-adult weight of about 0.317. These are
diagnostic class representatives, not measured historical size means or a
validated gamete-to-recruit conversion. The existing 3 COTS/ha fertilisation
half-saturation acts on **total** adults. The optional growth rule moves
surviving small adults into the large class at fixed probability or at a
probability reduced by low coral cover; new adults enter the small class.
Existing survival and per-adult coral consumption remain unchanged, so this
candidate can lower future recruitment but cannot by itself remove a standing
population of adult grazers.

The [bounded screen](screen_1991_size_weighted.jl) replays archived capacity
pair 20 at `theta=3` as a strict control, then tests fixed versus
food-mediated size-class growth and two initial large-adult fractions plus a
350 mm large-class representative. It freezes the coherent but provisional
`q=0.2281847 COTS/tow per adult COTS/ha`, 1991 initial stage and coral state,
three-year-smoothed peak objective, Owen connectivity and external boundary,
and forcing seed. The [audit](audit_1991_size_weighted.py) checks size
accounting, all four reef trajectories, peak counts/timing/height, three
inter-peak troughs, component recruitment fluxes and coral guardrails. This
one-seed, one-anchor factorial is a *qualitative* screen, not a fitted size
history or production optimization. The four-channel `cots_size_log` stores
post-transition small and large adults, pre-transition effective breeders,
and annual small-to-large growth flux (all COTS/ha, except growth is
COTS/ha/year). The existing 12-channel `cots_flow_log` is unchanged.

```powershell
$env:JULIA_DEPOT_PATH = (Join-Path (Get-Location) '.julia_cots_screen') + ';' + (Join-Path $env:USERPROFILE '.julia')
julia --compiled-modules=existing --project=sandbox ..\COTSMod.jl\test\runtests.jl
julia --compiled-modules=existing --project=sandbox sandbox/calibration/test_1991_initial_state.jl
$env:OWEN_SIZE_RUN_ID = '<unused_size_run_id>'
julia --compiled-modules=existing --project=sandbox sandbox/calibration/screen_1991_size_weighted.jl
python sandbox/calibration/audit_1991_size_weighted.py sandbox/calibration/runs/<unused_size_run_id>
```

Do not promote the switch based on a single-anchor result. A passing treatment
still needs replicate seeds, class-size/growth and tow-conversion sensitivity,
coral support assessment, held-out reefs/years and regional validation. If
every arm fails, retain the negative diagnostic and investigate an adult
survival/grazing mechanism or shared internal/external recruitment floor
without changing the supported Allee threshold.

The [completed one-seed pair-20 run](runs/20261008Tsize_weighted_pair020_v2/metadata.toml)
passed the archived-control replay and completed all five arms without model
exceptions. The [independent audit](runs/20261008Tsize_weighted_pair020_v2/size_audit.json)
and [COTS/coral figure](runs/20261008Tsize_weighted_pair020_v2/cots_coral_size_weighted_pair020.png)
show that **none** passes the qualitative cycle gate. Every arm matches at
most six of seven frozen observed peaks; all have incorrect 2/2/2/1 peak
counts, late first peaks, and inter-peak adult-density floors far above the
required <=0.10 trough-to-smaller-peak ratio. The worst of the three
two-wave reef trough ratios is 0.360 for control, 0.359 with food-mediated
size growth, and 0.383 when the large representative is 350 mm. Mean frozen
loss falls from 3.130 to 2.967 in the food-mediated arm, but the modest loss
change is **not** a mechanism pass. First modeled peaks remain 2002/2002/2001
at Lizard/MacGillivray/North Direction versus observed 1997/1997/1994, and
second modeled peaks remain near 2014-15.

The mechanism acts as intended on production. Across the 2004-09 gap,
food-mediated weighting changes effective breeder density at Lizard from
1.189 to 0.603 COTS/ha, local potential production from 6.253 to 4.977
COTS/ha/year, and internal settled immigration from 1.779 to 0.912
COTS/ha/year. At MacGillivray the corresponding production drops from
15.652 to 7.466 and internal settlement from 2.341 to 1.181. External
settlement is unchanged in each paired reef (Lizard 0.226,
MacGillivray 0.197 COTS/ha/year), because this treatment does not alter the
out-of-domain boundary. These reductions do not translate into deep enough
adult troughs: existing adults survive and graze, and settlement from the
external boundary persists. Modeled coral still drops far below 9 m manta
observations, though the supports differ. Thus the size proxy is retained
**experimental and off by default**. Do not compensate by lowering the
supported Allee threshold. A next discriminating test should combine an
explicitly bounded post-peak adult survival/grazing response with size
weighting, with external supply controlled in a separate arm, and require
both peaks and coral guardrails before optimization.

### 1991 cohort-anchored crash factorial (2026-10-08)

The next [preregistered screen](crash_mechanism_protocol.md) returns to the
1991 full-cohort state that recovers Lizard's first-peak year. It crosses
the existing opt-in coral-dependent juvenile-storage switch with adult
background mortality `m3=0.30/year`, retaining the archived control, Owen
boundary, `q=0.015` observation conversion sensitivity, Allee `3 COTS/ha`,
reef mapping, hydrodynamic seed, and score. The aim is to see whether delayed
maturation and adult attrition interact to deepen the gap *without* losing
either peak; it is not an optimization or a validated cohort reconstruction.
The [runner](run_1991_initialization_test.jl) uses
`OWEN_1991_EXPERIMENT=crash_factorial`; the [independent audit](audit_1991_crash_factorial.py)
checks the archived control/state replay, unchanged external immigration,
peak/trough/coral diagnostics and component fluxes, and makes a paired plot.
Results and the promotion decision are added below only after the audit passes.

```powershell
$env:JULIA_DEPOT_PATH = (Join-Path (Get-Location) '.julia_cots_screen') + ';' + (Join-Path $env:USERPROFILE '.julia')
$env:OWEN_1991_EXPERIMENT = 'crash_factorial'
$env:OWEN_1991_RUN_ID = '<unused_crash_factorial_run_id>'
julia --compiled-modules=existing --project=sandbox sandbox/calibration/run_1991_initialization_test.jl
python sandbox/calibration/audit_1991_crash_factorial.py sandbox/calibration/runs/<unused_crash_factorial_run_id>
```

The [completed four-arm run](runs/20261008T_crash_factorial_1991_v1/metadata.toml)
and [independent audit](runs/20261008T_crash_factorial_1991_v1/crash_factorial_audit.json)
replayed the archived cohort control and all its logged fluxes exactly
(`max absolute difference = 0`) and held external immigration identical
across arms. The Julia metadata recorded Git revisions as `unknown`; a
[post-run provenance supplement](runs/20261008T_crash_factorial_1991_v1/provenance_supplement.json)
records actual revisions, dirty flags and source/manifest hashes. This is a
dirty-worktree diagnostic—both repositories contained uncommitted changes,
not corrupted model files—so revision IDs alone cannot reconstruct it. See the
[paired COTS/coral figure](runs/20261008T_crash_factorial_1991_v1/cots_coral_crash_factorial.png)
and [per-reef readout](runs/20261008T_crash_factorial_1991_v1/crash_factorial_audit.csv).

| 1991 full-cohort arm | Worst available two-wave trough/smaller-peak | Observed peaks matched / 7 | Lizard first detected peak | Lizard first height, COTS/tow |
| --- | ---: | ---: | ---: | ---: |
| Control | 0.594 | 6 | 1997 | 0.847 |
| Juvenile stage gate | 0.691 | 5 | 2007 | 0.401 |
| Adult `m3=0.30` | 0.515 | 5 | 2005 | 0.370 |
| Stage gate + `m3=0.30` | 0.474* | 4 | no first-wave peak | — |

`*` The combined-arm ratio excludes Lizard because it no longer has two
detected peaks; it is **not** a passing three-reef trough result. The stage
gate retains juveniles and later recruitment rather than removing adults:
at Lizard its 2004–12 mean adult density rises from `1.028` to `1.126
COTS/ha`, mean retained juveniles from `0` to `0.309 COTS/ha/year`, and
2012 model coral falls from `0.093` to `0.078`. Constant mortality slightly
spares coral but sacrifices the early peak; even the combined arm remains
far from the `<=0.10` trough criterion. The paired plot also shows that
whole-reef model coral generally remains below 9 m manta coral observations,
while those two supports are not directly interchangeable. **No arm is
promoted.** The next discriminating opt-in candidate is an adult-specific
lagged density burden, with density-only versus density × coral-stress arms
and independently bounded hazard parameters, while preserving the same
archived control and upstream boundary. It must not be implemented as a
calendar-wave timer or mislabelled as confirmed pathogen mortality. The
bounded implementation and its negative result are recorded below.

### Opt-in lagged adult-hazard screen (2026-10-08)

The [protocol](crash_mechanism_protocol.md) preregisters a site-level burden
`B[t+1]=rho*B[t]+(1-rho)*A[t]`, with `B[0]=0` in adult COTS/ha. A density
factor switches on only above its separate `1.5 COTS/ha` threshold; optional
coral stress multiplies it. The resulting one-year hazard removes a fraction
`1-exp(-h)` of adults that survived existing background and food mortality.
Newly matured adults are not subject to that additional hazard until the next
step. The switch `COTS_LAGGED_ADULT_HAZARD` is **false by default**, and
`COTS_HAZARD_MAX=0` also reproduces the legacy transition. All five arms hold
the supported `3 COTS/ha` Allee half-saturation, Owen boundary, 1991 cohort
seed, provisional `q=0.015`, hydrodynamic seed, area-weighted reef mapping,
and three-year-smoothed score fixed. Diagnostic strengths are `h_max=0.75`
and `1.5`; these are **not** estimated pathogen-mortality rates.

The [completed run](runs/20261008T_lagged_hazard_1991_v1/metadata.toml),
[independent audit](runs/20261008T_lagged_hazard_1991_v1/lagged_hazard_audit.json),
[reef readout](runs/20261008T_lagged_hazard_1991_v1/lagged_hazard_audit.csv),
[site-year hazard summary](runs/20261008T_lagged_hazard_1991_v1/site_hazard_summary.csv),
and [paired COTS/coral plot](runs/20261008T_lagged_hazard_1991_v1/cots_coral_lagged_hazard.png)
show that **none passes the qualitative gate**. The archived control and all
12 recruitment-flux channels replay exactly; external immigration is
identical in every arm. The site-level burden recurrence and extra-adult-death
mass bounds pass. A [provenance supplement](runs/20261008T_lagged_hazard_1991_v1/provenance_supplement.json)
records dirty-worktree/source hashes alongside the run metadata's revisions,
seeds, Manifest hash, settings and output hashes.

| Arm | Worst available two-wave trough/smaller-peak | Matched observed peaks / 7 | Lizard first detected peak |
| --- | ---: | ---: | ---: |
| Control | 0.593935 | 6 | 1997 |
| Density-only, `h_max=0.75` | 0.616370* | 5 | 2017* |
| Density-only, `h_max=1.5` | 0.600000 | 6 | 2001 |
| Density × food, `h_max=0.75` | 0.605740* | 5 | 2017* |
| Density × food, `h_max=1.5` | 0.593493 | 6 | 2001 |

`*` The weak-hazard arms lose Lizard's first detected wave; their worst ratio
excludes that reef. The strongest food-coupled arm's `0.000442` change in
worst trough ratio is negligible and comes with delayed/weaker first peaks.
The hazard acts early, then shuts off: in the strong density-only arm it is
active at about `33%` of sites in 1995 and `51%` in 2005, but only `0.5%`
in 2008 and `0.3%` in 2012. At North Direction, mean added adult deaths in
the 2004–12 gap are effectively zero despite ongoing external settlement
(`0.529 COTS/ha/year`). Coral improves slightly but remains low relative to
9 m manta observations, with different observation support. This is a
**negative structural diagnostic, not a calibrated mortality mechanism**.
Changing only the memory weight to `0.75` or `0.9` in an *offline recurrence*
on the fixed control trajectory leaves no more than `1%` of sites above the
same burden threshold by 2012; that calculation is not a model treatment.
The next gate is a separately justified post-peak persistence/recovery or
boundary mechanism, not a broad optimization of this hazard.

```powershell
$env:JULIA_DEPOT_PATH = (Join-Path (Get-Location) '.julia_cots_screen') + ';' + (Join-Path $env:USERPROFILE '.julia')
julia --compiled-modules=existing --project=sandbox ..\COTSMod.jl\test\runtests.jl
julia --compiled-modules=existing --project=sandbox sandbox/calibration/test_1991_initial_state.jl
$env:OWEN_1991_EXPERIMENT = 'lagged_hazard'
$env:OWEN_1991_RUN_ID = '<unused_lagged_hazard_run_id>'
$env:COTS_RUN_ADRIA_REVISION = (git rev-parse HEAD).Trim()
$env:COTS_RUN_COTSMOD_REVISION = (git -C ..\COTSMod.jl rev-parse HEAD).Trim()
julia --compiled-modules=existing --project=sandbox sandbox/calibration/run_1991_initialization_test.jl
python sandbox/calibration/audit_1991_lagged_hazard.py sandbox/calibration/runs/<unused_lagged_hazard_run_id>
```

### Readout-only cohort floor audit (2026-10-09)

The [exact-control replay](runs/20261009T_cohort_audit_1991_v1/metadata.toml)
adds [area-weighted annual N1/N2/adult and transition diagnostics](runs/20261009T_cohort_audit_1991_v1/cohorts.csv)
without changing COTS dynamics. Adults, coral, and all 12 recruitment-flow
channels agree with the archived 1991 full-cohort control exactly (`max absolute
difference = 0`). The [Lizard cohort plot](runs/20261009T_cohort_audit_1991_v1/cohort_floor_diagnostics.png)
shows why the adult series sawtooths after the first peak: food survival
alternates sharply during 1999–2007, affecting both N2 maturation and
established adults. It is not an adult-only mortality factor. N1 recruits are
not food-limited by that factor. By 2008–12 the mean realized adult food factor
is `0.986`; Lizard adults average `0.805 COTS/ha`, split into `0.667` surviving
adults and `0.138` newly matured. External settlement still averages `0.542
COTS/ha/year` and maintains the later cohort pipeline. The frozen Lizard
observation scale is `0.2942 COTS/tow per adult COTS/ha`, so `0.1 COTS/tow`
corresponds to `0.340 adult COTS/ha` **under that provisional mapping**. The
`q=0.015` factor used to seed 1991 adults is not the plotting scale.

[ReefMod 7.0's settings](https://github.com/ymbozec/REEFMOD.7.0_GBR/blob/main/settings/settings_COTS.m)
and [transition](https://github.com/ymbozec/REEFMOD.7.0_GBR/blob/main/functions/f_runmodel.m)
use a density-gated `4–6`-year checked outbreak duration and reset older age
classes to low background abundance while leaving younger recruits. A
different-version [regional report](https://www.barrierreef.org/uploads/CCIP-R-04-Final-Report-Regional-Modelling.pdf)
describes a `2–5`-year disease window and separate low-preferred-coral reset.
These justify an opt-in *phenomenological benchmark*, not a claim that disease
caused the observed Lizard crash or that ReefMod parameters transfer directly
between `COTS/400 m²`, `COTS/ha`, and manta-tow units. The bounded next-test
design and failure gates are in the [crash protocol](crash_mechanism_protocol.md);
no timer was implemented or promoted in this audit.

To reproduce the cohort readout with a new run ID:

```powershell
$env:JULIA_DEPOT_PATH = (Join-Path (Get-Location) '.julia_cots_screen') + ';' + (Join-Path $env:USERPROFILE '.julia')
$env:OWEN_1991_EXPERIMENT = 'cohort_audit'
$env:OWEN_1991_RUN_ID = '<unused_cohort_audit_run_id>'
julia --compiled-modules=no --pkgimages=no --project=sandbox sandbox/calibration/run_1991_initialization_test.jl
python sandbox/calibration/plot_1991_cohort_audit.py sandbox/calibration/runs/<unused_cohort_audit_run_id>
```

### Preferred-prey starvation × adult senescence screen (2026-10-09)

The [preregistered four-arm screen](crash_mechanism_protocol.md) crossed two
opt-in COTSMod switches on the 1991 full-cohort anchor:
- `preferred_prey_starvation`: food survival from remembered fast cover, threshold `0.05`, memory `0.5` yr;
- `adult_senescence`: adult ages 2–5 at the archived `m3`, and a 6+ class at `0.8` annual mortality.

Everything else was held fixed: seed, `q=0.015` preload, 3 COTS/ha Allee threshold, external outbreak production, connectivity, fecundity and consumption, frozen observation scales and score. `per_capita_consumption` is implemented but was **not run**, pending the units decision below. See [metadata](runs/20261009T_structural_crash_1991_v1/metadata.toml), [audit](runs/20261009T_structural_crash_1991_v1/structural_crash_audit.json), [per-reef table](runs/20261009T_structural_crash_1991_v1/structural_crash_audit.csv) and [COTS/coral plot](runs/20261009T_structural_crash_1991_v1/cots_coral_structural_crash.png).

Integrity checks:
- The control replays the archived trajectories and all 12 flux channels exactly (max delta `0`).
- External immigration is identical in every arm.
- Mechanism adults equal trajectory adults exactly.
- Switched-off diagnostics are inert.

Tests: COTSMod `205/205`; the 1991 adapter test file passes, including the new 20-assertion adapter set; parse tests `36/36`.

| Arm | Two-wave reefs with 2 peaks | Worst trough ratio | Matched peaks | Lizard first peak | Min first-peak height vs control | Min 2012 coral vs control |
| --- | --- | --- | --- | --- | --- | --- |
| Control | 3/3 | `0.594` | 6/7 | 1997 | 1.00 | 1.00 |
| Preferred starvation | 2/3 | `0.723` | 6/7 | 1995 | 0.14 | 2.94 |
| Senescence | 3/3 | `0.473` | 6/7 | 1997 | 0.83 | 1.09 |
| Both | 2/3 | `0.710` | 6/7 | 1995 | 0.14 | 3.00 |

**No arm passes the gate; nothing is promoted.**

**Preferred-prey starvation.** The new food channel exposes a structural fact the legacy rule hid: in the control, Lizard **fast coral is ~0 from 1992 to 2024**. The archived `a_F=1.18` with linear (`h=0`) grazing removes a fraction `a_F·A ≥ 1` of fast coral per year at the 1991–92 densities. Model COTS have therefore lived on massives for the whole run, which total-cover starvation could not see.

Under preferred-prey starvation, COTS starve immediately (Lizard food survival `0.20` in 1994). Peaks fall to 14% of control. Total coral recovers towards, and at MacGillivray and North Direction close to, the LTMP 9 m values.

The residual population (about `0.3 COTS/ha`) still removes about `a_F·0.3 ≈ 35%` of fast coral each year. Fast coral is held at `0.02–0.035`, just under the threshold, so the consumer–resource lock-in has moved from total to preferred coral rather than becoming boom–bust.

This is **not a fair test of H1**. The archived `a_F` was fitted jointly with total-cover starvation, so H1 cannot be judged until consumption is restructured.

**Adult senescence.** This arm behaves as specified. It lowers the worst trough ratio (`0.594→0.473`; Lizard `0.399→0.326`), keeps Lizard's 1997 first peak (height −17%), and slightly raises 2012 coral.

It is insufficient because the 2004–12 adults are mostly *young*. Lizard maturation averaged `0.32 COTS/ha/yr` in the gap, fed by internal immigration of `1.0–1.8 COTS/ha/yr` in 2000–04 from asynchronously outbreaking domain reefs, plus external supply. Senescent deaths average only `~0.05 COTS/ha/yr` there. Ageing removes the outbreak cohort, but continuous settlement replaces it.

**Decision and next gate.** Keep both switches experimental and off by default. The binding constraints are now:
1. Linear grazing with an `a_F` that removes all preferred coral at about one model adult/ha.
2. The unresolved density unit: the frozen `0.2942 COTS/tow per model COTS/ha`, versus `~0.015` in CoCoNet (Moran & De'ath 1992) and implied by ReefMod.
3. Settlement supply in the gap from asynchronously outbreaking domain reefs.

Testing `per_capita_consumption` (or a saturating `h>0`) with preferred-prey starvation requires an explicit, documented decision to fix the observation scale or otherwise justify the density unit. That is a change to the frozen observation model, not a mechanism tweak.

Reproduce:

```powershell
$env:OWEN_1991_EXPERIMENT = 'structural_crash'
$env:OWEN_1991_RUN_ID = '<unused_structural_run_id>'
julia --project=sandbox sandbox/calibration/run_1991_initialization_test.jl
python sandbox/calibration/audit_1991_structural_crash.py sandbox/calibration/runs/<unused_structural_run_id>
```
