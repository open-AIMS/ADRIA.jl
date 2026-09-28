# COTSMod Calibration Workflow

This folder is the working seed for the standalone COTSMod calibration study repository. The goal is to calibrate COTS outbreak dynamics against Lizard Island reef observations while keeping the ecological model in `COTSMod.jl` and the ecosystem orchestration in `ADRIA.jl`.

> [!IMPORTANT]
> This file is the authoritative status and execution guide for the active COTS
> calibration workflow. Parameter values and commands in `AGENT_HANDOFF.md`,
> `sandbox/README.md`, and archived scripts describe earlier experiments unless
> they are explicitly restated here.

## Current Status

- COTSMod is integrated into ADRIA through the compatibility adapter in
  `ADRIA/src/ecosystem/cots.jl`.
- Package, adapter, and synthetic metric tests pass.
- The focused, expanded, pulse-enabled, and COTS-connectivity bounded pilots
  all failed the biological promotion criteria. None is a calibrated result.
- Every expanded candidate produced at most one detected peak per reef. The
  post-2005 temporal holdout matched zero of the four observed later peaks.
- A dedicated Lizard COTS matrix is now built reproducibly from six ReefMod
  spawning seasons and selected by the scenario runner, with the coral matrix
  retained as an explicit control/fallback.
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
2. **`EXPANDED`**: 15 demographic, functional response, mortality, and dispersal parameters:
   `a_F`, `a_S`, `IMM`, `seed_mult`, `a_ricker`, `b_ricker`, `m1`, `m2`, `m3`, `p_tilde`, `C_max`, `tau_condition`, `allee_threshold`, `imm_threshold`, `eta_imm`.
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
The next bounded mechanism experiment should:

1. cross a much finer low Allee-threshold range with `a_ricker` and `IMM`;
2. log local fecundity, background immigration, and dispersed recruits
   separately rather than only total age-1 animals;
3. test an externally specified absolute larval-supply pulse distributed over
   biologically justified source reefs, rather than scaling a pulse only to the
   model's previous internal supply; and
4. reject the mechanism unless it produces a second peak before any optimizer
   expansion.

### Connectivity and validation utilities

```powershell
# Build and verify the site-level COTS matrix (large data outputs are ignored)
julia --project=sandbox sandbox\domain_building\build_lizard_cots_connectivity.jl
julia --project=sandbox sandbox\domain_building\validate_lizard_cots_connectivity.jl

# Compare same-seed coral and COTS-connectivity pilots
julia --project=sandbox sandbox\calibration\analyze_connectivity_factorial.jl

# Freeze the pre-2006 observation scale and score 2006 onward
$env:BBO_RUN_ID = 'expanded_cotsconn_pilot_seed20260930'
julia --project=sandbox sandbox\calibration\validate_temporal_holdout.jl

# Diagnose the post-outbreak age structure and Allee multiplier
julia --project=sandbox sandbox\calibration\diagnose_recurrence.jl
```

Promotion criteria, the locked observed-peak contract, expanded-search order,
sensitivity outputs, connectivity experiment, and validation design are defined
in `CALIBRATION_PROTOCOL.md`.
