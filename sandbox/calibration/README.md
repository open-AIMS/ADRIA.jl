# COTSMod Calibration Workflow

This folder is the working seed for the standalone COTSMod calibration study repository. The goal is to calibrate COTS outbreak dynamics against Lizard Island reef observations while keeping the ecological model in `COTSMod.jl` and the ecosystem orchestration in `ADRIA.jl`.

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

## Components

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

The candidate log files record full cycle metrics:
- `bbo_evaluated_candidates.csv`
- `bbo_best_summary.csv`

The by-reef evaluation includes:
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

The overall objective sorts candidates by:
```text
loss = mean_cycle_loss + 0.25 * legacy_loss
```
where `legacy_loss` maintains backward compatibility while making cycle quality the dominant objective.


## Validation Tests

Run the metric tests with:

```powershell
julia --project=sandbox sandbox\calibration\test_cots_cycle_metrics.jl
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
$env:BBO_MAX_TIME = '0.0'     # Time limit in seconds (0.0 = unlimited)
julia --project=sandbox sandbox\calibration\calibrate_cots_blackbox.jl
```

### Candidate Evaluation Logging

As specified, every evaluated candidate is saved into a detailed CSV file:

```text
sandbox/data/bbo_evaluated_candidates.csv
sandbox/data/bbo_best_summary.csv
```

Each logged evaluation contains the exact parameter candidate vector along with the overall loss and full metric breakdown (`cycle_loss`, `legacy_loss`, `mean_rmse`, `mean_abs_percent_bias`, `mean_peak_count_penalty`, `mean_peak_timing_penalty`, `mean_period_penalty`, `mean_amplitude_penalty`, `mean_flatline_penalty`, `mean_lag_correlation_penalty`, `mean_best_lag_years`). This allows full diagnosis of whether BlackBoxOptim is improving true cycle mechanics rather than gaming individual metric components.

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
$env:COTS_N_STOCHASTIC_SCENS = '1'
julia --project=sandbox sandbox\calibration\simulate_best_calibration.jl

# Step 3: Render publication-quality plots
python sandbox\calibration\render_calibration_figures.py
```