# Lizard Island COTS Calibration Sandbox

> [!WARNING]
> The calibration sections in this document describe the archived LHS-era
> workflow. Use `sandbox/calibration/README.md` for the active BlackBoxOptim
> workflow, current status, and development gates. Use
> `sandbox/domain_building/README.md` for the versioned V2 domain build and
> connectivity contracts.

## Overview

This sandbox contains the scripts, data, and outputs for calibrating the ADRIA
Crown-of-Thorns Starfish (COTS) predation submodel against historical AIMS manta
tow survey data for the **Lizard Island cluster** (1985–2024).

The goal is to reproduce two empirical COTS outbreak peaks:
- **Peak 1 ≈ 1996–1998** — the first major outbreak at Lizard Island
- **Peak 2 ≈ 2012–2014** — the second outbreak, approximately 15 years later

---

## Current domain versions

The reproducible control remains `Lizard_Historical_v0.1`. The additive
`Lizard_Historical_v0.2` package contains the final V2 subreef polygons,
`indexV2`-labelled interim site connectivity, corrected historical DHW, a legacy
COTS-connectivity control, and an explicitly opt-in Owen 2018-2023 COTS option.
See `domain_building/README.md` before comparing either domain.

## Directory Structure

```
sandbox/
├── README.md                     ← This file
├── Project.toml / Manifest.toml  ← Julia sandbox environment
│
├── calibration/                  ← Parameter calibration scripts
│   ├── calibrate_lizard_cots.jl      BlackBoxOptim single-best optimisation
│   ├── formal_calibration_sweep.jl   Latin Hypercube Sampling (250 scenarios)
│   ├── analyze_calibration.jl        Extract & re-run top-5 ensemble
│   └── save_best_calibration.jl      Run the best param set and export CSVs
│
├── domain_building/              ← Domain construction scripts
│   ├── build_lizard_domain.jl        Build Lizard domain from RME data
│   └── build_historical_dhw.jl       Build historical DHW NetCDF (1985–2024)
│
├── diagnostics/                  ← Stand-alone model experiments
│   └── toy_cots_cycles.jl            Toy model comparing 5 recruitment variants
│
├── plotting/                     ← Visualisation scripts and outputs
│   ├── plot_calibration.py           Single best-fit trajectory vs empirical
│   ├── plot_ensemble.py              Top-5 ensemble overlay vs empirical
│   ├── plot_disturbances.py          DHW disturbance timeline per reef
│   ├── calibration_plot.png          Latest single-best plot
│   ├── ensemble_calibration_plot.png Latest ensemble plot
│   └── disturbance_plot.png          Latest disturbance plot
│
├── data/                         ← Input data and intermediate outputs
│   ├── Lizard_Historical_v0.1/       Lizard domain package (sites, DHW, conn)
│   ├── reef_cots.csv                 AIMS COTS manta tow data (all GBR reefs)
│   ├── reef_manta.csv                AIMS hard coral cover manta tow data
│   ├── dhw_historical.csv            Historical DHW per site (1985–2024)
│   ├── formal_calibration_ensemble.csv  LHS sweep results (250 runs)
│   ├── calibration_results_sim.csv   Simulated trajectories (single-best)
│   ├── calibration_results_emp.csv   Filtered empirical COTS data
│   ├── calibration_results_emp_coral.csv  Filtered empirical coral data
│   └── top5_trajectories.csv        Top-5 ensemble trajectories
│
└── archive/                      ← Obsolete/one-off scripts (kept for reference)
```

---

## COTS Population Model (`COTSMod.jl`)

The ecological transition lives in the sibling `COTSMod.jl` package
(`src/COTSMod.jl`). `ADRIA/src/ecosystem/cots.jl` is a thin adapter that maps
ADRIA parameters and `COTS_*` environment switches onto `COTSMod.COTSParams`,
converts coral tensors, and logs diagnostics. Densities are COTS ha^-1; coral
prey cover is a fraction of each site's habitable area. Defaults below are the
`COTSParams` defaults; calibrated values are recorded in each run's metadata.

### Architecture (legacy default path)

The model is a **stage-structured predator-prey system** with three classes
per site, stepped annually:

| Class | Description | Annual survival |
|-------|-------------|-----------------|
| `N[1]` | Settled recruits (age 0) | `1-m1` (default `m1=0.4`); not food-limited |
| `N[2]` | Juveniles (age 1) | `(1-m2)·f` (default `m2=0.2`) |
| `N[3]` | Adults, age 2+ plus group | `(1-m3)·f` (default `m3=0.1`) |

`f` is food survival. By default `f=1` while **total** prey cover `F+S`
exceeds `0.15·C_max`, then falls cubically to `1-p_tilde`. The exponent 3 is
hard-coded; the sampled `eta_starve` factor is currently **unused**, as are the
Beverton-Holt `a`, `b` and `fecundity_gate` factors.

1. **Ricker recruitment with Allee fertilisation** —
   `larvae = a_ricker · condition² · A · exp(-b_ricker·A) · A²/(allee² + A²)`.
   With `COTS_SEPARATED_RECRUITMENT=true` production is dispersed through the
   reef-scale connectivity before settlement.
2. **Maternal body condition** — an exponential moving average of total cover
   with timescale `tau_condition` (default 5 years). It scales fecundity only.
3. **Consumption** — `Cons_F = A·a_F·F^eta_F / (1 + h(a_F F^eta_F + a_S S^eta_S))`
   using post-transition adults. The default `h=0` (used in all calibrated
   runs to date) makes this **linear (Type I)**: the fraction of fast coral
   removed each year is `a_F·A`, independent of how much remains.
4. **Coral-gated background immigration** — `IMM` scaled by a cover gate.

### Opt-in experimental switches (all default off)

Each switch is set through `COTS_*` environment variables in `scenario.jl`.
None is promoted; see `calibration/MODEL_LOG.md` and
`calibration/crash_mechanism_protocol.md` for evidence and decisions.

| Switch | Mechanism |
|--------|-----------|
| `COTS_JUVENILE_STORAGE` | Cover-dependent maturation; unmatured juveniles are retained |
| `COTS_HABITAT_MEDIATION` | Settlement gated by open habitat |
| `COTS_LOW_COVER_MORTALITY` | Extra mortality below a total-cover threshold |
| `COTS_SIZE_WEIGHTED_FECUNDITY` | Two adult size classes with size-weighted fecundity |
| `COTS_LAGGED_ADULT_HAZARD` | Density/food-triggered adult hazard with lagged burden |
| `COTS_PREFERRED_PREY_STARVATION` | Food survival driven by remembered **fast** cover (threshold `COTS_STARVATION_PREFERRED_THRESHOLD`, memory `COTS_STARVATION_MEMORY_YEARS`) |
| `COTS_ADULT_SENESCENCE` | Adult age classes 2..`COTS_SENESCENCE_AGE`; the plus group dies at `COTS_SENESCENCE_MORTALITY` |
| `COTS_PER_CAPITA_CONSUMPTION` | Fixed m² of coral per adult per year (`COTS_CONSUMPTION_M2_PER_ADULT_YEAR`), preferred prey first. Requires a resolved density unit; not yet run |

### Prey Categories

| Category | Functional Groups | Consumption Rate |
|----------|-------------------|------------------|
| Fast (preferred) | Tabular Acropora, Corymbose Acropora, Corymbose non-Acropora | `a_F` |
| Slow | Small massives, Large massives | `a_S` |

### Spatial Components

- **Larval dispersal** — the legacy `mean` mode uses the domain's COTS (or
  coral-proxy) site matrix. The `cycle`/`sample` modes use annual reef-scale
  COTS matrices (e.g. Owen 2018–2023, source rows → sink columns), with
  optional source larval survival and an evidence-derived external outbreak
  boundary. See `domain_building/README.md`.
- **Initialization** — from the domain's COTS probability surface scaled by
  `COTS_INITIAL_MULTIPLIER`, or from an explicit site-by-stage CSV via
  `COTS_INITIAL_STATE_CSV`.

---

## Calibration Parameters (archived LHS-era workflow)

### Primary Calibration Targets (swept in LHS)

| Parameter | Symbol | Bounds | Best Candidate | Role |
|-----------|--------|--------|----------------|------|
| Starvation threshold | `a_F` | [0.1, 2.0] | 1.47 | Fast coral consumption rate; controls outbreak severity |
| Slow coral consumption | `a_S` | [0.01, 0.9] | 0.071 | Slow coral consumption; affects prey-switching timing |
| Background immigration | `IMM` | [0.0, 0.1] | 0.038 | External larval input; controls re-seeding between cycles |
| Initial seed multiplier | `seed_mult` | [0.5, 3.0] | 1.63 | Scales starting COTS density; controls Peak 1 timing |

### Loss Function

The calibration uses a hybrid loss combining point-wise SSE with a
**Dual-Peak Phase-Shift Penalty**:

```
loss = Σ_reef [ Σ_obs (sim_norm - emp_norm)²
              + |sim_peak1_year - emp_peak1_year| × 2.0
              + |sim_peak2_year - emp_peak2_year| × 2.0 ] / n_reefs
```

- **Peak 1** is identified as `argmax(sim_cots_norm[1:40])` and compared
  against the empirical peak (≈1997 for most Lizard reefs).
- **Peak 2** is identified as `argmax(sim_cots_norm[21:40])` and compared
  against empirical data from years > 2005 (≈2013).

### Current Calibration Status

| Metric | Status |
|--------|--------|
| Peak 1 timing (≈1997) | ✓ Well aligned |
| Peak 2 timing (≈2013) | ○ Lagging 6–10 years |
| Peak 1 amplitude | ✓ Reasonable |
| Peak 2 amplitude | ○ Under investigation |
| Inter-cycle period | ○ Currently ~20–22 yr (target: ~15 yr) |

---

## Known Issues & Limitations

1. **Second cycle lag** — The simulated second outbreak consistently peaks
   6–10 years later than observed. The model's natural oscillation period
   is ~20 years rather than the empirical ~15 years. Likely causes:
   - Coral recovery too slow (logistic growth in CoralBlox)
   - Starvation kills COTS too aggressively, extending the bust phase
   - No external larval pulse from upstream reefs (Cairns initiation box)

2. **Coral connectivity used as COTS proxy** — The dispersal matrix was
   calibrated for coral larvae. COTS produce ~60M eggs per female vs
   ~10K for corals. The `immigration_scalar` partially compensates but
   is not biologically rigorous.

3. **ENV-based initial density injection** — The `COTS_INITIAL_MULTIPLIER`
   environment variable is a calibration workaround. This should eventually
   be replaced with a proper initial condition parameter in the ADRIA
   parameter table.

4. **Single DHW scenario** — The historical domain uses duplicated DHW
   columns (same values in both "scenarios"). No stochastic DHW sampling.

5. **No cyclone forcing** — Cyclone tracks are zeroed in the historical
   domain. Real cyclone damage (e.g., Cyclone Ita 2014) likely influenced
   the COTS-coral dynamics.

---

## Empirical Data Sources

| Dataset | Source | Coverage |
|---------|--------|----------|
| `reef_cots.csv` | AIMS LTMP manta tow COTS counts | 1985–2024, all GBR reefs |
| `reef_manta.csv` | AIMS LTMP manta tow hard coral cover | 1985–2024, all GBR reefs |
| `dhw_historical.csv` | eReefs/CoralWatch historical DHW | 1985–2024, Lizard cluster sites |
| `COTS_prob_*.tif` | COTS habitat suitability raster | Spatial initialization |

### Reef Name Mapping (Simulation → AIMS)

| Simulation Name | AIMS Survey Name |
|----------------|------------------|
| Lizard Island Reef | Lizard Isles |
| MacGillivray Reef | Macgillivray Reef |
| North Direction Reef | North Direction Island |
| Eyrie Reef | Eyrie Reef |

---

## How to Run

The workflow below describes the older V1 calibration examples. Do **not**
use it for V2 production optimization: the V2 optimizer is gated. The current
V1/V2/Owen area-weighted, three-seed comparison and updated COTS peak plot are
documented in `sandbox/calibration/README.md` under "Paired domain gate, Owen
coverage, and observation smoothing". Source-build and missing-survival audits
are documented in `sandbox/domain_building/README.md`.

### Prerequisites
- Julia 1.10+ with ADRIA.jl activated (`Pkg.activate(".")` from repo root)
- Python 3.10+ with `pandas` and `matplotlib` for plotting
- `BlackBoxOptim.jl` and `LatinHypercubeSampling.jl` for calibration

### Workflow

```bash
# 1. Rebuild or validate the versioned V2 domain using the commands in
#    sandbox/domain_building/README.md

# 2. Run BBO optimisation (single best candidate)
julia sandbox/calibration/calibrate_lizard_cots.jl

# 3. Export best candidate trajectories for plotting
julia sandbox/calibration/save_best_calibration.jl

# 4. Plot single-best calibration
python sandbox/plotting/plot_calibration.py

# 5. Run 250-scenario LHS ensemble sweep
julia sandbox/calibration/formal_calibration_sweep.jl

# 6. Analyse & extract top-5 from ensemble
julia sandbox/calibration/analyze_calibration.jl

# 7. Plot ensemble overlay
python sandbox/plotting/plot_ensemble.py
```
