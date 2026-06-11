# Lizard Island COTS Calibration Sandbox

## Overview

This sandbox contains the scripts, data, and outputs for calibrating the ADRIA
Crown-of-Thorns Starfish (COTS) predation submodel against historical AIMS manta
tow survey data for the **Lizard Island cluster** (1985–2024).

The goal is to reproduce two empirical COTS outbreak peaks:
- **Peak 1 ≈ 1996–1998** — the first major outbreak at Lizard Island
- **Peak 2 ≈ 2012–2014** — the second outbreak, approximately 15 years later

---

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

## COTS Population Model (`src/ecosystem/cots.jl`)

### Architecture

The model is a **stage-structured predator-prey system** with three age classes
tracked per spatial location:

| Age Class | Description | Key Mortality |
|-----------|-------------|---------------|
| Age 0 | Recruits (settled larvae) | `m1` (default 0.4) |
| Age 1 | Juveniles | `m2` (default 0.2) |
| Age 2+ | Adults (predatory stage) | `m3` (default 0.08) + starvation |

### Key Biological Mechanisms

1. **Ricker Recruitment** — Allows overcompensation (massive larval pulses)
   gated by maternal body condition. At high adult density, density-dependent
   effects reduce per-capita recruitment, but the *total* larval output can
   still spike dramatically:
   ```
   R = a_ricker × body_condition² × N_adults × exp(-b_ricker × N_adults)
   ```

2. **Maternal Body Condition** — An exponential moving average of food
   availability with memory timescale `tau_condition` (default 3 years).
   This creates a critical *lag*: good conditions now → larval explosion
   2–3 years later when recruits mature to adults.

3. **Threshold Starvation** — Adult mortality stays near-zero while coral
   prey exceeds 15% of `C_max`. Below that threshold, survival crashes
   steeply (cubic drop-off). This is biologically accurate: COTS can
   survive for extended periods with moderate food, then starve rapidly.

4. **Holling Type II/III Functional Response** — Consumption saturates at
   high prey density via handling time `h`. The generalized form supports
   both Type II (linear numerator) and Type III (sigmoidal) via exponents
   `eta_F` and `eta_S`:
   ```
   Cons_F = N_adults × (a_F × F^eta_F) / (1 + h × (a_F × F^eta_F + a_S × S^eta_S))
   ```

5. **Coral-Gated Immigration** — External larval immigration is modulated
   by local food availability. Immigration shuts off when coral cover is
   too low to support incoming settlers.

6. **Allee Effect** — Fertilisation success crashes when the adult
   population is too sparse, preventing runaway recruitment from near-zero
   populations.

### Prey Categories

| Category | Functional Groups | Consumption Rate |
|----------|-------------------|------------------|
| Fast (preferred) | Tabular Acropora, Corymbose Acropora, Corymbose non-Acropora | `a_F` |
| Slow | Small massives, Large massives | `a_S` |

### Spatial Components

- **Larval Dispersal** (`disperse_cots_larvae!`) — Uses the coral
  connectivity matrix as a proxy for COTS larval transport. An
  `immigration_scalar` bridges the orders-of-magnitude difference in
  larval production between corals and COTS.
- **Spatial Initialization** (`init_cots_from_spatial`) — Seeded from a
  COTS probability raster (habitat suitability), scaled by an
  `ENV["COTS_INITIAL_MULTIPLIER"]` for calibration tuning.

---

## Calibration Parameters

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

### Prerequisites
- Julia 1.10+ with ADRIA.jl activated (`Pkg.activate(".")` from repo root)
- Python 3.10+ with `pandas` and `matplotlib` for plotting
- `BlackBoxOptim.jl` and `LatinHypercubeSampling.jl` for calibration

### Workflow

```bash
# 1. (Optional) Rebuild domain from RME data
julia sandbox/domain_building/build_lizard_domain.jl
julia sandbox/domain_building/build_historical_dhw.jl

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
