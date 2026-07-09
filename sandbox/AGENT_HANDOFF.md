# ADRIA COTS Submodel & Lizard Island Calibration: Agent Handoff Document

> [!IMPORTANT]
> **To Future AI Agents (Codex, Antigravity, etc.) and Developers:**
> This document (`AGENT_HANDOFF.md`) is the authoritative **Point of Truth** for the Crown-of-Thorns Starfish (COTS) submodel calibration and integration within `ADRIA.jl`. Please read this entire document carefully before making code edits or running optimization sweeps.

---

## 1. Executive Summary & Objective
Our objective is to calibrate and refine the ADRIA COTS submodel on the **Lizard Island Cluster domain** (`Lizard_Historical_v0.1`, covering 1985–2024) to accurately reproduce historical ~15-year predator-prey boom-bust cycles and empirical observations (`reef_cots.csv` and `reef_manta.csv`).

We have successfully resolved previous lag and flatlining issues by combining:
1. A **4-parameter Latin Hypercube Sweep (LHS)** of core ecological parameters (`a_F`, `a_S`, `IMM`, `seed_mult`).
2. The **Cairns Initiation Box** external larval pulse boundary condition, which injects incoming larvae from upstream northern GBR outbreaks every 15 years starting at `t = 1` (1985).

---

## 2. Directory Structure & Key Files (`sandbox/`)
To maintain clean separation and avoid ad-hoc script sprawl, all COTS calibration and visualization workflows live inside organized subdirectories under `sandbox/`:

```
sandbox/
├── calibration/
│   ├── formal_calibration_sweep.jl    # Sequential/Parallel LHS sweep generating parameter ensembles
│   ├── analyze_calibration.jl         # Post-sweep analysis extracting top-N parameter candidates
│   ├── simulate_best_stochastic.jl    # Runs N stochastic simulations for a single fixed parameter candidate
│   └── save_best_calibration.jl       # Utility to run & export CSVs for single parameter candidates
├── domain_building/
│   ├── build_lizard_domain.jl         # Builds Lizard_Historical_v0.1 domain from RME spatial/connectivity data
│   └── build_historical_dhw.jl        # Prepares NetCDF DHW historical cubes (dhw_RCPhistorical.nc / dhw_RCP45.nc)
├── plotting/
│   ├── plot_ensemble.py               # Python Matplotlib script for multi-candidate or multi-run plotting
│   ├── plot_best_stochastic.py        # Python script plotting median + 10th-90th %ile uncertainty bands
│   ├── plot_calibration.py            # Plots single trajectory validation curves vs empirical data
│   └── plot_disturbances.py           # Plots DHW thermal stress and cyclone disturbance timelines
├── data/
│   ├── Lizard_Historical_v0.1/        # The compiled historical domain data package (spatial, connectivity, DHWs)
│   ├── reef_cots.csv                  # Empirical COTS manta/tow survey observations (1985–2024)
│   ├── reef_manta.csv                 # Empirical coral cover survey observations (1985–2024)
│   ├── formal_calibration_ensemble.csv# Full table of LHS sweep parameter candidates and their loss scores
│   ├── top5_trajectories.csv          # Simulated COTS & coral trajectories for the top 5 candidates
│   └── best_stochastic_trajectories.csv # Simulated trajectories from N stochastic runs of the single best candidate
└── AGENT_HANDOFF.md                   # Reference copy of this handoff document
```

---

## 3. Core Architecture & Integration Points (`src/`)

### A. Relative Cover Translation Layer (`src/scenario.jl` & `src/ecosystem/cots.jl`)
* **How COTS Predation Operates:** COTS predation operates on **relative coral cover (`C_cover_t`)**, not absolute cover directly. In `scenario.jl`, environmental stressors (bleaching/DHW, cyclones) are applied first to `C_cover_t`, and `cots_mortality!(C_cover_t, cots_models, cots_prey_map)` runs *afterward*. This ensures COTS correctly sense and respond to food-limited environments following bleaching or storms.
* **Initial Density Injection:** `cots.jl` reads `ENV["COTS_INITIAL_MULTIPLIER"]` during model initialization to scale the baseline density (`N[1]`) at `t = 1` across all sites.

### B. The Cairns Initiation Box (`src/scenario.jl`)
To model external larval ingress riding the East Australian Current southward from the Cairns initiation zone every ~15 years, `scenario.jl` checks `ENV` variables during every time step:
```julia
# Mimic incoming larvae/recruitment from upstream outbreak (Cairns Initiation Box)
_cots_pulse_enabled = cots_enabled && get(ENV, "COTS_EXTERNAL_PULSE", "false") == "true"
if _cots_pulse_enabled && !isnothing(_cots_seed_locs)
    _pulse_period = parse(Int, get(ENV, "COTS_PULSE_PERIOD", "15"))
    _pulse_offset = parse(Int, get(ENV, "COTS_PULSE_OFFSET", "1"))
    _pulse_val = parse(Float64, get(ENV, "COTS_PULSE_VAL", "2.0"))
    
    if mod(tstep - _pulse_offset, _pulse_period) == 0
        inject_upstream_pulse!(cots_models, _cots_seed_locs; pulse_val=_pulse_val)
    end
end
```
* **Required ENV Setup for Phase 9b & Best Candidate:**
  ```julia
  ENV["COTS_EXTERNAL_PULSE"] = "true"
  ENV["COTS_PULSE_PERIOD"] = "15"
  ENV["COTS_PULSE_OFFSET"] = "1"
  ENV["COTS_PULSE_VAL"] = "1.5"
  ```

---

## 4. Current Best-Fit Parameterization (Phase 9b Results)

Our lowest loss achieved to date is **`Loss = 23.0897`** (Candidate **`run_id = 2`** from the 250-scenario LHS sweep).

| Parameter | Symbol | Best Candidate (`run_id = 2`) Value | LHS Search Bounds | Description |
| :--- | :---: | :---: | :---: | :--- |
| **Fast-Coral Starvation Threshold** | `a_F` | **`1.34378`** | `[0.1, 2.0]` | Critical consumption rate / threshold for Acropora / fast growers |
| **Slow-Coral Starvation Threshold** | `a_S` | **`0.0957831`** | `[0.01, 0.9]` | Consumption rate / threshold for massive / slow growers |
| **Background Immigration Rate** | `IMM` | **`0.0947791`** | `[0.0, 0.1]` | Baseline low-level background immigration across reefs |
| **Initial Seed Multiplier** | `seed_mult` | **`2.23695`** | `[0.5, 3.0]` | Multiplier on starting COTS population density via `COTS_INITIAL_MULTIPLIER` |
| **Cairns Pulse Amplitude** | `pulse_val` | **`1.5`** | Fixed / ENV | Added to `Age-0` (`N[1]`) at `t=1, 16, 31` (1985, 2000, 2015) |

> [!TIP]
> **Why `run_id = 2` works so well:** The relatively high `a_F` (`1.34`) combined with low `a_S` (`0.096`) allows COTS to rapidly consume fast-growing corals during an outbreak and then experience sharp, deep starvation collapses (**busts**) down to ~7% normalized density (`2000–2008`). The external Cairns pulse arriving around `t=31` (`2015`) perfectly ignites the second historical outbreak wave, peaking around `2018–2019`.

---

## 5. Known Gotchas & Troubleshooting Guide

### A. `switch_RCPs!` on `LizardDomain`
When running `ADRIA.run_scenarios(dom, scens, "45")` multi-threaded, internal ADRIA logic calls `switch_RCPs!(dom, rcp)`. We have explicitly added the method to `src/ExtInterface/Lizard/LizardDomain.jl`:
```julia
function ADRIA.switch_RCPs!(d::LizardDomain, RCP::String)::LizardDomain
    return d
end
```
*Never remove this method from `LizardDomain.jl` unless `LizardDomain` is refactored to inherit dynamic multi-RCP paths.*

### B. Sampling Counterfactual Stochastic Scenarios (`guided = 0`)
If you call `scens = ADRIA.sample(dom, N)` to generate stochastic runs across environmental/demographic variables, `sample()` will by default perturb **Intervention / Seeding** parameters (`N_seed_TA`, `N_seed_CA`, `guided`, `fogging`, etc.).
If left active, those seeded corals will flood the domain and cause:
`Warning: Number of seeding devices exceeds available space. Excluding devices...`
* **Fix:** Always zero out intervention columns when testing counterfactual COTS dynamics:
  ```julia
  scens.guided .= 0
  for col in names(scens)
      if startswith(col, "N_seed_") || col in ["fogging", "SRM", "a_adapt"]
          scens[!, col] .= 0.0
      end
  end
  ```

### C. Serial vs Parallel Execution & Memory
* **`run_scenario(dom, p)`** runs a single scenario sequentially in ~10–15 seconds. For small ensembles (`N <= 20`), a simple `for` loop in `simulate_best_stochastic.jl` (`20 * 12s ≈ 4 minutes`) is completely robust and avoids Zarr batch I/O overhead.
* **`run_scenarios(dom, scens, "45")`** runs multi-threaded and uses Zarr chunked execution. Use this when sampling `N >= 64` scenarios.

---

## 6. Next Phases & Action Roadmap

### Phase 9a/c: COTS-Specific Larval Dispersal Matrix
Currently, COTS larval dispersal (`disperse_cots_larvae!`) uses coral connectivity as a proxy.
* **Action:** Source or generate an empirical hydrodynamic COTS connectivity matrix for the Lizard Island cluster (`Lizard_COTS_Connectivity.csv`). Update `build_lizard_domain.jl` and `scenario.jl` to pass this COTS-specific connectivity into `disperse_cots_larvae!` separately from coral connectivity.

### Phase 10: Formal Global Sensitivity Analysis (Sobol / PAWN)
* **Action:** Execute a formal Sobol or PAWN sensitivity sweep using `ADRIA.sensitivity` across the `[a_F, a_S, IMM, seed_mult, pulse_val]` parameter space around the `run_id = 2` regime to quantify the exact first-order and total-order sensitivity indices driving cycle period and trough depth.

### Phase 11+: Regional & Full-Scale GBR Cross-Validation
* **Action:** Validate the calibrated COTS parameter set against empirical AIMS tow datasets across broader Moore/Cairns regional domains (`RMEDomain`) to confirm generalizability outside the Lizard Island cluster.
