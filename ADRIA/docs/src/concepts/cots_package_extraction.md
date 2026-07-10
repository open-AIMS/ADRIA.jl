# COTS Package Extraction Plan

This document defines the migration path from the current ADRIA-owned COTS submodel to a separate `COTSMod.jl` package, plus a separate reproducible calibration-study repository.

## Goals

The extraction has two goals:

1. Move the COTS ecological model into a small Julia package that ADRIA can call in the same style as CoralBlox.
2. Move Lizard Island calibration scripts, plotting, and validation workflows into a standalone reproducible study repository.

The migration should preserve current ADRIA behavior at each step. Do not move code before the behavior contract is covered by tests.

## Current Ownership

COTS is currently split across these ADRIA files:

| Area | Current file | Should move to COTSMod? |
| :--- | :--- | :---: |
| COTS age-class state | `ADRIA/src/ecosystem/cots.jl` | Yes |
| COTS local timestep | `ADRIA/src/ecosystem/cots.jl` | Yes |
| Starvation, body condition, recruitment | `ADRIA/src/ecosystem/cots.jl` | Yes |
| Prey mapping type | `ADRIA/src/ecosystem/cots.jl` | Yes |
| Predation calculation | `ADRIA/src/ecosystem/cots.jl` | Yes |
| Direct mutation of ADRIA coral cover tensor | `ADRIA/src/ecosystem/cots.jl` | Maybe, but likely keep adapter in ADRIA |
| Larval dispersal | `ADRIA/src/ecosystem/cots.jl` | Yes |
| External pulse injection | `ADRIA/src/ecosystem/cots.jl` and `scenario.jl` | Yes, but with typed config |
| ADRIA parameter factors | `ADRIA/src/ecosystem/cots_factors.jl` | Split: package defaults plus ADRIA factor adapter |
| Scenario parameter extraction | `ADRIA/src/scenario.jl` | No |
| Domain-specific initialization | `ADRIA/src/scenario.jl`, `LizardDomain.jl` | No |
| Result logging and Zarr output | `ADRIA/src/scenario.jl`, `io/*` | No |
| Calibration scripts | `sandbox/calibration` | Study repo |
| Plots and validation metrics | `sandbox/plotting` | Study repo |

## Proposed COTSMod API

A first package API should be intentionally small:

```julia
module COTSMod

export CotsParams, CotsState, CotsPreyMap, CotsPreyState
export initialize_cots, step_cots!, apply_predation!, disperse_larvae!, apply_external_supply!

end
```

Recommended core types:

```julia
struct CotsParams
    a::Float64
    b::Float64
    IMM::Float64
    p_tilde::Float64
    C_max::Float64
    m1::Float64
    m2::Float64
    m3::Float64
    a_F::Float64
    a_S::Float64
    h::Float64
    eta_F::Float64
    eta_S::Float64
    eta_starve::Float64
    eta_imm::Float64
    imm_threshold::Float64
    fecundity_gate::Bool
    a_ricker::Float64
    b_ricker::Float64
    tau_condition::Float64
    allee_threshold::Float64
end

struct CotsPreyMap
    fast_group_indices::Vector{Int}
    slow_group_indices::Vector{Int}
end

struct CotsPreyState
    fast_cover::Vector{Float64}
    slow_cover::Vector{Float64}
end
```

The package should not depend on ADRIA, YAXArrays, DataFrames, Documenter, Zarr, or domain objects. Initial dependencies should be limited to `StaticArrays` plus standard libraries such as `Random` and `SparseArrays`.

## ADRIA Adapter Boundary

ADRIA should keep the orchestration layer. In `scenario.jl`, ADRIA should:

1. Read sampled COTS parameters from `param_set`.
2. Convert them into `COTSMod.CotsParams`.
3. Build a `COTSMod.CotsPreyMap` from the active coral functional-group ordering.
4. Initialize COTS state using domain initial density or seed locations.
5. At the correct timestep point, pass the current relative coral cover to COTSMod.
6. Apply returned predation to the ADRIA/CoralBlox cover tensor.
7. Run larval dispersal using selected COTS connectivity.
8. Apply optional external supply.
9. Log COTS age classes and body condition into ADRIA outputs.

The important behavior to preserve is event order: COTS predation currently runs after bleaching and cyclone mortality and before COTS larval dispersal and external supply.

## Test Contract Before Moving Code

Before extracting a package, these tests should pass and be copied into the new package:

- `CotsHuman` or `CotsState` construction preserves age classes and body condition.
- One local timestep keeps age classes finite and non-negative.
- Starvation reduces adults under zero coral cover.
- Allee effect suppresses recruitment at low adult density.
- Predation reduces fast and slow coral cover without producing negative cover.
- Larval dispersal skips self-connectivity and adds cross-reef recruits.
- Explicit seed locations initialize only intended locations.
- Scalar external pulse adds recruits only to target locations.
- Vector external pulse adds location-specific recruits.

ADRIA should retain an integration smoke test that runs the Lizard historical domain for a tiny scenario count and verifies that `cots_log` and `cots_condition_log` are present.

## Migration Steps


Status update: ADRIA now includes a package-shaped internal COTS API: `initialize_cots`, `apply_predation!`, `disperse_larvae!`, and `apply_external_supply!`. `scenario.jl` calls these wrappers instead of the legacy implementation names. The legacy functions remain in place for compatibility and as implementation details until the external package exists.

Verification passed:

```powershell
julia --project=sandbox test\ecosystem\cots.jl
$env:COTS_N_STOCHASTIC_SCENS='2'; $env:COTS_OUTPUT_TAG='adapter_smoke'; $env:COTS_EXTERNAL_PULSE='false'; julia sandbox\calibration\simulate_best_stochastic.jl
```

Latest focused COTS test result: 54 passing assertions. The scenario smoke wrote `sandbox/data/best_stochastic_adapter_smoke_*` outputs.
### Step 1: Stabilize Contract Inside ADRIA

- Keep current files in ADRIA.
- Expand `test/ecosystem/cots.jl` to include vector pulse behavior and a small cover-mutation regression.
- Remove any direct `ENV` reads from core COTS functions where feasible; keep ENV parsing in `scenario.jl`.
- Introduce a package-like namespace or adapter boundary if useful, but avoid large code movement.

### Step 2: Create COTSMod.jl Skeleton

- Create a new Julia package with `Project.toml`, `src/COTSMod.jl`, and `test/runtests.jl`.
- Copy the tested biological code from ADRIA.
- Convert ADRIA factor definitions into package defaults plus an ADRIA-specific factor adapter.
- Run package tests without ADRIA loaded.

### Step 3: Make ADRIA Depend on COTSMod

- Add `COTSMod` as an ADRIA dependency.
- Replace unqualified calls in `scenario.jl` with `COTSMod` calls.
- Keep ADRIA result names unchanged: `cots_log`, `cots_condition_log`.
- Run the post-rebase smoke and compare trajectories against a saved baseline.

### Step 4: Remove Duplicated ADRIA COTS Internals

- Delete or deprecate `ADRIA/src/ecosystem/cots.jl` once ADRIA calls the package.
- Keep `ADRIA/src/ecosystem/cots_factors.jl` only if ADRIA still needs factors for sampling.
- Update docs to link to COTSMod for biological model details.

### Step 5: Extract Calibration Study Repo

The calibration repo should contain:

- `Project.toml` and manifest or lockfile.
- `scripts/calibration/*.jl` from current `sandbox/calibration`.
- `scripts/plotting/*.py` from current `sandbox/plotting`.
- `data/validation/reef_cots.csv`, `reef_manta.csv`, and reef/site mapping tables.
- Lightweight domain-building instructions or data-release download instructions.
- `README.md` with commands to regenerate calibration sweeps, stochastic plots, pulse scenarios, and validation metrics.
- `outputs/` ignored by default, with selected figures committed only when intended.

Do not move large raw RME/vendor data into the study repo until licensing and storage expectations are explicit.

## Decision Points

Before creating the external repos, decide:

- Package name is currently `COTSMod.jl`.
- Whether `CotsParams` in the package is a plain struct or remains compatible with ADRIA `Factor` definitions.
- Whether predation should mutate coral cover inside COTSMod or return a survival/consumption object for ADRIA to apply.
- Whether external pulses are a generic external supply model or a Lizard-specific calibration diagnostic.
- How large input data are versioned for the calibration study.