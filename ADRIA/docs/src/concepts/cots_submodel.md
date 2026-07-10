# The ADRIA Crown-of-Thorns Starfish (COTS) Predation Submodel
## Technical Reference and Calibration Document

This document provides a comprehensive overview of the stage-structured Crown-of-Thorns Starfish (COTS, *Acanthaster planci*) predation submodel implemented in **ADRIA.jl**. It covers the model's ecological design, mathematical equations, data dependencies, and code integration details.

---

## 1. Executive Summary & Objective

The primary objective of the COTS submodel in ADRIA is to simulate realistic predator-prey dynamics between COTS and corals on a reef network. In particular, it is calibrated to reproduce the historical **~15-year boom-bust cycles** observed empirically on the Great Barrier Reef (e.g., at Lizard Island between 1985 and 2024).

The submodel tracks the local population density of COTS across three age classes at each site, resolves their feeding preference for fast-growing tabular/corymbose *Acropora* corals over slow-growing massive corals, and captures starvation feedbacks, larval connectivity, and external outbreak pulses (the "Cairns Initiation Box").

---

## 2. Model Structure and Mechanics

The COTS population is modeled using a discrete-time, stage-structured model with a yearly timestep.

### A. Stage Structure (Age Classes)
At each location, the local COTS density vector $N_t = [N_1, N_2, N_3]^T$ tracks three age classes:
1. **Age 0 (Recruits):** Larvae that have successfully settled on the reef.
2. **Age 1 (Juveniles):** Non-predatory or minimally predatory transitioning stage.
3. **Age 2+ (Adults):** The primary predatory stage that feeds on corals and reproduces.

### B. Maternal Nutrition & Body Condition
Adult fecundity and larval quality depend heavily on the availability of coral prey. The submodel tracks a maternal nutrition index called **body condition** ($BC \in [0, 1]$), modeled as an exponential moving average of food availability:

$$\text{food signal}_t = \min\left(1.0, \frac{F_t + S_t}{C_{max} \times 0.5}\right)$$

$$BC_t = (1 - \alpha) \cdot BC_{t-1} + \alpha \cdot \text{food signal}_t$$

*   $F_t$ and $S_t$ are the relative covers of fast-growing and slow-growing corals, respectively.
*   $C_{max}$ is the carrying capacity/reference coral cover level at which starvation mortality is zero.
*   $\alpha = 1/\tau_{\text{condition}}$ is the update rate (where $\tau_{\text{condition}}$ is the moving average timescale, typically set to 3.0 years).

### C. Starvation Mortality
When coral cover falls below a critical starvation threshold ($15\%$ of $C_{max}$), adult COTS experience a rapid, non-linear increase in mortality:

$$\text{If } (F_t + S_t) > 0.15 \cdot C_{max}: \quad f = 1.0$$

$$\text{If } (F_t + S_t) \le 0.15 \cdot C_{max}: \quad f = (1.0 - \tilde{p}) + \tilde{p} \left(\frac{F_t + S_t}{0.15 \cdot C_{max}}\right)^3$$

*   $f \in [1-\tilde{p}, 1]$ is the starvation survival modifier applied to juveniles and adults.
*   $\tilde{p}$ (`p_tilde`) is the maximum fraction of mortality attributable to starvation (range: $[0, 1]$).

### D. Reproduction (Ricker Dynamics & Allee Effect)
COTS recruitment combines Ricker density-dependent overcompensation with an Allee effect (representing fertilization failure when adult densities are too low to spawn successfully):

$$\text{fecundity} = a_{\text{ricker}} \cdot BC_t^2$$

$$\text{allee effect} = \frac{N_3^2}{A^2 + N_3^2}$$

$$R = \text{fecundity} \cdot N_3 \cdot \exp(-b_{\text{ricker}} \cdot N_3) \cdot \text{allee effect}$$

*   $a_{\text{ricker}}$ is the Ricker reproductive rate parameter.
*   $b_{\text{ricker}}$ is the Ricker density-dependent suppression parameter.
*   $A$ (`allee_threshold`) is the adult density at which fertilization success drops to $50\%$.

### E. Transition Equations
The COTS population vector is updated at each timestep according to:

$$\begin{aligned}
N_{t+1}[3] &= N_t[2] \cdot (s_2 \cdot f) + N_t[3] \cdot (s_3 \cdot f) \\
N_{t+1}[2] &= N_t[1] \cdot s_1 \\
N_{t+1}[1] &= R + \text{imm}_t
\end{aligned}$$

*   $s_i = 1 - m_i$ are the baseline survival rates for age class $i$.
*   $\text{imm}_t$ is the background larval immigration from external sources, gated by coral cover:
    $$\text{imm}_t = \text{IMM} \cdot \min\left(1.0, \left(\frac{F_t + S_t}{\text{imm threshold} \times C_{max}}\right)^{\eta_{\text{imm}}}\right)$$

---

## 3. Coral Predation & Prey Preference

COTS feeding is modeled via a **generalized Holling Type II/III Functional Response** that targets two broad coral functional groups:
1.  **Fast-Growing Corals ($F$):** Tabular *Acropora*, Corymbose *Acropora*, and Corymbose non-*Acropora* (Groups 1, 2, 3 in ADRIA).
2.  **Slow-Growing Corals ($S$):** Small and large massive corals (Groups 4, 5 in ADRIA).

### A. Consumption Rates
The relative coral cover consumed per adult COTS is:

$$\text{Cons}_F = N_3 \cdot \frac{a_F \cdot F^{\eta_F}}{1.0 + h \cdot (a_F \cdot F^{\eta_F} + a_S \cdot S^{\eta_S})}$$

$$\text{Cons}_S = N_3 \cdot \frac{a_S \cdot S^{\eta_S}}{1.0 + h \cdot (a_F \cdot F^{\eta_F} + a_S \cdot S^{\eta_S})}$$

*   $a_F$ and $a_S$ are the consumption rate coefficients for fast and slow corals, respectively.
*   $\eta_F$ and $\eta_S$ are exponents allowing for sigmoidal Type III responses (default is 1.0, representing Type II).
*   $h$ is the handling time coefficient (typically set to 0.0 in the current calibration).

### B. Proportional Cover Mortality
Predation reduces the relative coral cover of each size class $s$ and group $g$ proportionally:

$$\text{surv}_F = \max\left(0.0, 1.0 - \frac{\min(\text{Cons}_F, F)}{F}\right)$$

$$\text{surv}_S = \max\left(0.0, 1.0 - \frac{\min(\text{Cons}_S, S)}{S}\right)$$

$$\begin{aligned}
C_{\text{cover}, t}[g, s, \text{loc}] &\leftarrow C_{\text{cover}, t}[g, s, \text{loc}] \times \text{surv}_F \quad (\text{for } g \in \text{fast groups}) \\
C_{\text{cover}, t}[g, s, \text{loc}] &\leftarrow C_{\text{cover}, t}[g, s, \text{loc}] \times \text{surv}_S \quad (\text{for } g \in \text{slow groups})
\end{aligned}$$

---

## 4. Input Data Expectations

To run the COTS submodel, the following inputs are required or configured:

| Data Type | Description | File Path / Source |
| :--- | :--- | :--- |
| **Initial Coral Cover** | Starting relative cover for the 5 functional groups across all size classes. | Domain bundle (`init_coral_cover` matrix) |
| **Site Connectivity** | Sparse transition probability matrix representing hydrodynamic larval transport. | Domain bundle (`conn` data cube) |
| **Initial COTS Density** | Spatial suitability vector used to seed initial COTS densities. | Domain bundle (`cots_init_density` or random initialization fraction) |
| **Environmental Drivers** | Degree Heating Weeks (DHW) and Cyclone mortality scenarios. | Domain bundle (`dhw_scens` and `cyclone_mortality_scens`) |
| **COTS Parameters** | Named tuple containing rates ($m_1, m_2, m_3, a_F, a_S, \text{IMM}$, etc.). | Scenario parameter DataFrame |

---

## 5. Implementation Details in ADRIA.jl

The COTS submodel code is modular and resides in two main files:
*   **Model Definitions & Timestep:** [src/ecosystem/cots.jl](../../../src/ecosystem/cots.jl)
*   **Simulation Loop Integration:** [src/scenario.jl](../../../src/scenario.jl)

### A. Integration Hooks into ADRIA and CoralBlox

The COTS submodel is currently implemented inside ADRIA rather than as an independent package. It is still structured around a small set of explicit hooks into the main ecosystem loop:

| Hook | ADRIA location | Purpose |
| :--- | :--- | :--- |
| Model construction | `src/scenario.jl` before the timestep loop | Build one `CotsHuman` state object per site using scenario parameters and optional spatial initial densities. |
| Coral state read | `C_cover_t` inside `run_model` | Read the current relative coral cover by functional group, size class, and location. |
| Coral mortality write | `cots_mortality!(C_cover_t, cots_models, cots_prey_map)` | Apply COTS predation directly to the mutable coral cover tensor. |
| Larval dispersal | `disperse_cots_larvae!(cots_models, cots_conn)` | Move age-0 COTS recruits between sites using the domain connectivity matrix. |
| External supply | pulse block in `src/scenario.jl` | Optionally inject external larvae into configured pulse locations. |
| Logging | `Ycots`, `Ycots_bc` | Store COTS age-class densities and body condition alongside ADRIA outputs. |

The integration point is intentionally late in each ADRIA timestep. Coral growth, recruitment, interventions, bleaching, and cyclones are applied first. COTS then consume the post-disturbance coral state, so starvation and predation respond to the reef condition that remains after acute environmental stress.

The core shared state is:

```julia
C_cover_t[g, s, loc]
```

where `g` is coral functional group, `s` is size class, and `loc` is site. COTS treats this as relative coral cover. During earlier CoralBlox growth operations, ADRIA temporarily converts cover to absolute area by multiplying by habitable site area. Before COTS predation runs, ADRIA converts `C_cover_t` back to relative cover, which is the scale expected by `cots_mortality!`.

### B. How COTS Interacts with CoralBlox State

CoralBlox owns the coral demographic mechanics used by ADRIA: colony growth, size-class transitions, coral cover calculations, fecundity scope, recruitment, and helper functions such as `coral_cover`, `loc_coral_cover`, `linear_extension_scale_factors`, and `relative_leftover_space`. COTS does not call most CoralBlox APIs directly. Instead, it interacts with the CoralBlox-managed state through the already materialized ADRIA cover tensor.

The interaction is therefore a predator-prey feedback rather than a deeply coupled implementation:

1. CoralBlox/ADRIA update coral cover for growth, recruitment, interventions, bleaching, and cyclones.
2. ADRIA passes the updated relative cover tensor to COTS.
3. COTS aggregates cover into fast and slow prey pools using `CotsPreyMap`.
4. COTS updates its own age classes and body condition from those prey pools.
5. COTS returns consumption of fast and slow coral cover.
6. ADRIA applies that consumption by multiplying the relevant CoralBlox cover groups and size classes by proportional survival.
7. The modified cover tensor continues through the rest of the ADRIA timestep and is written back to `C_cover`.

This means COTS currently depends on CoralBlox concepts, but mostly through a narrow data contract:

```julia
relative_cover::Array{Float64,3}  # groups x sizes x locations
habitable_area::Vector{Float64}
connectivity::SparseMatrixCSC{Float64,Int64}
functional_group_mapping::CotsPreyMap
```

The most important convention is group ordering. The current prey map assumes ADRIA's five coral functional groups:

- fast prey: groups `1, 2, 3`
- slow prey: groups `4, 5`

If CoralBlox group definitions change, or if another domain uses different taxa ordering, the `CotsPreyMap` must be supplied explicitly rather than assumed.

### C. State Ownership and Side Effects

COTS owns these state variables:

- `N[1]`: age-0 recruits
- `N[2]`: age-1 juveniles
- `N[3]`: age-2+ adults
- `body_condition`: maternal nutrition state
- COTS-specific parameter tuple (`m1`, `m2`, `m3`, `a_F`, `a_S`, `IMM`, `a_ricker`, etc.)

ADRIA/CoralBlox owns the coral state. COTS mutates coral cover only through `cots_mortality!`, which is deliberately analogous to other ADRIA disturbance functions such as cyclone mortality: it takes the mutable cover tensor, calculates proportional survival, applies mortality, clamps to valid bounds, and returns `nothing`.

The current side-effect boundaries are:

- `cots_timestep!` mutates only a single `CotsHuman` object and returns consumption values.
- `cots_mortality!` mutates `C_cover_t` and all `cots_models` because predation and COTS population dynamics are coupled in the same pass.
- `disperse_cots_larvae!` mutates age-0 recruits after local recruitment.
- `inject_upstream_pulse!` mutates age-0 recruits at selected pulse locations.

Because these are mutable operations, any future parallelization over locations should preserve the current ordering: local COTS timestep first, larval dispersal from a snapshot of age-0 recruits second, external pulse third.

### D. The Simulation Timestep Order of Events
Within the main time loop in `scenario.jl` (lines 1192–1959), the nominal sequence of events for a single timestep is as follows:
1.  **Coral Spawn and Recruit:** Coral recruits settle based on October/November spawning.
2.  **Solar Radiation Management (SRM) & Shading:** Temperature stresses are adjusted.
3.  **Coral Seeding and Relocation (Moving Corals):** Interventions are deployed.
4.  **Bleaching Stress:** Bleaching mortality is applied (typically November–February).
5.  **Cyclones:** Cyclone mortality is applied (typically January–March).
6.  **COTS Predation Mortality:** COTS feed on remaining coral cover via `cots_mortality!`.
7.  **COTS Larval Dispersal:** Recruits disperse between reefs via `disperse_cots_larvae!`.
8.  **External Larval Pulse:** If enabled, external Cairns initiation pulses are injected.
9.  **Logging:** COTS populations and body conditions are logged for analysis.

### E. COTS Larval Dispersal Method
Larval transport between reefs is modeled using the sparse hydrodynamic connectivity matrix. To prevent double-counting of local settlement, self-connectivity is excluded, and cross-reef dispersal is scaled by an immigration multiplier to account for COTS' massive fecundity:

```julia
function disperse_cots_larvae!(
    cots_models::Vector{CotsHuman},
    conn::SparseMatrixCSC{Float64, Int64};
    immigration_scalar::Float64=1.0
)
    n_locs = length(cots_models)
    local_recruits = [cots_models[loc].N[1] for loc in 1:n_locs]

    for sink in 1:n_locs
        immigration = 0.0
        for k in SparseArrays.nzrange(conn, sink)
            src = SparseArrays.rowvals(conn)[k]
            if src == sink
                continue  # Skip self-connectivity (handled by local Ricker model)
            end
            p = SparseArrays.nonzeros(conn)[k]
            immigration += local_recruits[src] * p * immigration_scalar
        end
        cots_models[sink].N[1] = max(0.0, local_recruits[sink] + immigration)
    end
end
```

### F. External Larval Pulse Boundary Condition
External larval pulses are **off by default**. They are a diagnostic boundary
condition for mimicking upstream COTS larval supply before an explicit external
source model is available.

When enabled, pulse timing is configured with a start timestep, duration, and
optional repeat interval:

| ENV variable | Default | Meaning |
| :--- | :---: | :--- |
| `COTS_EXTERNAL_PULSE` | `false` | Enables/disables the pulse mechanism. |
| `COTS_PULSE_START` | `1` | First model timestep where the pulse window can be active. |
| `COTS_PULSE_DURATION` | `1` | Number of consecutive timesteps in each active pulse window. |
| `COTS_PULSE_REPEAT_INTERVAL` | `0` | Repeat interval in timesteps. `0` means no repeat. |
| `COTS_PULSE_RELATIVE_MAGNITUDE` | `0.0` | Pulse size as a proportion of each location's running maximum internal age-0 larval supply. |

The relative scaling is site-specific. After local recruitment and larval
dispersal, ADRIA tracks the running maximum age-0 larval supply at each site.
During an active pulse window, the external pulse adds:

```julia
pulse_vals = COTS_PULSE_RELATIVE_MAGNITUDE .* cots_max_larval_supply
inject_upstream_pulse!(cots_models, pulse_locs, pulse_vals)
```

This lets calibration experiments vary external supply timing and magnitude
without hard-coding absolute larval densities. It is intended as a temporary
proxy for upstream supply until a real external COTS source/connectivity layer is
implemented.

---

## 6. Current Best-Fit Calibration Parameters

Through a 250-scenario Latin Hypercube Sweep (LHS) calibrated against empirical tow data at Lizard Island, the following parameters have been identified as the optimal combination (`run_id = 2`) to capture both the peak and the subsequent population bust:

*   **Fast-Coral Starvation Threshold ($a_F$):** `1.34378`
*   **Slow-Coral Starvation Threshold ($a_S$):** `0.0957831`
*   **Background Immigration Rate ($IMM$):** `0.0947791`
*   **Initial Seed Multiplier:** `2.23695`
*   **External Pulse:** off by default; configure with `COTS_EXTERNAL_PULSE=true` and `COTS_PULSE_RELATIVE_MAGNITUDE` for diagnostics.

> [!NOTE]
> The higher consumption rate of fast corals ($a_F \approx 1.34$) combined with a low slow coral consumption rate ($a_S \approx 0.096$) ensures COTS rapidly deplete Acropora during outbreaks, triggering a sharp and realistic starvation crash.

---

## 7. If COTS Became a Separate Package

A future `COTSMod.jl` package, analogous to `CoralBlox`, should separate COTS ecology from ADRIA scenario orchestration. ADRIA would remain responsible for domain loading, scenario sampling, intervention policy, environmental drivers, and output collation. The COTS package would own COTS state transitions, predation calculations, larval dispersal, and external supply boundary conditions.

### A. Proposed Package Boundary

A clean package boundary would expose a small API:

```julia
cots = COTSMod.initialize(n_locs, params; init_density, seed_locs, rng)
consumption = COTSMod.step!(cots, prey_state, env_state)
COTSMod.apply_predation!(coral_state, consumption, prey_map)
COTSMod.disperse_larvae!(cots, cots_connectivity; scalar=1.0)
COTSMod.apply_external_supply!(cots, supply_config, supply_state)
```

ADRIA would call those functions from `scenario.jl`, but would not need to know the details of Ricker recruitment, starvation functions, Allee effects, or external pulse scaling.

The package should define explicit data structures for:

- `CotsParams`: biological parameters and bounds.
- `CotsState`: age-class densities and body condition by location.
- `CotsPreyMap`: mapping from coral groups to fast/slow prey categories.
- `CotsPreyState`: fast and slow coral cover by location.
- `CotsConnectivity`: larval dispersal matrix and any COTS-specific scaling.
- `ExternalSupplyConfig`: start, duration, repeat interval, magnitude, and target locations.
- `CotsOutput`: age-class logs, condition logs, and optional diagnostic supply logs.

### B. Data Contract with CoralBlox

A separate package should not depend on ADRIA internals such as `Domain`, `DataFrameRow`, or `YAXArray`. It should also avoid assuming a fixed five-group coral layout. Instead, ADRIA should translate CoralBlox/Domain state into a minimal prey-state interface:

```julia
struct CotsPreyState
    fast_cover::Vector{Float64}
    slow_cover::Vector{Float64}
    total_cover::Vector{Float64}
end
```

or, when size-class detail is needed:

```julia
relative_cover::AbstractArray{Float64,3}  # groups x sizes x locations
prey_map::CotsPreyMap
```

This keeps CoralBlox as the owner of coral demography while allowing COTS to operate on generic coral-cover data. ADRIA would be the adapter layer that knows how to convert CoralBlox functional groups into COTS prey categories.

### C. What ADRIA Would Still Own

Even if COTS became `COTSMod.jl`, ADRIA should still own:

- domain loading and validation;
- scenario parameter tables and sampling;
- time-loop orchestration and event ordering;
- coupling to bleaching, cyclones, SRM, seeding, and moving corals;
- conversion between absolute cover and relative cover;
- mapping of reef/site IDs to external pulse target locations;
- result storage, Zarr batching, and post-processing metrics.

In that design, ADRIA's timestep would read more like an orchestration layer:

```julia
CoralBlox.step!(coral_state, coral_env, interventions)
Disturbances.apply!(coral_state, disturbance_state)
prey_state = ADRIA.to_cots_prey_state(coral_state, prey_map)
consumption = COTSMod.step!(cots_state, prey_state)
ADRIA.apply_cots_consumption!(coral_state, consumption, prey_map)
COTSMod.disperse_larvae!(cots_state, cots_connectivity)
COTSMod.apply_external_supply!(cots_state, external_supply)
```

### D. Why This Separation Would Help

Separating COTS would make it easier to:

- test COTS population dynamics independently of a full ADRIA run;
- run toy predator-prey systems for calibration and sensitivity analysis;
- swap coral connectivity for COTS-specific hydrodynamic connectivity;
- support multiple external source models instead of ENV-driven diagnostics;
- document biological assumptions separately from ADRIA orchestration;
- reuse the COTS model with other coral models that can provide a compatible prey-state interface.

### E. Risks and Design Constraints

The main risk is hiding important coupling details. COTS predation timing matters because it occurs after bleaching and cyclone mortality in the current ADRIA loop. A package boundary should therefore avoid owning the global timestep order. It should provide deterministic, side-effect-clear operations that ADRIA can call at the correct point.

The package should also avoid hard-coded assumptions about:

- five coral functional groups;
- Lizard Island-specific seed locations;
- coral connectivity being a valid COTS larval dispersal proxy;
- annual timesteps only;
- ENV variables as the primary configuration mechanism.

A practical migration path would be:

1. Move `CotsHuman`, `CotsParams`, `CotsPreyMap`, and pure COTS functions into a package module.
2. Keep ADRIA's `scenario.jl` integration unchanged except for qualified calls such as `COTSMod.cots_mortality!`.
3. Replace ENV pulse controls with a typed `ExternalSupplyConfig` passed from ADRIA.
4. Add package-level unit tests for population transitions, starvation, predation, larval dispersal, and external supply.
5. Add ADRIA integration tests that verify event ordering and coral-cover feedbacks are unchanged.

---

## 8. Multi-Scale Application: Local vs. GBR-Wide

The COTS submodel in ADRIA.jl is designed to be **domain-agnostic**. The core ecosystem model loop reads parameters and dimensions directly from whichever Domain object is loaded. This allows the submodel to run on both localized clusters (like Lizard Island) and GBR-wide domains.

### A. How the GBR-Wide Model Works

When loading the GBR-wide model via `RMEDomain` (ReefMod Engine Domain), the submodel scales as follows:

1.  **Scale-Free Loop Execution:** The simulation loop in `scenario.jl` determines the number of reefs dynamically using `n_locs = domain.coral_growth.n_locs`. For `RMEDomain`, this covers all **3,800+ reefs** of the Great Barrier Reef.
2.  **Sparse Dispersal Scaling:** Larval transport is driven by a GBR-scale connectivity matrix. Since ADRIA uses sparse matrix multiplication (`SparseMatrixCSC` from Julia's `SparseArrays`), COTS larval dispersal (`disperse_cots_larvae!`) only calculates active connections, avoiding the $O(N^2)$ penalty and maintaining high performance at GBR scale.
3.  **Spatial Initializations:** Initial COTS densities can be seeded using the GBR-wide observed spatial suitability vector `domain.cots_init_density` via `init_cots_from_spatial()`, scaling initial cohort densities based on reef-specific habitat quality.

### B. GBR-Wide Calibration Challenges (Not the Current Focus)

While the code fully supports GBR-wide simulation, it is not the primary focus of the current phase and requires extensive calibration due to several key factors:

*   **Outbreak Waves and EAC Propagation:** Outbreaks on the GBR do not happen simultaneously. Historically, primary outbreaks initiate in the northern GBR (Cairns/Moore reefs) and propagate southward as larval waves carried by the East Australian Current (EAC) over 10-15 years. Modeling this requires calibrating the spatial path of these waves rather than assuming synchronous local dynamics.
*   **Connectivity and Dispersal Scaling:** The `immigration_scalar` multiplier (representing COTS' high larval production relative to corals) is highly sensitive to domain resolution. A scalar calibrated for a local 10-reef cluster will over-saturate or under-represent larval supply when scaled up to 3,800+ reefs with larger physical distances.
*   **External Larval Pulse Targeting:** The Cairns Initiation Box external pulse must be configured to inject larvae only at the specific northern initiator reefs, requiring GBR-specific mapping of `_cots_seed_locs` instead of local domain index sets.
*   **Regional Parameter Heterogeneity:** Starvation thresholds ($a_F, a_S$) and baseline survival rates may vary along the latitudinal gradient of the GBR due to water temperature, background predator density (e.g., giant triton snails), and coral community composition, requiring a spatially varying parameterization.

