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

### A. The Simulation Timestep Order of Events
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

### B. COTS Larval Dispersal Method
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

### C. Cairns Initiation Box External Larval Pulse
To simulate large-scale outbreak waves that originate in the northern GBR (Cairns region) and propagate southward along the East Australian Current, an external pulse boundary condition is checked at each timestep:

```julia
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

---

## 6. Current Best-Fit Calibration Parameters

Through a 250-scenario Latin Hypercube Sweep (LHS) calibrated against empirical tow data at Lizard Island, the following parameters have been identified as the optimal combination (`run_id = 2`) to capture both the peak and the subsequent population bust:

*   **Fast-Coral Starvation Threshold ($a_F$):** `1.34378`
*   **Slow-Coral Starvation Threshold ($a_S$):** `0.0957831`
*   **Background Immigration Rate ($IMM$):** `0.0947791`
*   **Initial Seed Multiplier:** `2.23695`
*   **Cairns Pulse Amplitude (`pulse_val`):** `1.5`

> [!NOTE]
> The higher consumption rate of fast corals ($a_F \approx 1.34$) combined with a low slow coral consumption rate ($a_S \approx 0.096$) ensures COTS rapidly deplete Acropora during outbreaks, triggering a sharp and realistic starvation crash.
