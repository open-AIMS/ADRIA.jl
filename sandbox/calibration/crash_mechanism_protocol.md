# Preregistered 1991 peak–crash mechanism screen (2026-10-08)

Status: **both bounded factorials completed and failed; no default or calibration-objective change**. Use the archived [1991 full-cohort diagnostic](runs/20261007T199111_owen_1991_cohort_history_final) as the exact control, with 1991 IDW adults at the provisional `q=0.015 COTS/tow per adult COTS/ha`, assumed juvenile:subadult:adult seed `6:3:1`, inherited coral, Owen annual sampled connectivity and boundary, seed `20260930`, `3 adult COTS/ha` Allee half-saturation, and the frozen 1992–2024 score/reef aggregation. `q` and stage preload are **not** biological estimates. The first year of `cots_log` is an initialization placeholder and is excluded from scoring.

## Hypotheses and mathematical changes

The existing adult state survives and grazes during the inter-wave gap. Existing food dependence acts in the same year as coral cover, and its mortality can spare enough coral to extinguish its own trigger. We need to distinguish delayed juvenile release from a loss of established adults, and only then test a lagged outbreak burden.

Let `J_t` be the model's immature stage `N[2]` (COTS/ha), `A_t=N[3]` (adult COTS/ha), `C_t` the model's food-coral fraction, `f(C_t)` existing food survival, and `m2,m3` annual background mortality fractions. The existing opt-in storage arm uses

`g(C_t)=clamp(g_min+(g_max-g_min)*sqrt(clamp(C_t/C_g,0,1)),0,1)`,

`J^surv_t=J_t(1-m2)f(C_t)`, `M_t=g(C_t)J^surv_t`,

`J_{t+1}=N[1]_t(1-m1)+J^surv_t-M_t`, `A_{t+1}=M_t+A_t(1-m3)f(C_t)`.

The switch is already implemented in COTSMod and defaults **off**. It delays stage transition but does not remove adults; stored juveniles can generate a later rebound. For this diagnostic use existing `g_min=0.05`, `g_max=1`, `C_g=0.40` coral fraction. These are screen settings, not evidence-based bounds. The comparator sets `m3=0.30/year` versus the archived candidate's `m3`, without altering production or consumption. A high constant mortality may attenuate the first peak, which is precisely why it is crossed with stage delay rather than treated as an answer.

If the 2×2 fails, a *separate, opt-in* candidate may add a lagged outbreak burden `B` and an adult hazard, with independently named parameters:

`B_{t+1}=rho*B_t+(1-rho)*A_t` (COTS/ha),

`d_t=max(B_t-B_crit,0)/(w_B+max(B_t-B_crit,0))`,

`h_t=h_max*d_t*[1 or sigma((C_crit-C_t)/w_C)]` (dimensionless one-year cumulative hazard),

`A_{t+1}=M_t+A_t(1-m3)f(C_t)*exp(-h_t)`.

`sigma` is logistic. The density factor is zero at and below its independent threshold, so the preload cannot incur an artificial background hazard before burden accumulates. Cross `h_t` with a density-only arm that omits the coral multiplier, so food interaction is identifiable. `B_crit` is **not** the Allee threshold. Record burden, hazard, added adult deaths and source/settlement fluxes by reef and site. No calendar-year or time-since-1991 trigger is eligible for promotion: it could reproduce observed peaks by construction. A *state-triggered outbreak-duration clock* is eligible as an explicitly phenomenological diagnostic comparator, subject to the same promotion gates as other mechanisms. Experimental disease induction in COTS supports biological plausibility of disease-related mortality, but does not identify a natural density threshold or prove that pathogens caused the observed Lizard trough ([Rivera-Posada et al., 2012](https://pubmed.ncbi.nlm.nih.gov/22303625/)). This proposed hazard is therefore **pathogen-compatible phenomenology**, not a claim of observed epizootics.

## Smallest factorial and frozen readout

First screen (one seed, one historical anchor; no optimizer):

| Arm | Juvenile storage | `m3` annual fraction | Distinguishes |
| --- | --- | --- | --- |
| Archived-equivalent control | off | archived value | Exact replay / provenance |
| Stage gate | on | archived value | Delayed maturation alone |
| Constant mortality | off | 0.30 | Established-adult attrition alone |
| Both | on | 0.30 | Whether immature storage rescues a mortality-damped second peak |

All arms retain the same `COTS/tow per adult/ha` observation scale, raw and three-year-smoothed observations, area-weighted reef mapping, weather and sampled connectivity years, local/external production, fertilization, food and coral transitions. Record annual adults (ha^-1 and CPUE), coral fraction, maturation, retained juveniles, local production, internal immigration, external immigration, and observed/treatment peak diagnostics. Do not refit `q`, scale, Allee, initial cohorts or boundary per arm.

The qualitative screen requires a first and second modeled peak within the existing detection/timing rules on Lizard, MacGillivray and North Direction; all seven observed peaks matched overall; worst inter-peak smoothed minimum divided by the smaller adjacent modeled peak `<=0.10`; and no material reduction in peak height or worsening of coral against the paired control. Report each component rather than hiding a tradeoff in a scalar. Coral comparisons remain a guardrail because the 9 m manta observations and whole-reef model coral are not identical supports. Any invalid control replay, failed run, non-finite state, mass-balance error, lost second peak, or improvement caused only by extinguishing external supply is a failure. An improved one-seed arm is only a *candidate*: repeat at multiple recorded hydro/stochastic seeds, provisional `q` and cohort sensitivities, and held-out reefs/years before promotion.

## First-factorial decision and next gate

The [completed 1991 cohort factorial](runs/20261008T_07_crash_factorial_1991_v1/metadata.toml) passed an exact archived state/flux control replay and invariant external-supply audit. None of the three treatments met the qualitative gate. Stage storage delayed Lizard's first peak ten years and worsened the worst trough ratio (`0.594→0.691`); `m3=0.30` delayed it eight years and cut its height by 56% while the worst trough remained `0.515`. The combination lost Lizard's first detected wave altogether. See [per-reef diagnostics](runs/20261008T_07_crash_factorial_1991_v1/crash_factorial_audit.csv) and [paired plot](runs/20261008T_07_crash_factorial_1991_v1/cots_coral_crash_factorial.png). **No arm is promotable.** The observed mechanism is compensation: storage retains immature COTS and raises later local/internal recruitment; constant mortality acts before as well as after the desired first peak. The model coral discrepancy persists, with the stated observation-support caveat.

For this **diagnostic capacity test**, hold `rho=0.5` (one-year memory weight), `B_crit=1.5 adult COTS/ha`, `w_B=0.5 adult COTS/ha`, `C_crit=0.15` coral fraction and `w_C=0.03` coral fraction. Cross density-only and density×food with `h_max` of `0.75` and `1.5` per annual step, plus a zero-hazard control: five arms, one recorded hydrodynamic seed, no optimizer. The higher hazard corresponds to at most `1-exp(-1.5)=0.777` extra loss of *background-and-food survivors* in one year, not 77.7% of all adults. These are **stress-test bounds, not field-estimated pathogen parameters**. They were selected before implementation: the archived 1991 site-level adult seed has median `0.363`, 95th percentile `0.993`, and 99th percentile `3.039 COTS/ha`; modeled control reef adult 90th percentiles over 1992–2024 are `2.410` at Lizard and `2.891` at MacGillivray. Thus `B_crit` is above most seed sites but within simulated wave densities, and is deliberately unrelated to the `3 COTS/ha` fertilization half-saturation. The new run must retain site-level burden/hazard/death summaries because archived site-level *dynamic* abundance was not saved; these settings cannot be called biological bounds.

Before implementing the lagged hazard, write pure core tests for (a) zero hazard reproducing legacy transitions bit-for-bit; (b) burden recurrence and units, including `B_0=0`; (c) increasing hazard only when the configured density/food trigger increases; (d) nonnegative, mass-accounted adult survivors and added-death diagnostics; and (e) deterministic replay. Add ADRIA adapter tests for switch, validation, dimensions and per-site logging. Hold the frozen score and boundary; record first-wave peak preservation, duration of the gap, second-wave reseeding, coral and fluxes. If the hazard fails because external settlement imposes a floor, test a separate labelled boundary sensitivity rather than silently suppressing outside-domain production. No generic outbreak timer is eligible as a biological mechanism.

## Lagged-hazard decision and next structural question

The [five-arm hazard run](runs/20261008T_08_lagged_hazard_1991_v1/metadata.toml) and [independent audit](runs/20261008T_08_lagged_hazard_1991_v1/lagged_hazard_audit.json) replay the archived control COTS/coral and all 12 recruitment-flow channels exactly. All arms retain identical external immigration. The core unit tests passed `128/128`; the ADRIA focused tests passed `23/23`. The site-level audit checks `B_0=0`, the annual recurrence to `7.1e-15 COTS/ha`, nonnegative states, and that extra deaths do not exceed pre-transition adults. See [reef outcomes](runs/20261008T_08_lagged_hazard_1991_v1/lagged_hazard_audit.csv), [site-year activation](runs/20261008T_08_lagged_hazard_1991_v1/site_hazard_summary.csv), and [COTS/coral plot](runs/20261008T_08_lagged_hazard_1991_v1/cots_coral_lagged_hazard.png).

No arm passes. Control worst two-wave trough/smaller-peak ratio is `0.593935` with 6/7 observed peaks matched. Density-only `h_max=1.5` gives `0.600000` and 6/7, delaying Lizard's first detected peak from 1997 to 2001. Density×food at `1.5` gives `0.593493` and 6/7, also moving Lizard's first peak to 2001; the tiny trough change is not an improvement. At `0.75`, both mode variants lose Lizard's first detected wave and match only 5/7. Coral is spared slightly, but remains low relative to 9 m observations, with the acknowledged support difference.

This is a *phase failure*, not a silent implementation failure. In the strong density-only treatment, hazard-active sites comprise roughly `33%` in 1995, `51%` in 2005, but `0.5%` in 2008 and `0.3%` in 2012. The mechanism acts near the first wave and switches off in the gap; adult density and continuing settlement then maintain the floor. The food multiplier does not rescue timing because the density trigger has already collapsed. An offline recurrence check on the **fixed control trajectory** (not a new ecological simulation) shows that changing only `rho` from `0.5` to `0.75` or `0.9` at `B_crit=1.5` leaves just `0.7%` or `1.0%` of sites above threshold by 2012; `rho=0.9` also delays substantial activation past the first wave. Do not launch a broad optimizer on this formulation. A future bounded test needs a separately justified persistence/recovery state or boundary mechanism that operates *after* the first peak while allowing later reseeding.

## Cohort audit and ReefMod duration precedent (2026-10-09)

The [readout-only 1991 replay](runs/20261009T_01_cohort_audit_1991_v1/metadata.toml) reproduces archived adult/coral trajectories and all 12 flux channels exactly. Its [annual area-weighted cohorts](runs/20261009T_01_cohort_audit_1991_v1/cohorts.csv) and [Lizard plot](runs/20261009T_01_cohort_audit_1991_v1/cohort_floor_diagnostics.png) show that the 1999–2007 sawtooth coincides with the realized food-survival factor alternating between about `0.27` and `1.0`. Food survival already multiplies N2-to-adult maturation and adult carryover, but does not multiply N1 recruits. In 2008–12, Lizard adult food survival averages `0.986`; adult abundance averages `0.805 COTS/ha`, of which `0.667` is carryover and `0.138` is newly matured. External settlement averages `0.542 COTS/ha/year` in those years. Thus maturation contributes to the floor but does not dominate it, and ongoing settlement maintains the future cohort pipeline. The fixed Lizard plotting scale is `0.2942 COTS/tow per adult COTS/ha`; `0.1 COTS/tow` corresponds to `0.340 adult COTS/ha` under that mapping. The provisional `q=0.015` used for 1991 initialization is distinct.

The published [ReefMod 7.0 settings](https://github.com/ymbozec/REEFMOD.7.0_GBR/blob/main/settings/settings_COTS.m) specify a randomly checked `4–6` **year** outbreak duration and comment a fixed two-year Lizard calibration option. The [transition code](https://github.com/ymbozec/REEFMOD.7.0_GBR/blob/main/functions/f_runmodel.m) confirms that adult density must remain above the threshold over the checked duration; it redraws the candidate duration at a qualifying check, then resets older age classes to background while leaving younger recruits. A [regional report](https://www.barrierreef.org/uploads/CCIP-R-04-Final-Report-Regional-Modelling.pdf) describes a different-version `2–5`-year disease window, an eight-year maximum age, and a separate `5%` preferential-prey starvation trigger. These are version-specific modelling conventions, not field estimates of natural pathogen timing. ReefMod's densities and reset values are in COTS per `400 m²` grid, so numerical transfer to this model's `COTS/ha` state and provisional CPUE mapping needs explicit justification.

The next bounded comparator should be opt-in and state-triggered, with a frozen density threshold in adult `COTS/ha`: control versus adult-only reset versus adult-plus-N2 reset at two preregistered durations. Log clock age, threshold crossings, pre/post state, removed adults and N2, surviving N1, later maturation, all settlement sources, coral, and peak/trough metrics. The adult-plus-N2 arm tests whether immature refill defeats adult-only collapse; keep external N1 supply unchanged so reseeding remains a constraint. Before implementation, specify whether the trigger requires **continuous** above-threshold density as in ReefMod or simply starts a countdown on outbreak onset; those are different hypotheses. Do not fit the clock to observed peak years or call a successful timer proof of epizootics. Failure includes losing the 1997 Lizard peak, a one-year-only trough followed by immediate N2 refill, an absent second peak, or unacceptable coral/flux degradation. Any one-seed improvement remains experimental pending replicated and held-out gates.

Rollback is `COTS_LAGGED_ADULT_HAZARD=false` (the default) or `COTS_HAZARD_MAX=0`; core COTS states, observation units, the Allee threshold, and legacy defaults are unchanged. The mechanism remains experimental and **must not** be described as confirmed pathogen mortality.

## Preregistered structural comparison: preferred-prey starvation × adult senescence (2026-10-09)

Status: **preregistered before implementation; no default or objective change.** This supersedes the ReefMod timer as the *next* screen. The timer remains eligible as a later phenomenological comparator. The reason is that ReefMod adds its timer on top of an eight-year maximum age and preferred-prey starvation. All three published GBR comparators share those two foundations, and this model lacks both: ReefMod 7.0, [CoCoNet v3.4](https://research.csiro.au/coconet/wp-content/uploads/sites/486/2025/11/CoCoNet-user-guide-and-technical-summary-v3.4.pdf) Table A.1, and MICE ([Rogers & Plagányi 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9085818/)). Porting only the timer would borrow the least mechanistic part.

### Structural diagnosis motivating the test

Values below are from the archived anchor (`a_F=1.18`, `a_S=0.165`, `C_max=0.498`, `m3=0.195`) and the [cohort audit](runs/20261009T_01_cohort_audit_1991_v1/cohorts.csv).

1. **Starvation responds to total cover.** Legacy food survival uses total prey cover `F+S` with a threshold at `0.15·C_max = 0.075`. Massive (`S`) cover is grazed at `a_S≪a_F`, so `F+S` stays above 0.075 after fast coral is depleted. Lizard adult food survival is `0.986` in 2008–12. ReefMod (5% preferred cover), CoCoNet (mortality ∝ `1/C^f`) and MICE (logistic in preferred cover relative to `K`) all trigger on *preferred* coral.
2. **Adults never age out.** Adults are one class with constant `m3`. Without inflow, a decline to 10% of peak takes `ln 0.1 / ln 0.805 ≈ 10.6` years. The floor approximates `maturation/(1-s3·f)`, about `0.14/0.2 ≈ 0.7 COTS/ha`, matching the audited floor. CoCoNet applies `M=0.8` to ages 6+ (and zero natural mortality to ages 2–5); ReefMod has an eight-year maximum age.
3. **Feeding is linear (Type I).** The fraction of fast coral removed per year is `a_F·A`, independent of how much remains. At `A≈0.85` this removes essentially all fast coral each year, so Lizard model coral stays at 0.035–0.10 from 1992 to 2024. This is addressed by switch 3 below, but that arm is **deferred** (see "Units caveat").

### Hypotheses and equations

Let `F_t`, `S_t` be fast (groups 1–3) and slow (groups 4–5) prey cover as fractions of the site's **habitable area** (the same support as the legacy rule). Adults enter at age 2. `p̃` and `m3` are the archived values.

**H1 — preferred-prey starvation with short memory** (`preferred_prey_starvation`):

`F̄_t = (1-α)F̄_{t-1} + αF_t`, with `α = 1-exp(-1/τ)` and `F̄` initialized to `F` at the first transition.

`f_t = 1` if `F̄_t > θ`; otherwise `f_t = (1-p̃) + p̃(F̄_t/θ)^3`.

This replaces only the argument and threshold of the legacy cubic. `f` still multiplies N2 survival and adult survival, and N1 is still exempt (age-0 COTS eat crustose coralline algae; this is consistent with all three comparators).

**H2 — adult age structure with senescence** (`adult_senescence`):

Adults are tracked in ages `2,…,a_s-1`, plus a senescent plus group `a_s+`. Ages below `a_s` survive at `(1-m3)f`; the plus group survives at `(1-m_s)f`. New adults enter age 2. `N[3]` remains the summed adult density, so fecundity, Allee, consumption, logging and observation mapping are unchanged in form.

Initial 1991 adults are spread over the stable age distribution under background survival with `f=1`. With `s3=0.805` and `m_s=0.8` the weights are about `0.286, 0.230, 0.185, 0.149, 0.150`. This is an **assumption**, not data.

**Expected observable consequences**

| Arm | Expected effect | Risk |
| --- | --- | --- |
| H1 | Steep post-peak decline once fast coral falls below `θ`, despite remaining massives. Deeper 2004–12 adult floor. Possible earlier coral release. | Model fast coral is pinned low from 1992 (switch 3 deferred), so H1 may also suppress the build-up to the 1997 wave. That would be an informative failure. |
| H2 | The first-wave cohort is removed about 4 years after maturation regardless of density or food. Lower adult carryover and a shorter peak. | Peak height may fall because old adults no longer accumulate. External N1 supply continues, so the floor becomes `maturation × ~3–4 yr` residence. |
| H1×H2 | The deepest trough, if the first wave survives. | Losing the first wave. |

### Preregistered settings and evidence

| Parameter | Value | Units | Basis | Status |
| --- | --- | --- | --- | --- |
| `θ` | `0.05` | preferred cover, fraction of habitable area | ReefMod 7.0 `COTS_coral_threshold` (fraction of a 400 m² cell). Habitable area ≤ cell area, so this is at most as strict. | Transferred convention, not a field estimate |
| `τ` | `0.5` | years (gives `α≈0.865`) | Starvation tolerance on the order of months. Near-contemporaneous at an annual step. | Fixed; not a factor |
| `a_s` | `6` | years | CoCoNet 6+ senescent class; Pratchett et al. 2014 senescent phase; ReefMod maximum age 8 | Transferred convention |
| `m_s` | `0.8` | annual mortality fraction | CoCoNet `M=0.8` at 6+, applied there as a rate scaled by `1/C^f` (≥0.93 annual fraction when `C^f≤0.3`). The fraction used here is weaker, hence conservative. | Transferred convention |
| `m3` (ages 2–5) | archived `0.1946` | annual mortality fraction | Archived candidate | Unchanged |

### Factorial

Four arms on the [archived 1991 full-cohort anchor](runs/20261007T199111_owen_1991_cohort_history_final): control (both off, exact replay), H1 only, H2 only, and H1×H2.

Held fixed across all arms:
- one hydrodynamic seed (`20260930`) and `q=0.015` initialization;
- the 6:3:1 preload and the 3 COTS/ha Allee threshold;
- external outbreak production, internal connectivity, fecundity and consumption parameters;
- the frozen per-reef observation scales and the frozen three-year-smoothed score.

Juvenile storage, the lagged hazard, low-cover mortality and size weighting are all **off**. No optimizer.

Record:
- reef-level stages;
- maturation, all settlement sources, fast and total coral;
- the preferred food memory and the realized food survival factor;
- senescent-class density and excess senescent deaths;
- the score table and the peak/trough audit.

### Gate, failures and rollback

The qualitative gate is unchanged:
- a first and a second detected modeled peak on Lizard, MacGillivray and North Direction;
- 7/7 observed peaks matched;
- worst two-wave trough/smaller-peak ratio `≤0.10`;
- no material peak-height or coral degradation against the paired control.

Each of the following is a failure:
- an inexact control replay (trajectories, all 12 flux channels), or any change to external immigration;
- a non-finite or negative state, or adult age classes that do not sum to `N[3]`;
- losing Lizard's 1997 peak;
- a trough that refills immediately from N2;
- an improvement caused only by extinguishing external supply.

A one-seed pass is only a *candidate* for replicated-seed, `q`/cohort-sensitivity and held-out testing.

Rollback is `COTS_PREFERRED_PREY_STARVATION=false` and `COTS_ADULT_SENESCENCE=false`, which are the defaults.

### Units caveat and deferred switch 3

The frozen Lizard observation scale is `0.2942 COTS/tow` per model adult/ha. CoCoNet uses `0.015` (Moran & De'ath 1992), and ReefMod's `2.7 per 400 m²` disease threshold is commented as "~1 per tow", which also implies about `0.015`. The fitted scale therefore implies one model "COTS/ha" ≈ 20 field COTS/ha. As a result, the modeled Lizard peak (`≈2.96`) lies *below* the 3 COTS/ha Allee threshold, and transferred density thresholds (e.g. ReefMod's `≈67 COTS/ha`) would never trigger.

H1 and H2 are expressed in cover fractions and years, so they are insensitive to this ambiguity, which is why they are tested first.

A third opt-in switch, `per_capita_consumption`, is implemented and unit-tested but **not run** here. It removes `A·c/10⁴` cover per year (`c` in m² per adult per year, default `10`, the Keesing & Lucas 1992 order of magnitude) preferred-first, capped by available prey. It interprets `N[3]` strictly as adult COTS per habitable hectare, so it is only meaningful once the observation scale is fixed near `0.015` or otherwise justified. That decision is a separate, explicit change to the frozen observation model and has not been made.

### Decision (2026-10-09)

[Run](runs/20261009T_02_structural_crash_1991_v1/metadata.toml) and [audit](runs/20261009T_02_structural_crash_1991_v1/structural_crash_audit.json): the control replay is exact, and external immigration is identical in every arm. **No arm passes.**

| Arm | Worst trough ratio | Matched peaks | Lizard first peak |
| --- | --- | --- | --- |
| Control | `0.594` | 6/7 | 1997 |
| H1 (preferred starvation) | `0.723` | 6/7 | 1995, 14% of control height |
| H2 (senescence) | `0.473` | 6/7 | 1997, −17% |
| H1×H2 | `0.710` | 6/7 | 1995 |

The H1 risk stated above materialised. Control fast coral is about 0 from 1992 because linear grazing at `a_F=1.18` removes all of it. H1 therefore starves the first wave and then pins COTS at about `0.3 COTS/ha`, holding fast coral just under `θ`. That makes it a test of the archived consumption calibration rather than of H1.

H2 deepens the trough as predicted, but the gap adults are young recruits from continued settlement, so senescent deaths are small there.

Both switches stay experimental and off. Next gate: an explicit, documented observation-unit decision, then a bounded test of preferred-first consumption (`per_capita_consumption` or `h>0`) combined with H1 and H2. The Allee threshold and external supply stay fixed throughout.

## Observation-unit migration to q = 0.015 and consumption factorial (2026-10-09)

Status: **observation-model change approved by the project lead on 2026-10-09; preregistered before implementation.** Earlier runs keep their per-reef frozen scales and are not rescored or overwritten. The new objective is labelled `fixed_q0015`.

### Observation-unit decision

Simulated manta CPUE is `q·A`, where `A` is adult COTS per hectare of site area and `q = 0.015 COTS/tow per COTS/ha` is used **for every reef**. This replaces the per-reef fitted nuisance scales (Lizard `0.2942`, MacGillivray `0.1642`, North Direction `0.1258`, Eyrie `0.0647`).

Basis:
- CoCoNet v3.4 Table A.5 uses `ε = 0.015` ("control program CoTS ha⁻¹ to CoTS per manta tow", citing Moran & De'ath 1992).
- ReefMod 7.0 comments its `2.7 COTS/400 m²` (≈67.5/ha) dieback threshold as "~1 per tow", which implies ≈0.0148.
- Both follow from pairing two outbreak benchmarks, `0.22 COTS/tow` and `15 COTS/ha` (Moran & De'ath 1992): `0.22/15 ≈ 0.0147`.

This is a **shared modelling convention, not a measured conversion**. No source states it as a fitted relationship. Tow-area assumptions differ between sources, and manta-tow detectability is low and variable: Fernandes et al. 1990 suggest under 5% detection; reported detectability across methods ranges from about 23% to 88%. `q` is therefore fixed rather than fitted, and its uncertainty remains an open observation-model limitation.

Consequences:
- Per-reef amplitude differences that the fitted scales used to absorb now count against the model.
- Literature density thresholds (the 3 COTS/ha Allee half-saturation, the 15 COTS/ha outbreak benchmark) become comparable to the model state.

### Parameter-unit migration and compatibility test

The archived candidate was fitted in a compressed density unit. Under Lizard's frozen scale, one model unit equals `k = 0.294195/0.015 ≈ 19.61` COTS/ha at `q`. Lizard is the primary calibration reef, so the migration uses Lizard's `k`. The other reefs' residual factors (`10.9`, `8.4`, `4.3`) are reported as amplitude mismatches, not hidden.

Every density-dimensioned quantity is transformed, which is an exact reparameterisation:

| Quantity | Transform |
| --- | --- |
| Adult, juvenile and recruit state; 1991 initial CSV | `×k` |
| Background immigration `IMM`; external outbreak production | `×k` |
| Ricker `b_ricker`; linear consumption `a_F`, `a_S` | `÷k` |
| Allee threshold | `×k` (compatibility arm only) |
| Per-adult fecundity, survival, cover thresholds, connectivity, settlement | unchanged |

**Compatibility test:**
- **Core:** a COTSMod unit test shows that scaling state and density parameters by `k` reproduces the trajectory `×k` and the same coral consumption.
- **Integration:** the run's compatibility arm (Allee `3k`, initial state `×k`) must reproduce the archived control's adults and flows `×k` to relative tolerance `1e-9`, with identical Lizard simulated CPUE.

The **migrated control** then differs from the archived control in two respects only:

1. **Allee half-saturation.** It is held at the evidence-based `3 COTS/ha` in true units, which is `k` times weaker relative to the dynamics than in the archived run.
2. **1991 initial state.** The `IDW CPUE/q` preload is now consistent with the observation unit. Under the archived scale it was effectively `k` times too large.

The archived 1997 Lizard first-peak match depended on that inconsistency. The migrated control is expected to shift the first wave, likely later from the lower start, although the weaker Allee pushes the other way. This must be reported, not tuned away.

### Consumption factorial (2×2, plus the compatibility arm)

Units: per-capita grazing removes `A·c/(10⁴·h)` of habitable prey cover per year, where:
- `A` is adults per site hectare;
- `c = 10 m²` per adult per year, the Keesing & Lucas 1992 order of magnitude for adult feeding;
- `h` is the site habitable fraction (`location_k`).

Fast prey is taken first, then slow, capped by availability.

| Arm | Consumption | Starvation / adult structure |
| --- | --- | --- |
| `compat_rescaled_exact` | migrated linear | legacy, with Allee `3k` and initial `×k` |
| `migrated_control` | migrated linear | legacy |
| `migrated_structural` | migrated linear | H1 (`θ=0.05`, memory `0.5` yr) + H2 (`a_s=6`, `m_s=0.8`) |
| `per_capita_legacy` | per-capita, `c=10` | legacy |
| `per_capita_structural` | per-capita, `c=10` | H1 + H2 |

H1 and H2 are combined because the previous screen separated their effects. The new factor is the consumption form.

All arms keep the following unchanged: seed `20260930`, Allee `3 COTS/ha` (except the compatibility arm), migrated `IMM` and external production, connectivity, fecundity, frozen three-year-smoothed score with `q=0.015` on every reef, and the reef mapping. There is no optimizer.

**Expected:**
- Per-capita grazing removes about `0.0017·A` habitable cover per year at `h≈0.6`. Fast coral should therefore persist through low and moderate densities and be stripped only at outbreak densities (`A ≳ 30–60/ha`).
- With H1 this should produce a starvation-driven crash *after* the outbreak peak, followed by coral release.
- Risk: with grazing per adult about 10× weaker than the migrated linear fit, outbreaks may no longer limit coral, and peaks may overshoot or fail to crash.

**Gate and failures.** The qualitative gate is unchanged, but amplitudes are now absolute under `q`. Failures include:
- a compatibility mismatch;
- any change to external immigration other than the documented `×k` migration;
- a non-finite state;
- losing all first waves;
- an improvement achieved only by extinguishing external supply.

A one-seed pass remains a candidate pending replicated seeds, `q` sensitivity (e.g. `0.0075`–`0.03`), held-out years/reefs, and a refit under the new unit.

**Rollback.** The per-reef frozen-scale objective and archived parameters remain valid for legacy comparisons. All switches default off.

### Unit-migration decision (2026-10-09)

**Deviation from preregistration (disclosed).** v1 failed its compatibility arm. The cause is ADRIA's parameter-derived coral RNG seed (`floor(sum(parameter values))`), not the migration. A one-ulp replay, a k=1 identity and a k=2 power-of-two rescale isolated it. v2 pins `ADRIA_RNG_SEED` to the archived control's derived value (`90972102`) for every migration arm. That keeps the archived control's random stream and makes the comparison fair. No other setting changed.

[v2 run](runs/20261009T_07_unit_migration_q015_1991_v2/metadata.toml) and [audit](runs/20261009T_07_unit_migration_q015_1991_v2/unit_migration_audit.json): the compatibility replay is exact (max relative delta `3e-15`); external immigration is identical in every arm; simulated CPUE equals `q·A`. **No arm passes.**

| Arm | Two-wave reefs with 2 peaks | Worst trough | Matched | Lizard first peak |
| --- | --- | --- | --- | --- |
| Compatibility (archived ×k) | 3/3 | `0.594` | 6/7 | 1997, 0.85 COTS/tow |
| Migrated control | 3/3 | `0.381` | 6/7 | 2000, 1.07 |
| Migrated linear + H1 + H2 | 3/3 | `0.801` | 6/7 | 2000, 0.24 |
| Per-capita | 0/3 | — | 2/7 | 2003, 2.62 |
| Per-capita + H1 + H2 | 1/3 | `0.171` | 3/7 | 2001, 1.37 |

Per-capita grazing with H1 and H2 reproduces the expected sequence:
1. Fast coral persists while COTS build.
2. The outbreak strips it.
3. Remembered preferred food collapses (food survival 0.06–0.12) and senescent deaths peak (2.6 COTS/ha/yr in 2002).
4. Adults fall from 120 to about 3 COTS/ha.
5. Coral recovers.

Its failures are the absent second waves at MacGillivray and North Direction and first peaks around 2001. The migrated control's lateness (2000 versus 1997) shows that the archived timing match depended on the inconsistent preload.

None is promotable. The per-capita structural configuration is the first candidate that passes the *qualitative crash* part of the mechanism gate on Lizard. The next bounded diagnostic should decompose second-wave supply on the reefs that fail to reseed, holding the Allee threshold and external supply fixed, before any refit under the new unit.

## Preregistered second-wave factorial: prey switching × preferred-food condition (2026-10-12)

Status: **preregistered before implementation; no default or objective change.** The base model is the project lead's preferred configuration, `per_capita_structural` from [unit migration v2](runs/20261009T_07_unit_migration_q015_1991_v2/metadata.toml):
- per-capita grazing at `10 m²/adult/yr`, preferred-first, habitable-area corrected;
- H1 preferred-prey starvation (`θ=0.05`, memory `0.5` yr);
- H2 senescence (`a_s=6`, `m_s=0.8`);
- `q=0.015` on every reef, migrated parameters (`k=19.61`), Allee `3 COTS/ha`;
- `ADRIA_RNG_SEED=90972102`.

### Diagnosis motivating the test

The [second-wave readout](runs/20261009T_08_second_wave_audit_1991_v1/metadata.toml) shows that settlement during the second build-up equals wave 1, so supply is not limiting. Two problems are:

1. **Fast coral is driven almost extinct.** Exhaustive preferred-first grazing takes fast cover to `0.001–0.005`. It regrows at only about `+0.01/yr`, reaching about `0.035` before wave 2 against `0.13–0.14` before wave 1, so wave 2 starves early.
2. **Starving adults stay fertile.** Body condition tracks *total* cover over `tau_condition=8.3` yr, so starving adults keep breeding. 63% of local larval production comes from sites with food survival `<0.5`. The fertile residual population keeps grazing fast coral.

### Hypotheses and equations

Fast and slow cover `F`, `S` are fractions of habitable area. `D = A·c/(10⁴·h)` is the per-capita grazing demand.

**S1 — prey switching** (`prey_switching`, requires `per_capita_consumption`). Following MICE (Rogers & Plagányi 2022), the share of demand directed at preferred prey falls as preferred cover is depleted:

`ρ = 1 / (1 + exp(-s·(F/K - F_mid)))`, with `K = 1` (cover is already relative to habitable area).

`Cons_F = min(ρ·D, F)`, `Cons_S = min(D - Cons_F, S)`.

Unmet preferred demand still falls on slow prey. Preregistered: `F_mid = 0.10`, `s = 70`, taken directly from the MICE prey-switching logistic.

*Expected:* remnant fast coral survives the first wave (no near-extinction), so fast coral recovers sooner and to a higher level, giving a larger, earlier second wave.

*Risk:* weaker coral loss in wave 1. H1 starvation still keys on fast cover, so COTS can feed on massives yet still be recorded as starving.

**S2 — preferred-food condition** (`preferred_condition`). Body condition (the multiplier on fecundity via `bc²`) is driven by fast cover rather than total cover:

`signal = min(1, F/F_c)`, `bc_t = (1-α)·bc_{t-1} + α·signal`, with `α = 1 - exp(-1/τ_c)`.

Preregistered: `F_c = 0.10` (the same preferred cover at which diet becomes preferred-dominated) and `τ_c = 1` yr. Supported by Caballes et al. maternal-diet trials: Acropora-fed females had heavier gonads and larger oocytes than Porites-fed or starved females (60-day trials, three females per treatment, so qualitative support only).

*Expected:* residual starving adults stop producing larvae, internal immigration falls in troughs, and troughs deepen. With less residual grazing, fast coral recovers sooner.

*Risk:* the second wave weakens unless external supply and recovered condition carry it, and the first wave builds more slowly.

### Factorial

2×2 on the base:

| Arm | Prey switching | Preferred-food condition |
| --- | --- | --- |
| `base` | off | off |
| `switching` | on | off |
| `condition` | off | on |
| `switching_condition` | on | on |

`base` must replay v2 `per_capita_structural` exactly (trajectories, flows and mechanisms), and external immigration must be identical in every arm.

Held fixed: seed, `q`, `k`, Allee, external production, connectivity, fecundity scale, `c=10`, H1/H2 settings, frozen three-year-smoothed score. No optimizer.

Record (in addition to the existing diagnostics): reef-level body condition and the fast-diet share `ρ`.

### Gate, failures and rollback

The qualitative gate is unchanged and in absolute CPUE under `q`:
- two detected peaks on Lizard, MacGillivray and North Direction;
- 7/7 observed peaks matched;
- worst trough ratio `≤0.10`;
- no material coral degradation against `base`.

Also report second-wave height against the observed 2013–15 peaks (0.67–0.85 COTS/tow), fast coral at second-wave onset, and trough adult density.

Failures:
- inexact base replay, or changed external supply;
- non-finite states;
- losing all first waves;
- improvement only through extinguishing external supply.

A one-seed pass is a candidate only, pending replicated seeds, `q` sensitivity, held-out tests and a refit.

Rollback: both switches default off.

Known confounders deferred to the next step:
- the uniform `0.5` habitable fraction on all target sites;
- zeroed cyclone and wave forcing.

### Switching × condition decision (2026-10-12)

[Run](runs/20261012T_01_switching_condition_1991_v1/metadata.toml) and [audit](runs/20261012T_01_switching_condition_1991_v1/switching_condition_audit.json): base replays v2 `per_capita_structural` exactly, and external immigration is identical across arms. **No arm passes.**

| Arm | Two-wave reefs | Worst trough | Matched | Mean 2010–17 max CPUE | Mean fast coral at onset | Mean gap adults/ha |
| --- | --- | --- | --- | --- | --- | --- |
| Base | 1/3 | `0.171` | 3/7 | 0.33 | 0.007 | 12.6 |
| Switching | 0/3 | — | 2/7 | 2.22 | 0.019 | 90.0 |
| Condition | 2/3 | `0.366` | 5/7 | 0.38 | 0.027 | 7.2 |
| Both | 2/3 | `0.306` | 3/7 | 2.67 | 0.019 | 73.9 |

**S2 (preferred-food condition)** behaved as predicted:
- starving adults lose fecundity (gap body condition `0.17–0.37` vs `0.63–0.97` in base);
- trough internal immigration falls about 70%;
- trough adults fall;
- fast coral recovers sooner;
- Lizard produces a trough of `0.103` and a 2014 second wave.

Its risk (weaker second waves) did not materialise on Lizard. On the other reefs second waves exist but are small (`0.14–0.22 COTS/tow` vs `0.71–0.85` observed), because their fast coral stays at 0–0.03.

**S1 (prey switching)** failed through the stated risk, amplified by two design facts:
- switching moves the full per-capita demand onto massives and treats them as full food, so slow coral is eaten to near zero;
- the switch midpoint (`0.10`) lies above the starvation threshold (`0.05`), so the protected fast-coral band keeps remembered food above threshold.

The result is persistent high-density populations and recurrent outbreaks, not boom–bust.

**Decision.** Keep S2 as the working candidate (opt-in, not default) on the per-capita structural base. Retain S1 as failed. In MICE, food limitation is driven by the *same* prey-switching function, so a coherent revisit would tie starvation to diet share or discount non-preferred nutrition. That requires a new preregistration, not a parameter retune.

**Next gate.** Fast-coral recovery on MacGillivray, North Direction and Eyrie: test the uniform `0.5` habitable fraction against site data, and the zeroed cyclone and wave forcing, before any refit. First-wave timing (around 2001 vs observed 1996–98) also remains open.

## Preregistered prey-switching repair: diet-linked starvation × non-preferred feeding efficiency (2026-10-12)

Status: **preregistered before implementation; no default or objective change.** The [switching × condition screen](runs/20261012T_01_switching_condition_1991_v1/metadata.toml) failed for two identified reasons:
1. Diet switching treated massives as full food at the full per-capita rate, so slow coral was eaten to near zero.
2. The switch midpoint (`0.10`) sat above the starvation threshold (`0.05`), so the protected fast-coral band kept COTS fed.

Both fixes below follow the published comparators and are not retunes of the failed parameters.

### Hypotheses and equations

Notation is as before: `F̄` is remembered fast cover (H1 memory `0.5` yr) and `ρ(F) = 1/(1+exp(-70(F-0.10)))` is the MICE switching logistic.

**A — diet-linked starvation** (`diet_linked_starvation`; requires `prey_switching` and `preferred_prey_starvation`). Food survival uses the same logistic that drives the diet switch, evaluated on remembered fast cover:

`f = (1 - p̃) + p̃·ρ(F̄)`, replacing the cubic below `θ`.

This is the MICE structure, where food limitation of 1+ COTS is a function of the same preferred-prey depletion term (Rogers & Plagányi 2022). The H1 threshold `θ` is not used when A is on.

*Expected:* starvation begins as soon as the diet shifts to non-preferred prey, so high-density persistence on massives ends and boom–bust returns while switching still protects remnant fast coral.

*Risk:* starvation now starts at fast cover around `0.10`, higher than `θ=0.05`. Reefs with low initial fast cover (Eyrie about `0.05`, North Direction about `0.10`) may lose their first wave.

**B — reduced feeding on non-preferred prey** (`nonpreferred_feeding_efficiency = e`, applies to per-capita grazing). Demand not met from fast prey is taken from slow prey at reduced efficiency:

`Cons_S = min(e·(D - Cons_F), S)`.

Preregistered `e = 0.2`: the order of magnitude of ReefMod 7.0's feeding-preference ratio for massive to Acropora groups (about `0.055/0.257 ≈ 0.21`). This is a **transferred convention**, not a field feeding rate. The default `e = 1` preserves the current behaviour.

*Expected:* massive coral is no longer stripped, and coral totals stay close to the condition arm.

*Risk:* B alone does not stop COTS persisting while starvation still keys on `θ` below the switch band.

### Factorial

All arms keep:
- per-capita grazing (`c = 10`), H1 memory, H2 senescence, preferred-food condition (saturation `0.10`, memory `1` yr);
- `q = 0.015`, `k = 19.61`, Allee `3`, external production, connectivity, frozen score;
- `ADRIA_RNG_SEED = 90972102`.

| Arm | Prey switching | A | B (`e`) | Purpose |
| --- | --- | --- | --- | --- |
| `condition_ref` | off | off | 1 | must replay the 2026-10-12 `condition` arm exactly |
| `switching_ref` | on | off | 1 | must replay the 2026-10-12 `switching_condition` arm exactly |
| `switch_diet_starvation` | on | on | 1 | A only |
| `switch_low_efficiency` | on | off | 0.2 | B only |
| `switch_diet_starvation_low_efficiency` | on | on | 0.2 | A + B |

### Gate, failures and rollback

The qualitative gate is unchanged:
- two peaks on Lizard, MacGillivray and North Direction;
- 7/7 matched;
- worst trough ratio `≤0.10`;
- no material coral degradation.

In addition, report against `condition_ref`:
- massive (slow) coral retained;
- trough adult density;
- second-wave height.

Failures:
- inexact reference replays, or changed external supply;
- non-finite states;
- losing all first waves;
- improvement only by extinguishing external supply.

Rollback: all switches default off, with `e = 1`.

### Prey-switching repair decision (2026-10-12)

[Run](runs/20261012T_02_switching_repair_1991_v1/metadata.toml) and [audit](runs/20261012T_02_switching_repair_1991_v1/switching_repair_audit.json): both references replay the previous run exactly, and external supply is identical. **No arm passes.**

| Arm | Two-wave reefs | Worst trough | Matched | Gap adults/ha | Gap slow coral | 2012 coral |
| --- | --- | --- | --- | --- | --- | --- |
| condition_ref | 2/3 | `0.366` | 5/7 | 7.2 | 0.237 | 0.347 |
| switching_ref | 2/3 | `0.306` | 3/7 | 73.9 | 0.001 | 0.030 |
| A diet-linked starvation | 2/3 | `0.261` | 1/7 | 37.8 | 0.087 | 0.076 |
| B efficiency 0.2 | 2/3 | `0.121` | 1/7 | 159.2 | 0.023 | 0.004 |
| A + B | 1/3 | `0.686` | 3/7 | 33.2 | 0.285 | 0.374 |

Against the predictions:
- **A ended the full-food persistence on massives,** but starvation along the switching curve is partial: with fast coral held near `0.07`, survival is about `0.4`.
- **B protected massives only in combination with A.** Alone it let starving-but-surviving COTS reach 4–5 COTS/tow.
- **A + B achieved the intended coral outcome but not boom–bust.** Fast coral is held in the switching band, starfish are partly starved, and condition stays high because condition saturates at fast cover `0.10`. The result is about 30 adults/ha that persist and seed one another.

**Decision.** Prey switching is retired for the outbreak-dynamics model. A prey refuge stabilises the consumer–resource interaction, which is the opposite of the observed boom–bust. Retuning the midpoint, steepness or efficiency to restore cycles would be fitting rather than mechanism. Both fixes remain opt-in in COTSMod for future comparison. The working candidate stays per-capita grazing + H1 + H2 + preferred-food condition.

**Next gate (unchanged).** Fast-coral recovery and second waves on MacGillivray, North Direction and Eyrie via the habitable fraction and disturbance forcing; first-wave timing.

## Preregistered habitat × bleaching screen on the condition model (2026-10-12)

Status: **preregistered before implementation; no default, domain-package or objective change.**

The base is the working candidate: per-capita grazing (`c=10`), H1 starvation, H2 senescence and preferred-food condition, **without** prey switching. It runs at `q=0.015`, `k_unit=19.61`, Allee `3`, with `ADRIA_RNG_SEED=90972102`. The base must replay `condition_ref` in [switching repair](runs/20261012T_02_switching_repair_1991_v1/metadata.toml) exactly.

### Motivation

Second waves on MacGillivray, North Direction and Eyrie are small because fast coral stays at 0–0.03 while total coral recovers. Two coral-side assumptions are untested.

1. **Habitable fraction.** `build_lizard_domain_v2.py` hard-codes `k=0.5` for every site and records this as a known limitation. Yet the V2 initial cover is stored as a fraction of each reef's *RME* habitable area, computed as `1 - sand_proportion` from `rme_ml_2025_06_05/data_files/id/id_list_2024_12_01.csv`. With `k=0.5` the absolute initial coral is understated by `0.5/k_RME` (about 0.7×). Per-capita grazing is also converted with `k`.
   - Site values: `k_RME` mean `0.715`, range `0.605–0.805`, all 2,895 sites matched.
   - Target reefs: Lizard `0.722–0.773`, MacGillivray `0.751`, North Direction `0.781`, Eyrie `0.605`.
2. **Disturbances.** Bleaching (historical DHW) is active. Cyclone and wave forcing are zero because the V2 package has no cyclone or wave data. The local RME cyclone files are 100 *hypothetical* sequences with random start years, so no historical cyclone hindcast is available locally. Cyclones therefore stay out of this screen and are recorded as a data gap. Bleaching is tested as a **labelled counterfactual** only, to quantify how much it suppresses fast-coral recovery. Removing it is not a candidate model.

### Hypotheses

- **H-k.** Using `k_RME` makes coral cover consistent with its source and enlarges absolute habitat. This raises absolute coral and per-site coral growth capacity, and lowers per-capita grazing in cover terms (÷`k`). *Expected:* faster fast-coral recovery and larger second waves. *Risk:* larger first outbreaks and shifted timing.
- **H-b (counterfactual).** Bleaching mortality in 1998, 2002, 2016 and 2017 suppresses fast-coral recovery. *Expected:* without bleaching, fast coral at wave-2 onset is higher and second waves are larger. The difference measures bleaching's share of the second-wave deficit.

### Implementation and factorial

- `k_RME` is applied at run time from the checked-in RME id list via each site's `GBRMPA_ID`, with the source hash recorded. The domain package is not rebuilt.
- The initial relative cover is unchanged, which is the consistent interpretation of the stored values.
- Bleaching is removed with ADRIA's existing `ADRIA_DEBUG_BLEACHING_SCALAR=0`.

| Arm | `k` | Bleaching |
| --- | --- | --- |
| `k05_bleaching` | 0.5 (must replay `condition_ref`) | on |
| `krme_bleaching` | RME | on |
| `k05_no_bleaching` | 0.5 | off (counterfactual) |
| `krme_no_bleaching` | RME | off (counterfactual) |

Held fixed: all COTS settings, seed, `q`, external supply, connectivity, frozen score, and zero cyclone and wave forcing.

Report for each arm:
- the qualitative gate;
- second-wave height against observed;
- fast coral at wave-2 onset;
- 2012 total coral against LTMP 9 m (an approximate support comparison: model cover is per habitable area, manta is per reef bottom).

### Gate, failures and rollback

The gate is unchanged. A no-bleaching arm cannot be promoted under any result.

Failures:
- inexact base replay, or changed external supply;
- non-finite states;
- loss of all first waves.

Rollback: `k=0.5` and bleaching on (the defaults).

If `k_RME` improves second waves without degrading first waves, the next step is a checked-in domain rebuild (`k_RME`, with provenance), not a run-time override.

### Habitat × bleaching decision (2026-10-12)

[v2 run](runs/20261012T_04_habitat_disturbance_1991_v2/metadata.toml) and [audit](runs/20261012T_04_habitat_disturbance_1991_v2/habitat_disturbance_audit.json): the base replays `condition_ref` exactly and external supply is identical. [v1](runs/20261012T_03_habitat_disturbance_1991_v1) completed all arms, but the runner failed while writing metadata because a top-level variable named `values` shadowed `Base.values`. v1 is invalid and retained; v2 is the identical design after the rename. **No arm passes.**

| Arm | Two-wave reefs | Worst trough | Matched | Gap adults/ha | Fast coral at onset | Mean 2012 coral |
| --- | --- | --- | --- | --- | --- | --- |
| `k05_bleaching` (base) | 2/3 | `0.366` | 5/7 | 7.2 | 0.027 | 0.347 |
| `krme_bleaching` | 1/3 | `0.128` | 3/7 | 15.8 | 0.017 | 0.301 |
| `k05_no_bleaching` (counterfactual) | 1/3 | `0.211` | 3/7 | 12.3 | 0.022 | 0.293 |
| `krme_no_bleaching` (counterfactual) | 1/3 | `0.118` | 3/7 | 31.2 | 0.014 | 0.188 |

Both expectations were wrong in direction:
- **H-k:** more habitat lowered per-capita grazing in cover terms and raised coral capacity, but this grew larger first outbreaks (Lizard 1.93 COTS/tow, about the observed 1.92, but in 2002) and doubled residual adults. Fast coral at onset fell.
- **H-b:** without bleaching, more coral again fed larger outbreaks and more residual adults. Bleaching is not the cause of the second-wave deficit.

**Decision.**
- Nothing is promoted.
- The weak second waves are controlled by the **trough adult density** (model 7–31/ha vs observed about 0–3/ha), which keeps fast coral below the level that can support a second outbreak.
- RME-k is the data-consistent habitable fraction. Adopting it through a checked-in domain rebuild is a data-correction decision for the project lead, to be taken before refitting.
- Historical cyclone forcing remains a data gap (no hindcast locally).

**Next gate.** Decompose what sustains trough adults on the condition model (external boundary supply, internal settlement, juvenile/adult survival) and address first-wave timing (around 2001–02 vs observed 1996–98).

## Preregistered benthic habitable fraction × historical cyclone screen (2026-10-12)

Status: **preregistered before implementation; no default, domain-package or objective change.** The base is the working condition model (no switching), as in the habitat screen. The base arm must replay `k05_bleaching` in [habitat v2](runs/20261012T_04_habitat_disturbance_1991_v2/metadata.toml) exactly. Bleaching is on in all arms.

### New derived inputs

Both are checked-in scripts with provenance. Their outputs live in the git-ignored `sandbox/data/Lizard_Historical_v0.2/derived/` because the sources are internal or licensed.

1. **Site habitable fraction** (`sandbox/domain_building/build_site_habitable_fraction.py`), from the GBR10 2025 benthic map (10 m; classes Sand 11, Rubble 12, Rock 13, Coral/Algae 15).
   - `k_total` = (Rock + Coral/Algae) / all cells in the site polygon. Mean `0.435`; 259 sites are zero.
   - `k_mapped` = (Rock + Coral/Algae) / *mapped* cells. It assumes unmapped area (mean 32% of sites, typically beyond the optically mapped depth) has the mapped composition. Sites that are less than 20% mapped use their reef's mapped composition (390 sites), then the RME value (2 sites). Mean `0.625`.
   - Target reefs (`k_total` / `k_mapped`): Lizard 0.62–0.64 / 0.71–0.97, MacGillivray 0.67 / 0.92, North Direction 0.80 / 1.0, Eyrie 0.64 / 0.74.
   - Sites with `k = 0` become non-habitable for coral in ADRIA.
2. **Historical cyclone exposure** (`sandbox/domain_building/build_cyclone_hours_rrn.py`), from the GBRMPA RRN 2025 workbook variable `cyc_maxHrs4mw`: maximum hours of ≥4 m significant waves in any one cyclone per summer, per reef. The 4 m damage threshold follows Puotinen et al. 2016.
   - All 113 domain reefs are matched (1985–2025); each summer is assigned to its ending year.
   - Within 1991–2024 the material events are 1991 (up to 35 h) and Ita 2014 (13–16 h) and Nathan 2015 (49–74 h).

### Cyclone mortality conversion

Annual coral survival from cyclone exposure is `exp(b·H)`, where `H` is hours of ≥4 m waves and `b = -0.00444 h⁻¹` (90% CrI −0.0058 to −0.0030). This is the posterior median of the LTMP dose-response model in COTS-Coral-OutcomesV2 (`q3_doseresponse_coefs.csv`: −0.0444 per 10 h on the annual log ratio of hard-coral cover).

It is applied uniformly to all functional groups and size classes, because the coefficient is estimated on total hard coral. Implied losses: 15 h ≈ 6%, 60 h ≈ 23%.

Group-specific vulnerability (e.g. CoCoNet's higher mortality for branching and tabular corals) is a deferred sensitivity.

### Hypotheses

- **H-benthic.** Site-resolved habitable fraction redistributes coral capacity. `k_total` lowers habitat on sand-dominated and unmapped sites, so total coral and outbreak food drop; this should give smaller first outbreaks and fewer trough adults. `k_mapped` raises habitat on most reef sites, with effects like the RME arm (larger outbreaks). *Expected:* `k_total` deepens troughs relative to base; `k_mapped` resembles RME-k (worse second waves).
- **H-cyc.** Historical cyclones mainly remove coral in 2014–15, during and after the observed second wave. *Expected:* lower 2015–24 coral and smaller or truncated late outbreaks; little effect on the 2004–12 trough.

### Factorial (3 × 2)

| Arm | `k` | Cyclones |
| --- | --- | --- |
| `k05_nocyc` | 0.5 (replay) | off |
| `k05_cyc` | 0.5 | RRN |
| `ktotal_nocyc` | benthic total | off |
| `ktotal_cyc` | benthic total | RRN |
| `kmapped_nocyc` | benthic mapped | off |
| `kmapped_cyc` | benthic mapped | RRN |

Held fixed: all COTS settings, seed, `q`, external supply, connectivity, frozen score, waves zero, bleaching on.

Report for each arm:
- the gate;
- second-wave height;
- trough adults/ha;
- fast coral at wave-2 onset;
- 2012 and 2016 coral against LTMP 9 m.

### Gate, failures and rollback

The gate is unchanged. Failures:
- inexact base replay;
- changed external supply;
- non-finite states;
- cyclone mortality outside [0, 1];
- unmatched sites.

Rollback: `k=0.5` and cyclones off (the defaults).

Adopting a benthic `k`, or cyclones, as the working model requires a separate checked-in domain-package rebuild and a decision by the project lead.

**Disclosed deviation (benthic arms, before results).** [v1](runs/20261012T_05_benthic_cyclone_1991_v1) ran both `k=0.5` arms. All four benthic arms failed with `DimensionMismatch`.

The cause is an upstream ADRIA indexing limitation. `scenario.jl` reshapes the habitable-site areas to `length(habitable_locs)`, which is the length of the BitVector, not the count of habitable sites. `constrain_recruitment!` also divides all-site cover by habitable-site areas. ADRIA therefore cannot run any site with `k = 0`.

No benthic result was seen. For v2, benthic `k` is floored at `0.01`, and sites whose benthic `k` is 0 keep zero initial cover. Coral is additionally logged as a fraction of site area (`relative × k`), so that low-`k` sites with high relative cover do not distort reef-level comparisons.

The upstream fix belongs in ADRIA core and is not made here.

### Benthic × cyclone decision (2026-10-12)

[v2](runs/20261012T_06_benthic_cyclone_1991_v2/metadata.toml) and [audit](runs/20261012T_06_benthic_cyclone_1991_v2/benthic_cyclone_audit.json): base replay is exact. **No arm passes.**

- Benthic `k` (total or mapped) raises site-area coral towards LTMP 9 m values, but raises trough adults to 11–17/ha and removes all second waves.
- Cyclones (Ita, Nathan) lower 2015–16 coral towards observations without changing COTS metrics.

Expectations held in direction for `k_total` coral, but not for troughs: `k_total` did not deepen them, because the larger habitat on reef-flat and slope sites outweighed the losses on sand sites.

**Decision.** Keep habitat and cyclone forcing as realism inputs to carry forward; adoption requires a domain rebuild by decision of the project lead. Do not expect habitat to fix troughs. Proceed to trough mechanisms on two bases: `k=0.5` (replayable) and geomorphic-infilled `k` with cyclones (realistic).

## Preregistered trough-mechanism screen: juvenile predator pit, preferred-prey maturation delay, gradual Allee fertilisation (2026-10-12)

Status: **preregistered before implementation; no default or objective change.**

### Motivation

The [decomposition](MODEL_LOG.md) shows that trough adults on the working model (4–10/ha, observed about 0–3/ha) are mostly new each year. Settlement of 15–60/ha/yr, 70–87% from other domain reefs, passes a coral-independent recruit stage (N1→N2 survival fixed at 0.59). Trough densities of 4–10/ha leave the 3/ha Allee term weak, so the domain recirculates. Realistic habitat makes this worse. The project lead asked for mechanisms that keep young adults down at low density.

### Evidence

- **Juvenile predation.** Keesing, Halford & Hall 2018 (*MEPS* 597:179) measured juvenile mortality of 2.6%/day at 3 mm, 73% of it from mobile predators, and 0.82%/day at 13 mm. That is about 95–99.99%/yr. Desbiens et al. 2023 (*Coral Reefs*) found 26 rubble-dwelling predators; *Schizophrys aspera* consumed more than 5 juveniles per day. A 2025 PNAS gut-DNA study found decapod consumption of juveniles 1.6–3× lower on reefs prone to outbreaks. CoCoNet represents benthic-invertebrate predation on age-1 COTS as a saturating term, `q = p·S·B/(S + p·B)`.
- **Delayed maturation.** Deaker et al. 2020 (*Biol Lett* 16:20190849) kept juveniles herbivorous on coralline algae for up to 6.5 years without coral; they stayed at 16–18 mm and grew normally once coral was offered ("Peter Pan" / "hidden army"). Wilmes et al. 2020 (*Proc R Soc B*) found the switch to coral is essential to reach maturity. CoCoNet makes the age-1→2 transition proportional to `(C^f)^0.5` (preferred coral cover) and lets non-transitioning age-1 COTS remain in that class.
- **Fertilisation versus spawner distance.** Fertilisation is up to 95% when spawners are adjacent, about 70% at 8 m and over 20% at 60 m (Babcock & Mundy 1992, as cited in later reviews). Models give 70–100% at 10 m (Benzie et al. 1994). For a random spatial pattern, the mean nearest-neighbour distance is about 29 m at 3/ha, 9 m at 30/ha and 5 m at 100/ha. The current Hill form `A²/(3²+A²)` gives 0.99 at 30/ha, whereas the distance data imply roughly 0.7. CoCoNet's form `exp(-(S_allee/A)^0.5)` with `S_allee = 3` (citing Rogers et al. 2017) gives 0.73 at 30/ha.
- **Juvenile background mortality (`m1`).** Field mortality far exceeds the model's `m1 = 0.41`/yr. However, N1 is an effective post-settlement density that enters linearly and is confounded with the fitted fecundity and external production. Raising `m1` alone would scale every wave down uniformly. The empirical evidence is therefore represented by the density-dependent predator pit, not by changing `m1`.

### Mechanisms (opt-in COTSMod switches)

Let `J = N2·(1 - m2)` be juvenile (age-1) survivors of background mortality, in COTS/ha.

**P — juvenile predator pit** (`juvenile_predation`). Predators remove `Q·J/(J + Q)` per year, a CoCoNet-form saturating response with capacity `Q` (COTS/ha/yr). The removed fraction `Q/(J + Q)` approaches 1 when `J ≪ Q` and falls as `Q/J` when juveniles swamp predators.

Preregistered `Q = 20 COTS/ha/yr`. This is a **phenomenological capacity** on the model's effective juvenile scale. Field consumption rates exist, but their mapping to effective model densities is unidentified. Q is chosen to sit in the trough range of `J` (about 10–35) and below outbreak-phase `J` (about 60–250); it is not fitted to peaks.

**D — preferred-prey maturation delay** (`preferred_maturation`). The fraction maturing is `g = (F̄)^0.5`, the CoCoNet form, where `F̄` is remembered fast cover (H1 memory), clamped to [0, 1]. Matured juveniles become corallivorous and face the food factor `f`. Non-maturing juveniles stay herbivorous in N2 and are **not** subject to coral starvation, only to `m2` and P.

- matured adults: `g·J·f`;
- retained juveniles: `(1-g)·J`.

This replaces the legacy total-cover storage path when on.

**A — gradual Allee fertilisation** (`allee_stretched`). The fertilisation factor is `exp(-(3/A)^0.5)` instead of `A²/(3² + A²)`. The scale stays at the evidence-based 3 COTS/ha; only the *shape* changes, following CoCoNet and the distance–fertilisation data. This is not a threshold change.

### Expected consequences

- **P:** deeper troughs, because trough juveniles are eaten while outbreak cohorts swamp predators. Outbreak peaks are slightly reduced. *Risk:* the first build-up is slower (later first peaks).
- **D:** juveniles accumulate in troughs as a herbivorous bank. Alone, this may raise later adult pulses and larger second waves (hidden army). With P, the bank is eroded. *Risk:* first waves are delayed, because `g ≈ 0.35` at initial fast cover of about 0.12.
- **A:** fertilisation drops 10–30% at 10–100/ha. This reduces internal settlement and deepens troughs. *Risk:* slower build-up and later first peaks.
- **P + D + A:** the deepest troughs with a preserved second wave, if the bank is released when coral recovers.

### Factorial (two bases × five treatments)

| Base | Description | Reference |
| --- | --- | --- |
| `b05` | `k = 0.5`, no cyclones | `b05_none` must replay `k05_nocyc` in [benthic v2](runs/20261012T_06_benthic_cyclone_1991_v2/metadata.toml) exactly |
| `bgeo` | geomorphic-infilled `k` (floor 0.01, zero initial cover where benthic and geomorphic give `k = 0`), RRN cyclones on | none (new realistic base) |

Treatments per base: `none`, `P`, `D`, `A`, `PDA` (10 arms).

Held fixed: working-model COTS settings, `q`, `k_unit`, Allee scale `3`, RNG seed, external supply, connectivity, frozen score, waves zero, bleaching on.

Report:
- the gate;
- trough adults/ha against the observed about 0–3/ha;
- second-wave height;
- first-peak timing;
- stored-juvenile bank and predation removals;
- site-area coral against LTMP.

### Gate, failures and rollback

The gate is unchanged. Failures:
- inexact `b05_none` replay, or changed external supply;
- non-finite or negative states;
- losing all first waves;
- an improvement achieved only by extinguishing external supply.

Rollback: all switches default off.

`Q` and the D exponent are not fitted. If P helps, its capacity is the first parameter to bound with field predation data before any refit.

### Trough-mechanism decision (2026-10-12)

[Run](runs/20261012T_07_trough_mechanisms_1991_v1/metadata.toml) and [audit](runs/20261012T_07_trough_mechanisms_1991_v1/trough_mechanisms_audit.json): `b05_none` replays exactly and external supply is identical. **No arm passes.**

Against the predictions:
- **P (predator pit)** removed 9–15 juveniles/ha/yr in the gap. It did not deepen troughs. Instead it suppressed the low-density build-up from 1991, so first peaks moved to 2004–07 and fell inside the "gap" window.
- **D (maturation delay)** created a herbivorous juvenile bank of 130–200/ha. Fast coral never recovered enough (`sqrt(F̄) ≈ 0.2`) to release it, so peaks were weaker and there were no second waves.
- **A (gradual Allee)** kept first-peak timing and reduced peaks about 20%. It deepened troughs modestly on Eyrie and North Direction but removed second waves on MacGillivray and North Direction.
- **PDA** gave the lowest trough minima (North Direction 1.3/ha) but flattened outbreaks.

On `bgeo`, model site-area coral tracks LTMP well at MacGillivray, North Direction and Eyrie.

**Decision.** Nothing is promoted. All three switches stay opt-in. The result identifies a structural coupling: low-density depensation is only compatible with the observed rapid 1994–98 rise if initiation is driven by a supply pulse large enough to swamp it. The model currently initiates from low internal densities with about 7 settlers/ha/yr in 1992–93.

**Next gate.** Preregister a first-wave initiation test: the external boundary pulse magnitude and timing in the early 1990s, and/or 1991 initial stock corrected for manta detectability. Cross it with P and A on the `bgeo` base. A large `Q` or a strong maturation delay must not be fitted to compensate.

## Adoption of v0.3 domain and preregistered first-wave initiation screen (2026-10-12)

Status: **domain adoption approved by the project lead (2026-10-12); initiation screen preregistered before implementation.** Not a promotion to the ADRIA default.

### v0.3 domain (adopted working base)

`sandbox/domain_building/build_lizard_domain_v3.py` builds `Lizard_Historical_v0.3` from v0.2 (git-ignored; internal and licensed sources). It changes three things:

- **Site `k`.** Benthic/geomorphic-infilled habitable fraction (mean `0.580`), floored at `0.01`. 20 sites have zero habitat; their raw value is kept in `k_raw`.
- **Initial cover.** Zero at those 20 sites.
- **Cyclones.** `cyclones/cyclone_mortality.csv` holds site × year mortality `1 - exp(-0.00444·h)` from RRN `cyc_maxHrs4mw`, uniform across groups.

The Lizard loader reads cyclone mortality **only** when `datapackage.json` declares `cyclone_mortality`, so v0.1/v0.2 are unchanged. A loader unit test covers alignment and validation.

**Compatibility test.** The v0.3 working-model arm must replay `bgeo_none` of the [trough run](runs/20261012T_07_trough_mechanisms_1991_v1/metadata.toml) exactly. That arm applied the same `k`, the same zero-cover rule and the same cyclone conversion at run time.

### Initiation hypothesis

The observed first wave rose about ×2.5/yr (Lizard 0.23 → 1.92 COTS/tow over 1994–98); the model rises about ×1.5/yr. Model 1992–93 settlement is about 7/ha/yr against roughly 35+/ha/yr needed. Low-density depensation (predator pit P, gradual Allee A) is plausible but delays initiation unless initiation is supply-driven ("predator pit + pulse"). Two candidate initiation sources are tested.

- **Ext×5 — larger external larval supply.** The external outbreak production (pre-pelagic recruits/ha/yr, migrated) is multiplied by 5 in all years. Its timing remains the evidence-derived calendar-year outbreak-probability boundary; only the magnitude, which is unidentified (it was matched, not measured), changes. ×5 is the factor implied by the settlement deficit.
  - *Risk:* trough supply rises proportionally, by about 13–29% of trough settlement × 5.
  - No calendar-specific pulse is used; that would be an oracle.
- **Juv×5 — undetected 1991 juvenile cohort.** The 1991 juvenile and subadult preload is ×5 (30:15:1 per adult instead of 6:3:1); adults are unchanged and consistent with `q`. Small and cryptic COTS are largely undetected by manta tows, so the 1991 juvenile stock is not observed.
  - *Risk:* a one-off early pulse that cannot sustain a 1994–98 rise.

### Factorial

| Factor | Levels |
| --- | --- |
| Initiation | none, Ext×5, Juv×5, both |
| Mechanisms | none, A (gradual Allee), P+A (predator pit `Q=20` + gradual Allee) |

That gives 12 arms, all on v0.3 with the working model (per-capita, H1, H2, condition), `q=0.015`, `k_unit=19.61`, Allee scale 3, and the pinned RNG seed. Reference: `none_none` must replay `bgeo_none` exactly.

Report:
- first-peak year and height (Lizard, MacGillivray observed 1996–98);
- 1994–98 rise rate;
- trough adults/ha (observed about 0–3);
- second-wave height (observed 0.67–0.85, 2013–15);
- external vs internal settlement;
- site-area coral against LTMP.

### Gate, failures and rollback

The gate is unchanged. Failures:
- inexact reference replay;
- external supply differing from the stated multiplier;
- non-finite states;
- loss of all first waves;
- an improvement achieved only by making the boundary dominate all settlement in troughs.

Rollback: multiplier 1 and preload 6:3:1 (defaults). Neither the multiplier nor `Q` may be tuned to observed peak years.

### Initiation decision (2026-10-12)

[v3 run](runs/20261012T_10_initiation_1991_v3/metadata.toml) and [audit](runs/20261012T_10_initiation_1991_v3/initiation_audit.json).

**Disclosed operational failures, both before any simulation:**
- v1: the new cyclone loader indexed a DataFrame with a DimensionalData year vector. Fixed; a test was added.
- v2: the runner's manta-rescaling guard rejected v0.3's zero-cover sites. The guard now skips zero-cover sites; v0.2 behaviour is unchanged.

**Compatibility.** `none_none` replays the run-time `bgeo_none` to a relative delta of `6e-15`, so the v0.3 rebuild is numerically equivalent. The audit uses a relative tolerance of `1e-9` because cyclone mortality is now computed with Python's `exp`. External settlement equals exactly 1× or 5× the reference. **No arm passes.**

| Arm | Mean first peak | Matched | Two-wave reefs | Mean 2010–17 max CPUE | Trough minimum (adults/ha) |
| --- | --- | --- | --- | --- | --- |
| none_none | 2001 | 2/7 | 0 | 0.30 | 7.5 |
| ext5_none | 1998 | 2/7 | 0 | 0.26 | 7.4 |
| ext5_PA | 1999.7 | 4/7 | 1 (+ Eyrie) | 0.54 | 8.7 |
| both_PA | 2004 | 5/7 | 1 (+ Eyrie) | 0.43 | 7.0 |
| juv5_none | 1998.7 | 3/7 | 1 | 0.21 | 5.7 |

Against the predictions:
- Ext×5 fixed the initiation deficit as hypothesised: early settlement rises from ~17 to ~76 COTS/ha/yr and the Lizard 1994 density matches observations. As risked, trough supply also rose, so the boundary provides 32–58% of gap settlement.
- Juv×5 gave only a one-off early pulse.
- With Ext×5, the predator pit no longer delays initiation, because wave-phase supply (60–190) swamps `Q=20`. It now helps produce correctly timed 2014 second waves on Lizard and Eyrie.

**Decision.** Nothing is promoted. Ext×5 + P + A on v0.3 is the leading candidate configuration. The remaining failures are the trough floor (6–10 vs about 0–3 adults/ha), the abrupt observed 1999 crash, and second waves on MacGillivray and North Direction.

**Next gate.** Preregister a predator-pit capacity sensitivity under Ext×5 (e.g. `Q` = 20, 50, 100 COTS/ha/yr), bounded by field juvenile-predation evidence. It must not be fitted to troughs. Its expected cost is slower or smaller second waves if `Q` exceeds wave-phase juvenile supply.


## Project-lead decisions and preregistered predator-pit capacity tuning (2026-10-12)

Status: **preregistered before implementation; no default or objective change.**

### Decisions by the project lead

- **External ×5 is the working initiation assumption.** It is labelled as an assumption: a magnitude on the never-measured external outbreak production term. Timing stays on the evidence-derived calendar boundary. The A and P+A arms with and without the juvenile preload behaved much like Ext×5, so Juv×5 is dropped.
- **`Q` is tuned, not evidence-bounded.** No field study gives the juvenile predation capacity on the model's effective juvenile scale. It is therefore a calibrated parameter and must be labelled as such in all outputs.

### What the field evidence does and does not give

Wolfe et al. 2025 (CCIP-P-05 technical report, *Resolving the impact of benthic and cryptic predation on the crown-of-thorns starfish*) and Desbiens et al. 2023 (*Coral Reefs* 42:579):

- **Consumption rates.** In rubble mesocosms with alternative prey, *Schizophrys aspera* ate 5.8 and portunid crabs 0.4 newly settled COTS (about 2 mm) per day. Predation largely stops at a size escape of about 10 mm (*S. aspera*) or 5 mm (portunids), roughly the first 150 days.
- **Predator densities near Lizard.** In the northern GBR, *S. aspera* was about 0.08 per 100 m² and portunids about 1.6 per 100 m² of rubble. Heron Reef had about 0.77 *S. aspera* per 100 m².
- **Implied capacity.** At Lizard this gives about 8 × 5.8 × 150 ≈ 7,000 *S. aspera* consumption plus about 160 × 0.4 × 150 ≈ 9,600 portunid consumption, so roughly 10⁴ settlers per ha of rubble per year. This is a mesocosm potential. The eDNA work found realised rates 1.6–3× lower in outbreak hotspots.
- **Mismatch with `Q`.** This predation acts on settlers in their first months, whereas `Q` acts on effective age-1 juveniles after background mortality. The field numbers therefore do not bound `Q` from above. They show only that an annual capacity well above the model's trough juvenile densities (about 10–35/ha) is ecologically plausible, and that Lizard sits in a low-predator region.
- Keesing, Halford & Hall 2018 (*MEPS* 597:179) remains the evidence that juvenile mortality is mostly predation.

### Hypothesis

Under Ext×5, wave-phase juvenile supply (about 60–250/ha/yr) swamps `Q = 20`, while gap supply (30–60/ha/yr) is only half removed. A larger `Q` in the range of gap supply should remove most gap juveniles and deepen troughs towards 0–3 adults/ha, while outbreak cohorts still escape.

**Expected cost:** as `Q` approaches wave-phase supply, first and second waves shrink, are delayed, or are lost. The tuned value should sit where troughs deepen before the waves are lost. If no such value exists, the predator pit cannot explain the troughs under this supply.

### Grid

All arms use v0.3 with the working model (per-capita, H1, H2, condition), External ×5, gradual Allee (A), `q = 0.015`, `k_unit = 19.61`, Allee scale 3, and the pinned seeds.

| Arm | P | `Q` (COTS/ha/yr) | Role |
| --- | --- | --- | --- |
| `ext5_A` | off | – | Reference; must replay `ext5_A` of the [initiation run](runs/20261012T_10_initiation_1991_v3/metadata.toml) exactly |
| `ext5_PA_q20` | on | 20 | Must replay `ext5_PA` of the initiation run exactly |
| `ext5_PA_q30` … `ext5_PA_q500` | on | 30, 45, 70, 100, 150, 220, 330, 500 | Log-spaced, about ×1.5 per step |

The grid is fixed in advance and is not refined after viewing results. A bounded one-dimensional grid is a screen, not a production optimisation; the AGENTS.md rule against optimisers before a gate pass still applies.

### Frozen selection rule

- **Objective.** The frozen `fixed_q0015` `smoothed_3y` cycle loss (`cots_cycle_score`, unchanged), summed over Lizard, MacGillivray, North Direction and Eyrie.
- **Selection.** `Q*` is the grid value with the lowest summed loss among arms that keep a detected first peak on every reef. Ties go to the smaller `Q`.
- **Held-out reefs.** For each reef, `Q` is reselected on the other three using the same rule. Report the held-out reef's loss and peaks at that `Q` against the `ext5_A` reference.
- **Sensitivity (gate 5).** Report the range of `Q` whose summed loss is within 10% of the minimum. A single-point optimum is flagged as narrow.

Secondary reports, not used for selection:
- the qualitative gate;
- trough minimum adults/ha (observed about 0–3);
- gap juvenile removal and the share of removed juveniles;
- second-wave year and height;
- first-peak year and height;
- gap external share;
- site-area coral against LTMP.

### Failures and rollback

Failures:
- an inexact replay of either reference arm;
- external supply other than exactly ×5 the base;
- non-finite or negative states;
- an improvement that comes only from extinguishing waves on some reefs (caught by the first-peak constraint).

Even the selected `Q*` is a single-seed, in-sample candidate. It requires replicated seeds and held-out years before any promotion.

Rollback: `COTS_JUVENILE_PREDATION=false` (the default).

### Predator-pit tuning decision (2026-10-12)

[Run](runs/20261012T_11_predation_capacity_1991_v1/metadata.toml) and [audit](runs/20261012T_11_predation_capacity_1991_v1/predation_capacity_audit.json). Both reference arms replay the initiation run exactly (max delta `0.0`), and external supply is exactly Ext×5 in every arm. The audit's first attempt failed while serialising numpy integers to JSON; its partial outputs were removed and the audit rerun. **No arm passes the gate.**

| `Q` | Summed loss | Matched | Reefs with 2 peaks | Mean first peak | Mean trough min (adults/ha) | Mean 2010–17 max CPUE |
| --- | --- | --- | --- | --- | --- | --- |
| off | 19.90 | 2/7 | 0 | 1998.0 | 7.7 | 0.28 |
| 20 | 15.28 | 4/7 | 2 | 2000.0 | 8.6 | 0.63 |
| **30** | **12.96** | 5/7 | 3 | 2000.5 | 7.0 | 0.94 |
| 45 | 14.55 | 5/7 | 3 | 2001.3 | 7.9 | 0.78 |
| 70 | 14.03 | 5/7 | 3 | 2001.5 | 8.0 | 0.65 |
| 100 | 14.63 | 5/7 | 3 | 2002.3 | 7.3 | 0.58 |
| 150 | 18.16 | 4/7 | 2 | 2003.0 | 7.1 | 0.35 |
| 220 | 18.06 | 3/7 | 3 | 2003.5 | 4.4 | 0.34 |
| 330 | 15.85 | 4/7 | 4 | 2003.5 | 3.3 | 0.26 |
| 500 | 16.63 | 4/7 | 4 | 2003.5 | 2.2 | 0.28 |

**Selection.** `Q* = 30` COTS/ha/yr; the plateau within 10% is `Q` 30–70. Leave-one-reef-out reselects `Q = 30`, or `Q = 20` with North Direction held out. Held-out loss beats the P-off reference on Lizard (1.69 vs 4.03) and Eyrie (2.73 vs 5.47), but not on MacGillivray (5.43 vs 4.83) or North Direction (5.81 vs 5.56). Per-reef loss at North Direction jumps between grid points (5.8, 3.1, 3.7, 2.7, 2.7, 5.8), because peak detection is near a threshold there.

Against the predictions:
- **Wave survival.** Outbreak cohorts escape at moderate `Q`, as hypothesised.
- **Gap depth.** Troughs did **not** deepen at moderate `Q`. Gap settlement *rises* with `Q` (46 → 70/ha/yr at `Q` 0 → 150), because the pit leaves more coral, which raises condition and fecundity.
- **Observed-depth troughs.** These appear only at `Q ≥ 330`, where the pit also suppresses the first wave (heights 0.2–0.9 COTS/tow in 2002–04). This is the predicted cost.
- **Initiation.** P delays the first peak by 2–6 years even under Ext×5, because early supply and gap supply have similar magnitudes.

**Decision.** Nothing is promoted. `Q = 30` is the **tuned working value** for the leading candidate (Ext×5 + P(`Q=30`) + A on v0.3), labelled as calibrated. A single per-hectare predator sink cannot reproduce both a fast 1994–98 rise and 0–3/ha troughs.

**Open structural gaps:**
- the trough floor of 5–10 adults/ha;
- the first peak about 3 years late with P on;
- the abrupt 1999 crash (no mechanism removes standing adults quickly);
- MacGillivray's second wave;
- North Direction's second wave, which is too high.

**Candidate next screens** (for the project lead to choose; none preregistered):
1. **Size-dependent manta detectability (observation model).** Gap adults are mostly young (age 2–3) and small. Manta detectability falls for small, cryptic COTS, so a constant `q` may overstate trough CPUE. This changes the frozen observation model and needs explicit approval and a new objective label.
2. **Post-outbreak adult loss.** For example, a disease or starvation pulse (a state-triggered hazard, as ReefMod uses) to reproduce the abrupt 1999 crash, which also lowers early-gap fecundity.
3. **Phase-dependent predator capacity.** For example, predators that track rubble, or a numerical response. This separates initiation from gap predation, but needs evidence to avoid becoming a free parameter.


## Preregistered post-outbreak collapse × fast-coral growth screen (2026-10-12)

Status: **preregistered before implementation; no default, domain-package or objective change.** The project lead asked for (1) the ReefMod-style post-outbreak adult die-off, aimed at the abrupt observed 1999 crash, and (2) a test of whether limited fast-growing coral prevents second peaks reaching >0.5 COTS/tow, including the supervisor's suggestion to adjust the CoralBlox growth scalar.

### Pre-implementation findings

**Coral growth factors in use.** All Lizard runs use **calibration group 1** for every site. The Lizard loader sets `CB_CALIB_GROUPS = 1` when the column is missing. In the GBR-wide RME layer (`sandbox/data/rme_ml_2025_06_05/data_files/region/reefmod_gbr.gpkg`), the 113 domain reefs belong to group 5 (58 reefs), group 1 (44) and group 4 (11). **All four calibration reefs are in group 5.**

| Factor (fast prey are groups 1–3) | Group 1 (in use) | Group 5 (RME assignment for Lizard, MacGillivray, North Direction, Eyrie) |
| --- | --- | --- |
| Linear-extension scale, tabular *Acropora* | 0.883 | 1.006 |
| Linear-extension scale, corymbose *Acropora* | 0.748 | 0.863 |
| Linear-extension scale, corymbose non-*Acropora* | 1.388 | 0.765 |
| `mb_rate` scale, tabular / corymbose *Acropora* / corymbose non-*Acropora* | −0.998 / −0.058 / −0.265 | −0.256 / −0.974 / −0.280 |
| Growth acceleration height / midpoint / steepness | 1.595 / 0.250 / −19.7 | 1.753 / 0.287 / −16.7 |

- **How the factors are applied.** The effective extension is base `linear_extensions(:calib)` × linear-extension scale × a site growth multiplier. That multiplier is the growth acceleration below 5% relative cover, a blend between 5% and 80%, and the space-limited factor above 80%. Survival is `s + 0.5(1 − s)·mb_rate_scale`.
- **Base extension of fast groups (cm/yr, size classes 1–6):**
  - tabular *Acropora*: 2.7–7.3;
  - corymbose *Acropora*: 2.7–3.9;
  - corymbose non-*Acropora*: 1.6–1.9.
- **Default calibration group in the CoralBlox factors.** The CB-group factors are defaults from ADRIA's coral factor tables. The COTS calibration never overrode them.

**Fast-coral diagnostic** (leading candidate, `ext5_PA_q30`):
- **Peak fast cover 2008–13:** 0.10–0.14 of habitable area on Lizard, Eyrie and North Direction, where body condition reaches 0.88–0.99 and second waves occur. MacGillivray peaks at only 0.038 (condition ≤ 0.51) and has no second wave.
- **Fast share of coral:** about 0.36 in 1991 on every reef. With `Q = 500`, which suppresses COTS, the share stays about 0.37 on MacGillivray. With `Q = 30` it is 0–0.09 from 2000 to 2013.
- **Total-coral recovery:** model site-area recovery already matches or exceeds LTMP on these reefs.

So MacGillivray's missing fast coral is a consequence of COTS grazing: the first wave removes it and trough adults (7–17/ha), feeding on preferred prey first, keep it down. Slow growth is not the cause. Faster fast-coral growth and a lower trough floor are therefore competing explanations, and the factorial must separate them.

### Collapse mechanism: ReefMod 7.0 precedent

[`settings_COTS.m`](https://github.com/ymbozec/REEFMOD.7.0_GBR/blob/main/settings/settings_COTS.m) and [`f_runmodel.m`](https://github.com/ymbozec/REEFMOD.7.0_GBR/blob/main/functions/f_runmodel.m):
- **Trigger.** Adult density must be at or above `COTS_density_threshold_for_disease = 2.7` per 400 m² ("~1 per tow") at every half-yearly step of a look-back of `2 × duration` steps. The duration is drawn from `[4, 6]` **years** when the threshold is hit.
- **Reset.** Classes 4–16 (ages ≥ 18 months) are rescaled, keeping their age proportions, so they sum to `COTS_background_density = 0.01` per 400 m². The larvae produced that step are kept.
- **Lizard option.** The earlier audit found a commented fixed two-year duration for Lizard calibration.

These are modelling conventions, not field estimates of pathogen timing. The disease link is plausible (Rivera-Posada et al. 2012) but unproven, so the mechanism is **pathogen-compatible phenomenology**.

Unit transfer is purely by area, because state is COTS/ha: 2.7 per 400 m² = **67.5 COTS/ha**, and 0.01 per 400 m² = **0.25 COTS/ha**. Neither depends on `q`. Under `q = 0.015`, 67.5/ha corresponds to about 1.0 COTS/tow, which matches ReefMod's comment.

### Mechanism (opt-in COTSMod switch `outbreak_collapse`)

After each annual transition, in COTS/ha:
1. The clock counts consecutive years in which post-transition adults `N3` are at or above `collapse_threshold_ha` (67.5). A year below the threshold resets it to 0.
2. When the clock reaches `collapse_duration_years` (`D`), adults are rescaled, keeping age-class proportions, to `min(N3, collapse_background_ha)` (0.25), and the clock resets.
3. With `collapse_juveniles = true`, `N2 + N3` together are rescaled to the background density. This is the ReefMod reset of classes ≥ 18 months, extended to all of N2, and it tests whether immature cohorts refill the trough.
4. N1, external supply and the larvae produced that year are unchanged. Grazing in that year uses the post-collapse adults, as in all post-transition grazing.

The clock is state-triggered, not calendar-triggered. It must not be tuned to observed peak years.

### Coral arms (run time; no domain rebuild)

| Arm | Calibration groups | Fast linear extension |
| --- | --- | --- |
| `g1` | group 1 everywhere (current) | default |
| `grme` | per reef from the RME layer, relabelled to contiguous IDs (RME 1, 4, 5 → 1, 2, 3) carrying each RME group's default factor values | default |
| `grme_fast120` | as `grme` | linear-extension scale ×1.2 for tabular *Acropora*, corymbose *Acropora* and corymbose non-*Acropora* |

- **Why relabel.** ADRIA's sampler keeps only the CB-group factors present in the domain, and `scenario.jl` indexes them by group number, so non-contiguous groups would be misindexed. This is an ADRIA limitation to report, not to fix here.
- **`grme` status.** It is a data correction candidate. Adoption requires a v0.4 domain rebuild by decision of the project lead.
- **`grme_fast120` status.** It is a labelled sensitivity, not a calibration: no fast-coral composition data are available locally to fit recovery rates. LTMP manta tows give hard coral only, and the photo-transect file has no data for these reefs.

### Factorial

Coral {`g1`, `grme`, `grme_fast120`} × collapse {`none`, `D2` (adults), `D2J` (adults + N2), `D4` (adults)} gives **12 arms**:
- the leading candidate on v0.3: working model + Ext×5 + P (`Q = 30`) + gradual Allee;
- `q = 0.015`, `k_unit = 19.61`, Allee scale 3, pinned seeds;
- threshold 67.5/ha and background 0.25/ha are fixed.

`D2` follows ReefMod's commented Lizard option; `D4` is ReefMod's lower bound.

**Reference.** `g1_none` must replay `ext5_PA_q30` of the [Q tuning run](runs/20261012T_11_predation_capacity_1991_v1/metadata.toml) exactly.

### Expected consequences

- **D2:**
  - an abrupt collapse 1–2 years after first-wave adults exceed 67.5/ha, close to the observed 1.92 → 0.05 COTS/tow drop;
  - lower early-gap adults and fecundity, faster fast-coral recovery, and possibly earlier second waves;
  - *risk:* N2 refills the trough within 1–2 years.
- **D2J:** removes that refill. *Risk:* second waves then rely on external and internal settlement alone.
- **D4:** rarely triggers, because modelled peaks exceed 67.5/ha for only about 2–4 years. It is a ReefMod-faithful negative control.
- **`grme` / `grme_fast120`:** if fast-coral growth limits second waves, MacGillivray's fast cover and second-wave height rise without the collapse. If grazing is the cause, they help mainly when combined with D2/D2J.
- **Coral against LTMP:** must not degrade materially, especially with `fast120`.

### Report

- the qualitative gate and the frozen summed loss (`fixed_q0015`, `smoothed_3y`);
- collapse years by reef and adults removed;
- 1998–2001 drop rate (observed Lizard 1.92 → 0.05 over one year);
- trough minimum adults/ha (observed about 0–3);
- second-wave year and height (observed 0.67–0.85, 2013–15);
- peak fast cover 2008–13 and body condition;
- site-area coral against LTMP.

### Failures and rollback

Failures:
- an inexact reference replay, or external supply other than Ext×5;
- non-finite or negative states, or adult age classes that do not sum to `N3`;
- a collapse that never triggers in D2 (wiring failure);
- losing all first waves;
- a trough that refills immediately from N2;
- an improvement achieved only by extinguishing waves.

Rollback: `COTS_OUTBREAK_COLLAPSE=false` (the default) and calibration group 1 (the current loader behaviour).

### Collapse × fast-coral decision (2026-10-12)

[Run](runs/20261012T_12_outbreak_collapse_1991_v1/metadata.toml) and [audit](runs/20261012T_12_outbreak_collapse_1991_v1/outbreak_collapse_audit.json).
- **Checks.** `g1_none` replays `ext5_PA_q30` exactly (max delta `0.0`). External supply equals Ext×5 in every arm. Collapse removals are zero with the switch off, and D2 triggers on all four reefs.
- **Tests.** COTSMod 326/326 tests pass, including a new collapse testset: opt-in default, clock, reset arithmetic, age-class mass balance, N2 option and validation. The ADRIA adapter test passes.
- **Result.** **No arm passes the gate.**

| Arm | Summed loss | Matched | Two-wave reefs (of 3) | Mean first peak (COTS/tow) | Steepest drop ratio | Trough min (adults/ha) | Mean 2010–17 max | Peak fast coral 2008–13 | Coral RMSE vs LTMP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| g1_none (ref) | 12.96 | 5/7 | 2 | 1.45 | 0.32 | 7.0 | 0.94 | 0.101 | 0.103 |
| g1_D2 | 14.65 | 5/7 | 2 | 0.94 | 0.13 | 7.7 | 0.48 | 0.044 | 0.086 |
| g1_D2J | 16.00 | 4/7 | 1 | 0.84 | 0.11 | 6.3 | 0.43 | 0.062 | 0.096 |
| g1_D4 | 15.51 | 4/7 | 1 | 1.32 | 0.21 | 8.6 | 0.49 | 0.065 | 0.092 |
| grme_none | 14.61 | 5/7 | 2 | 1.97 | 0.32 | 7.2 | 0.91 | 0.103 | 0.133 |
| grme_D2 / D2J / D4 | 18.71 / 16.32 / 18.93 | 3–4/7 | 0–1 | 0.86–1.62 | 0.03–0.08 | 4.5–7.6 | 0.31–0.48 | 0.04–0.06 | 0.13–0.17 |
| grmefast_none | 18.82 | 6/7 | **3** | 2.56 | 0.20 | 11.5 | 1.73 | 0.168 | 0.137 |
| grmefast_D2 / D2J / D4 | 19.61 / 15.50 / 19.81 | 3–4/7 | 0–1 | 0.87–1.80 | 0.01–0.04 | 5.0–8.1 | 0.35–0.47 | 0.04–0.07 | 0.13–0.18 |

**Against the predictions — collapse.**
- **Abrupt drop.** D2 produced the predicted abrupt first-wave drop (Lizard 2000; North Direction collapses in 1999, 2001 and 2003).
- **Rebound and lost second waves.** The predicted N2 refill risk occurred, and D2J showed it is not only N2. Adults are reset before coral is exhausted, so the next cohorts of settlers find food. The first wave becomes a 2–3 pulse sequence, and gap fast coral stays near the 0.05 starvation threshold instead of recovering to 0.10–0.14. Second waves therefore shrink or vanish.
- **D4.** It triggered on 2–4 reefs; peaks stay above 67.5/ha for longer than expected.
- **Mechanism.** In this model the second wave depends on a starvation-driven coral reset: COTS exhaust fast coral, starve, and the bare reef then regrows. A density-triggered disease reset decouples the crash from coral depletion and breaks that cycle.
- **Data.** Observed Lizard coral fell with the crash (LTMP 0.19 → 0.08 by 1999), consistent with a starvation-coupled crash.

**Against the predictions — coral.**
- **`grme`.** As expected, it changes first-wave height more than second-wave occurrence: higher peaks and higher coral, but worse coral RMSE.
- **`grmefast`.** If fast-coral growth limits second waves, MacGillivray should gain one; it did, for the first time in any screen. Second waves on all three gate reefs show that fast-coral availability is a binding constraint on second-wave occurrence. At ×1.2, though, waves overshoot (MacGillivray 2017 at 1.74, North Direction 2016 at 2.50 vs observed 0.67–0.71) and arrive 2–3 years late. Faster coral regrowth also feeds larger first peaks (up to 3.3) and higher troughs.

**Decision.** Nothing is promoted.
- **Collapse switch.** It stays opt-in and is **not recommended**: it is pathogen-compatible phenomenology that is counterproductive in this model.
- **Calibration-group mismatch.** The Lizard loader's group 1 versus RME group 5 is recorded as an input issue. Correcting it would need a v0.4 domain rebuild by decision of the project lead. Group 5 alone worsens the coral fit.
- **Fast-coral growth.** It is the strongest lever found for second-wave occurrence. It should be constrained by observed coral recovery, not by COTS peaks.

**Known ADRIA limitation.** Non-contiguous `CB_CALIB_GROUPS` are misindexed, because the sampler drops absent groups and `scenario.jl` indexes by group number. Runs relabel groups contiguously.

**Next gate (candidates for the project lead).**
1. **Fast-coral growth constrained by data.** A bounded fast linear-extension grid (e.g. ×1.0–1.2 on the chosen calibration group) scored against observed coral recovery, ideally with group-level (Acropora / fast) cover from LTMP photo transects or the COTS-Coral-OutcomesV2 data, not against COTS peaks. Then rerun the COTS screen with that fixed value.
2. **A faster starvation crash.** A sharper or shorter-memory preferred-prey starvation response (memory below 0.5 yr), so the crash follows coral exhaustion within a year, as observed in 1999.


## Preregistered fast-coral growth fit × faster starvation screen (2026-10-12)

Status: **preregistered before implementation; no default, domain-package or objective change.**

### Decisions by the project lead (2026-10-12)

1. **Calibration groups.** Use **CoralBlox calibration group 5** for the calibration reefs. It is implemented as the per-reef RME assignment (`grme` from the [collapse screen](runs/20261012T_12_outbreak_collapse_1991_v1/metadata.toml)): all four calibration reefs are group 5, and the other domain reefs keep their RME groups (5: 58 reefs, 1: 44, 4: 11). Because ADRIA misindexes non-contiguous groups, the runner relabels them contiguously at run time. A v0.4 package with a `CB_CALIB_GROUPS` column must wait for an ADRIA indexing fix.
2. **Fast-coral growth.** Fit the fast-coral linear-extension multiplier (×1.0–1.2) to observed fast-coral recovery, then rerun the COTS screen.
3. **Starvation.** Test a faster starvation-coupled crash.

### Observations for the coral fit

**Source.** LTMP photo transects (`sandbox/data/reef_photo_transect.csv`): `purpose = COMPOSITION`, `variable = HARD CORAL`, 9 m, `median`, a fraction of benthos.

**Model fast groups.** Tabular *Acropora*, corymbose *Acropora*, and corymbose non-*Acropora* ("and Pocillopora"). The observed equivalent is `reefpage_category` *Acropora* + Pocilloporidae. *Isopora* and *Montipora* are not reported separately on these reefs; they fall in "Other" and are excluded.

**Reef matching:**
- LTMP Lizard Isles → Lizard Island Reef;
- Macgillivray Reef → MacGillivray Reef;
- North Direction Island → North Direction Reef. This is a caveat: the LTMP site is on the island's fringing reef, not the separate North Direction Reef.

Eyrie has no photo transects.

| Observed fast cover | Low (year) | Recovery high (year) | 2017 → 2024 |
| --- | --- | --- | --- |
| Lizard Isles | 0.016 (2001) | 0.078 (2013) | 0.002 → 0.062 (peak 0.093 in 2023) |
| Macgillivray Reef | 0.009 (1998) | 0.042 (2011) | 0.004 → 0.064 |
| North Direction Island | 0.020 (1997) | 0.113 (2013) | 0.005 → 0.110 |

**Frozen coral criterion.**
- **Primary.** The least-squares slope of fast cover against year over **2017–2024**, observed (survey years) against modelled (annual values, model `fast_coral` as a fraction of habitable area). The window follows the 2016–17 bleaching and cyclone reset, when COTS were rare in both data and model, so it isolates coral recovery from COTS grazing as far as the data allow. Photo transects sit on hard substrate on the 9 m slope, so they are compared with fast cover per habitable area. The multiplier `m*` is the grid value minimising the summed squared slope error over the three reefs, using the reference starvation arm; ties go to the smaller `m`.
- **Secondary, not used for selection:**
  - the 2001–2013 recovery slope (COTS-confounded);
  - full-series RMSE of fast cover against photo transects;
  - total site-area coral RMSE against LTMP manta.

### Starvation diagnosis (leading candidate, Lizard)

The first-wave decline takes 2–3 years (adults 106 → 65 → 25 → 10/ha over 2002–05) for two reasons:
- **Lag.** Preferred-food memory is updated with fast cover *before* that year's grazing. In 2000 memory was 0.090 while end-of-year fast cover was 0.023, so adults starve the year *after* the coral is exhausted.
- **Refill.** A large juvenile pipeline (subadults 150/ha, recruits 240/ha in 2001) keeps maturing about 50 adults/ha/yr.

The memory already weights the current year by `1 - exp(-1/0.5) = 0.865`, so shortening it alone should change little.

### Starvation arms

- **`ref`:** H1 on remembered pre-grazing fast cover, memory 0.5 yr (current).
- **`M0`:** memory 0, so the current year has weight 1. This is the user's suggestion; the expected effect is small.
- **`PG` — post-grazing preferred food** (new opt-in COTSMod switch `starvation_post_grazing`). H1 memory is updated with the preferred cover left after this year's per-capita demand of the current adults:

  `F_food = max(F - N3·c / (10⁴·k), 0)`, in fractions of habitable area.

  Here `c` = 10 m²/adult/yr and `k` is the site habitable fraction. This is parameter-free and requires per-capita consumption and H1. Body condition still uses cover `F`. Starvation is then coupled to coral exhaustion within the same year, and its strength scales with COTS density: 10 adults/ha remove 0.017 of habitable cover per year at `k = 0.6`, and 100/ha remove 0.17.

### Factorial

Fast multiplier `m ∈ {1.0, 1.05, 1.10, 1.15, 1.20}` × starvation {`ref`, `M0`, `PG`} gives **15 arms**. All use the leading candidate (v0.3, working model, Ext×5, `Q = 30`, gradual Allee), per-reef RME calibration groups, `q = 0.015`, `k_unit = 19.61`, Allee scale 3 and pinned seeds.

**Reference.** `m100_ref` must replay `grme_none` of the collapse run exactly.

**Two-stage reading:**
1. Select `m*` from the `ref` arms with the coral criterion; COTS metrics are not used.
2. Report the COTS gate, frozen summed loss, trough minima, crash speed and second waves for `ref`, `M0` and `PG` at `m*`. The full grid is shown as sensitivity.

### Expected consequences

- **`m`:** a higher `m` raises fast-coral recovery slopes. At ×1.2 the collapse screen gave second waves on all reefs but overshot. `m*` is expected near the lower end if the model's 2017–24 recovery already matches.
- **`M0`:** about one year faster at most; small changes.
- **`PG`:**
  - the crash follows exhaustion within a year (steepest one-year drop ratio well below 0.3);
  - lower peak heights;
  - deeper troughs, because trough adults face post-grazing food below θ while fast cover is 0.03–0.06;
  - faster fast-coral recovery in the gap.
  - *Risk:* suppressed second-wave build-up if `F_food` stays below θ as COTS rebuild.

### Failures and rollback

Failures:
- an inexact reference replay, or external supply other than Ext×5;
- non-finite or negative states;
- `m*` at a grid boundary (reported as unidentified within the bound);
- losing all first waves;
- an improvement only by extinguishing waves.

Rollback: `COTS_STARVATION_POST_GRAZING=false` (default), memory 0.5, multiplier 1.0, and calibration group 1 (current loader default).

### Growth × starvation decision (2026-10-12)

[Run](runs/20261012T_13_growth_starvation_1991_v1/metadata.toml) and [audit](runs/20261012T_13_growth_starvation_1991_v1/growth_starvation_audit.json).
- **Checks.** `m100_ref` replays `grme_none` of the collapse run exactly (max delta `0.0`), and external supply is exactly Ext×5.
- **Tests.** COTSMod 333/333 pass, including a new post-grazing starvation testset. The ADRIA adapter test passes.
- **Result.** **No arm passes the gate.**

**Stage 1 — coral fit (`ref` arms).** The values are model / observed 2017–24 fast-cover slopes per year.

| `m` | Criterion (×10⁻⁴) | Lizard | MacGillivray | North Direction |
| --- | --- | --- | --- | --- |
| 1.00 | 1.39 | 0.0083 / 0.0117 | 0.0042 / 0.0087 | 0.0040 / 0.0143 |
| 1.05 | 1.22 | 0.0118 | 0.0050 | 0.0039 |
| **1.10** | **1.21** | 0.0093 | 0.0056 | 0.0041 |
| 1.15 | 1.80 | 0.0100 | 0.0217 | 0.0174 |
| 1.20 | 2.61 | 0.0102 | 0.0242 | 0.0188 |

- **Selection.** `m* = 1.10` is interior. The criterion is nearly flat from 1.00 to 1.10, so `m` is only weakly identified.
- **Above 1.10.** MacGillivray and North Direction have COTS outbreaks inside the fit window (mean adults 20–33/ha in 2017–24), so the low-COTS assumption fails.
- **North Direction.** Its recovery is underestimated at every `m ≤ 1.10`, which may partly reflect the island-versus-reef site mismatch.

**Stage 2 — COTS at `m* = 1.10`.**

| Arm | Summed loss | Matched | Two-wave reefs | Mean first peak (COTS/tow) | Steepest drop ratio | Trough min (adults/ha) | Mean 2010–17 max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| m100_ref (group 5, `m` 1.0) | 14.61 | 5/7 | 2 | 1.97 | 0.32 | 7.2 | 0.91 |
| **m110_ref** | 15.93 | 5/7 | 2 | 2.33 | 0.20 | 9.1 | 0.99 |
| m110_M0 | 15.65 | 5/7 | 2 | 2.09 | 0.22 | 8.0 | 0.81 |
| m110_PG | 17.73 | 3/7 | 0 | 0.91 | 0.31 | 7.0 | 0.32 |

Against the predictions:
- **`M0`.** Negligible, as predicted, because the 0.5-yr memory already weights the current year 0.865.
- **`PG`.** The density-scaled starvation acts before the peak, so peaks are a third to a half of the reference and broad. It removes second waves, the stated risk, and does not steepen the crash.
- **Crash speed.** The juvenile pipeline (subadults 150/ha, recruits 240/ha at Lizard's peak) keeps refilling adults, and neither arm addresses it.

**Diagnostic finding (not preregistered, from the photo-transect comparison).**

| 1993–95 | Observed fast cover / share of hard coral | Model fast cover / share |
| --- | --- | --- |
| Lizard | 0.023–0.039 / 0.18–0.24 | 0.18–0.21 / 0.48–0.52 |
| MacGillivray | 0.022–0.036 / 0.10–0.14 | 0.15–0.17 / 0.36–0.39 |
| North Direction | 0.036–0.064 / 0.19–0.31 | 0.17–0.18 / 0.47–0.51 |

- **What the transects show.** Observed fast coral stayed at 0.01–0.03 through the first outbreak and rose steadily through the 2000s.
- **What the model does.** It starts with 3–5× the fast coral, which fuels oversized first waves (2–3 COTS/tow). It then grazes fast coral to zero and rebounds to 0.11–0.16 by 2011.
- **Source of the excess.** The fast-coral excess is in the initial state: inherited 2023 RME composition, rescaled to 1991 manta total cover.

Caveat: photo transects are fixed 9 m slope sites, whereas the model value is an area-weighted reef aggregate.

**Decision.** Nothing is promoted.
- **Working values.** CoralBlox group 5 (project lead) and `m = 1.10` (fitted, weakly identified) become working values. The leading candidate is v0.3 + group 5 + `m` 1.10 + working model + Ext×5 + `Q = 30` + gradual Allee.
- **Re-tuning.** `Q = 30` was tuned under group 1 and should be re-checked after the next change.
- **Switches.** `COTS_STARVATION_POST_GRAZING` stays opt-in and is not recommended. Memory stays 0.5 yr.

**Next gate (proposed).** A data-derived 1991 coral composition: set the initial fast/slow split per calibration reef from LTMP photo-transect composition (*Acropora* + Pocilloporidae share of hard coral, 1993–95 mean), keeping total cover from the 1991 manta surface. Other reefs use the domain-wide mean observed share. Rescale within groups so size structure is preserved. Compare with the inherited composition on the leading candidate; this is an input correction, not a tuned parameter. Then re-check `Q`.


## Preregistered data-derived 1991 coral composition screen (2026-10-12)

Status: **approved by the project lead (2026-10-12); preregistered before implementation.** No default, domain-package or objective change.

### Correction to the previous decision text

The recent arms (`coral = "inherited_2023"`) start from the domain's initial cover **as-is**. That cover is the ReefMod **2023** mean cover, divided by the ReefMod habitable fraction and split into 35 group-size bins (`sandbox/domain_building/README.md`). It is *not* rescaled to 1991 manta totals; that rescaling applies only to the `manta_1989_1991` arms. The 2023 composition is a post-recovery, *Acropora*-rich state, with a fast share of 0.36–0.52 on the calibration reefs. Inherited totals are roughly consistent with 1991 manta (site-area fraction: Lizard 0.16 vs 0.17; MacGillivray 0.30 vs 0.31).

### Hypothesis

The oversized first waves (2–3 COTS/tow) and the fast-coral crash-to-zero then rebound come from an initial fast-coral pool 3–5× the observed. Setting the 1991 fast share from observations should:
- lower first-wave heights towards 1.4–1.9 COTS/tow;
- keep fast coral above zero through the first outbreak, as observed;
- slow its 2000s rebound, as observed.

### Correction (input, not a tuned parameter)

**Target fast share.** For each calibration reef, the target is the 1993–95 mean of the LTMP photo-transect fast share: *Acropora* + Pocilloporidae divided by the sum of all `HARD CORAL` composition categories (`purpose = COMPOSITION`, 9 m).

| Reef | Years | Target fast share |
| --- | --- | --- |
| Lizard Island Reef (from Lizard Isles) | 1994–95 | about 0.21 |
| MacGillivray Reef | 1993–95 | about 0.13 |
| North Direction Reef (from North Direction Island) | 1994–95 | about 0.25 |
| Eyrie Reef and every other domain reef | — | mean of the three, about 0.20 |

The runner computes these from the CSV and records the values and the file hash.

**Applying it at each site.** Total cover and the relative size/group structure within the fast block (groups 1–3) and the slow block (groups 4–5) are kept. Only the fast:slow split is set to the target share. Sites with zero total cover stay zero. A site whose inherited fast or slow block is zero cannot be rescaled; it is left unchanged and counted.

**Caveats:**
- Photo transects are fixed 9 m slope sites, while the model's reef is an area-weighted aggregate.
- North Direction Island is a different reef from North Direction Reef.
- **Starvation risk.** The observed fast cover (about 0.03–0.06 of habitable area) lies near the H1 starvation threshold of 0.05 (ReefMod's 5% preferred-cover convention). With the observed composition, the first wave may be starved or suppressed. **θ will not be tuned to compensate.**

### Factorial (8 arms)

All arms use v0.3, the working model, Ext×5, gradual Allee, CoralBlox group 5 (RME per reef), `q = 0.015`, `k_unit = 19.61`, Allee scale 3 and pinned seeds.

| Arm | Composition | Fast-growth multiplier `m` | `Q` |
| --- | --- | --- | --- |
| `inh_m110_q30` | inherited 2023 | 1.10 | 30 |
| `obs_m110_q20` … `obs_m110_q100` | observed 1993–95 share | 1.10 | 20, 30, 45, 70, 100 |
| `obs_m100_q30`, `obs_m120_q30` | observed | 1.00, 1.20 | 30 |

- **Reference.** `inh_m110_q30` must replay `m110_ref` of [run 13](runs/20261012T_13_growth_starvation_1991_v1/metadata.toml) exactly.
- **Re-checking `Q`.** Under the observed composition, `Q` is re-selected with the frozen rule: summed `smoothed_3y` loss, first peak on every reef, ties to the smaller `Q`.
- **The `m` arms.** These check whether the growth fit is sensitive to composition.

### Report

- the gate and summed loss;
- first-peak year and height;
- fast-cover RMSE against photo transects over 1993–2024;
- whether fast cover stays above zero in 1998–2003;
- trough minimum adults/ha;
- second waves;
- site-area coral RMSE against LTMP manta.

### Failures and rollback

Failures:
- an inexact reference replay, or external supply other than Ext×5;
- total cover not preserved at any site, or the target share not met at rescaled sites;
- non-finite states;
- losing all first waves (this would be reported as evidence of the starvation-threshold tension above, not "fixed" by retuning).

Rollback: inherited composition (`coral = "inherited_2023"`).

### Composition decision (2026-10-12)

[Run](runs/20261012T_14_initial_composition_1991_v1/metadata.toml) and [audit](runs/20261012T_14_initial_composition_1991_v1/initial_composition_audit.json).
- **Run ID.** Auto-sequenced as `20261012T_14_…`; earlier undated runs were renumbered (see `runs/RUN_RENAMES.csv`).
- **Checks.** `inh_m110_q30` replays `m110_ref` exactly (`0.0`), and 1991 total cover is identical across compositions (`0.0`). External supply is Ext×5. No site needed to be left unchanged.
- **Observed fast shares:** Lizard 0.2086, MacGillivray 0.125, North Direction 0.2502, others 0.1946.
- **Result.** **No arm passes the gate.**

| Arm | Summed loss | Matched | Two-wave reefs | Mean first peak (year / COTS/tow) | Trough min (adults/ha) | Mean 2010–17 max | Fast RMSE vs photo | Coral site RMSE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| inh_m110_q30 (ref) | 15.93 | 5/7 | 2 | 2001.5 / 2.33 | 9.1 | 0.99 | 0.063 | 0.124 |
| obs_m110_q20 (selected `Q`) | 14.10 | 4/7 | 1 | 1998.5 / 0.79 | 2.8 | 0.38 | 0.034 | 0.196 |
| obs_m110_q30 | 14.28 | 4/7 | 1 | 1999.0 / 0.80 | 2.8 | 0.37 | 0.035 | 0.195 |
| obs_m110_q45 / q70 / q100 | 14.50 / 15.13 / 15.50 | 4/7 | 1 | 1999.8–2001.0 / 0.78–0.84 | 3.0–3.1 | 0.32–0.35 | 0.037–0.045 | 0.20 |
| **obs_m100_q30** | **12.18** | 5/7 | 2 | 1998.3 / 0.66 | 2.9 | 0.35 | 0.032 | 0.202 |
| obs_m120_q30 | 16.94 | 3/7 | 0 | 1999.8 / 1.11 | 4.2 | 0.45 | 0.037 | 0.187 |

Against the predictions:
- **First waves.** They were lowered as predicted, but too far on Lizard (0.49 vs observed 1.92) and MacGillivray (0.51 vs 0.86). North Direction (1.1) and Eyrie (1.06) remain above observed. Timing moved to 1997–99, matching the observed 1996–98.
- **Fast coral through the first outbreak.** It stayed above zero only partly: the 1998–2003 minimum is 0.003–0.005 against observed 0.01–0.03. It then rebounded slowly, as predicted, and too slowly: Lizard about 0.04 in 2008–13 against observed 0.07–0.08, North Direction 0.035 against 0.10.
- **Troughs.** Unexpectedly, troughs reached the observed range (2–3 adults/ha), because trough fast cover (about 0.03–0.04) stays below the starvation threshold.
- **Starvation risk.** The preregistered risk materialised partly. First waves still occur, but COTS cannot build second waves on 0.03–0.04 fast cover, which is below θ = 0.05 and well below the condition saturation of 0.10.
- **Growth fit.** Under the observed composition, the 2017–24 recovery-slope criterion (×10⁻⁴) is 1.92 / 1.79 / 1.49 for `m` 1.00 / 1.10 / 1.20, against 1.21 for the inherited composition at 1.10. The model slope is about 0.004/yr against observed 0.009–0.014. The run-13 fit of ×1.10 therefore depended on the inherited composition, and extension alone cannot reproduce observed recovery from an observed-size pool.

**Decision.** Nothing is promoted.
- **Composition.** The observed composition is recommended as the data-based initial state for the leading candidate, pending the project lead's decision. It is an input correction with no tuned parameter. `Q` is re-selected as 20, but the loss differences across `Q` are small.

**Two linked gaps for the next gate:**
1. **Fast-coral recovery rate** from an observed-size pool. Candidates are fast-coral recruitment (larval supply and settlement for groups 1–3, including the external/domain coral larval boundary) and the growth acceleration at low cover. These should be fitted to photo-transect recovery (2001–13 and 2017–24) with COTS-off and COTS-on checks.
2. **Food scale.** Observed outbreaks reached 1.9 COTS/tow with transect fast cover at 0.02–0.04.
   - **Data questions first:** does "Other" contain preferred prey (*Montipora*, *Isopora*), and are the 9 m slope transects representative of reef-wide food?
   - **Then the model:** H1 θ (ReefMod's 5% convention) and the condition saturation are defined on model fast cover. Revisiting them needs evidence on the food scale, not a fit to COTS peaks.


## Preregistered food-scale screen: plausible 1991 composition range × starvation grace year (2026-10-12)

Status: **directed by the project lead (2026-10-12); preregistered before implementation.** No default, domain-package or objective change.

### Decisions and reasoning by the project lead

- **Composition.** The observed composition is an appropriate starting point, but the 1991 fast share is uncertain, for three reasons:
  - photo transects start in 1993–94, 2–3 years after the simulation;
  - "Other" may hide preferred *Montipora* and *Isopora*;
  - many *Acropora* plates sit on the upper slope and crest, above the 9 m transects.
  Plausible starting points should therefore be tested as a range.
- **Food scale first.** COTS may need to survive about one year at low food before crashing, so that second peaks can build.

### Evidence

- **Starvation tolerance.**
  - Kettle (1990, JCU thesis) found *A. planci* survives enforced starvation for about **six months**: large adults shrink and die quickly, while small ones keep growing, and gonad indices rose about fourfold after six months.
  - Yang et al. (2022, *Front. Mar. Sci.* 9:1021377) found little change in stomach primary metabolism after four months' starvation.
  - In an annual model, a one-year grace period is therefore the **generous end** of the evidence.
- **Plausible fast share** (1993–95 photo transects, share of all hard coral):

| Reef | *Acropora* + Pocilloporidae | + half of "Other" | + all of "Other" |
| --- | --- | --- | --- |
| Lizard | 0.209 | 0.249 | 0.290 |
| MacGillivray | 0.125 | 0.201 | 0.277 |
| North Direction | 0.250 | 0.383 | 0.516 |
| Mean (Eyrie, other reefs) | 0.195 | 0.278 | 0.361 |

"Other" is a residual category, so the upper bound is generous. It is nonetheless consistent with crest *Acropora* missing from 9 m transects, and it approaches the inherited 2023 share (0.34–0.46).

### Mechanism (opt-in COTSMod switch `starvation_grace_years`, integer years, default 0)

- **Counter.** H1 keeps a per-site count of consecutive annual steps in which remembered preferred food is at or below θ (0.05).
- **Grace.** While that count is ≤ `starvation_grace_years`, the food factor is `f = 1`. After that, the usual cubic H1 factor applies. Food above θ resets the count.
- **Default.** With grace 0 the result is identical to the current model.
- **What still applies.** Body condition, and therefore fecundity, still follows food, so starving COTS breed less while surviving. The factor applies to maturing N2 and adults alike.
- **Simplification.** Grace is a property of the site population, not of individual cohorts.
- **Rules.** It requires H1 and is incompatible with diet-linked starvation (retired).

### Factorial (6 arms)

Composition {`acrpoc`, `halfother`, `other`} × grace {0, 1}, all on the run-14 best configuration:
- v0.3, working model, Ext×5, gradual Allee, `Q = 30`;
- CoralBlox group 5 (RME), fast extension ×1.0;
- `q = 0.015`, `k_unit = 19.61`, Allee scale 3, pinned seeds.

Composition is applied as in run 14 (totals and within-block structure preserved).

**Reference.** `acrpoc_g0` must replay `obs_m100_q30` of [run 14](runs/20261012T_14_initial_composition_1991_v1/metadata.toml) exactly.

### Expected consequences

- **Grace 1:**
  - trough COTS survive one extra low-food year, so troughs may be shallower (risk: back towards 4–10/ha);
  - cohorts arriving as fast coral climbs from 0.03 to 0.05+ survive long enough to build second waves;
  - crashes come about 1 year later;
  - first peaks are higher.
- **Higher fast share:**
  - first peaks rise towards observed heights on Lizard and MacGillivray;
  - risk: the inherited-style overshoot (2+ COTS/tow) and late peaks.
- **Discrimination.** If grace raises second waves only together with deep troughs lost, the food-scale problem lies in the threshold scale, not in starvation timing.

### Report

- the gate and summed loss;
- first-peak year and height;
- trough minimum adults/ha;
- second-wave heights;
- fast RMSE against photo transects;
- coral site RMSE;
- the number of grace years used per reef.

### Failures and rollback

Failures:
- an inexact reference replay, or external supply other than Ext×5;
- totals not preserved;
- non-finite states;
- losing all first waves.

Rollback: `COTS_STARVATION_GRACE_YEARS=0` (default) and the run-14 composition.

### Food-scale decision (2026-10-12)

[Run](runs/20261012T_15_food_grace_1991_v1/metadata.toml) and [audit](runs/20261012T_15_food_grace_1991_v1/food_grace_audit.json).
- **Checks.** `acrpoc_g0` replays run 14 `obs_m100_q30` exactly (`0.0`). 1991 totals are identical across compositions, and external supply is Ext×5.
- **Tests.** COTSMod 350/350 pass (new grace testset); the ADRIA adapter test passes.
- **Operational failure.** The first audit attempt stopped on an unrenamed output key before writing anything.
- **Result.** **No arm passes the gate.**

| Arm | Summed loss | Matched | Two-wave reefs (of 3) | Mean first peak | Trough min (adults/ha) | Mean 2010–17 max | Fast RMSE vs photo | Coral site RMSE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| acrpoc_g0 (ref) | **12.18** | 5/7 | 2 | 1998.3 / 0.66 | 2.9 | 0.35 | 0.032 | 0.202 |
| acrpoc_g1 | 16.67 | 3/7 | 0 | 1999.3 / 0.89 | 3.0 | 0.40 | 0.033 | 0.183 |
| halfother_g0 | 17.75 | 3/7 | 0 | 1999.8 / 1.24 | 4.5 | 0.45 | 0.039 | 0.191 |
| halfother_g1 | 18.92 | 4/7 | 1 | 2000.3 / 1.57 | 4.9 | 1.01 | 0.044 | 0.185 |
| other_g0 | 15.96 | 4/7 | 1 | 2000.3 / 1.52 | 6.5 | 0.60 | 0.048 | 0.174 |
| other_g1 | 19.66 | 4/7 | 1 | 2000.5 / 1.80 | 6.4 | 1.09 | 0.049 | 0.175 |

**Against the predictions:**
- **Grace.** It raised first peaks and let more cohorts build, but second waves appeared only where fast coral recovers above θ (North Direction, overshooting at 2–3 COTS/tow, and Eyrie). Lizard and MacGillivray were unchanged.
- **Higher fast share.** It raised first peaks as expected, too far on North Direction, and raised troughs (the stated risk).

**Diagnosis.** In every arm, Lizard and MacGillivray fast cover in 2000–13 stays ≤ 0.04 (observed Lizard 0.065–0.083 in 2007–13), regardless of starting share or grace. Slow corals grow unchecked:

| Reef, observed `acrpoc` arm | Slow cover (habitable): 1991 → 2003 → 2012 | Total: 1991 → 2012 |
| --- | --- | --- |
| Lizard | 0.21 → 0.50 → 0.62 | 0.26 → 0.66 |
| MacGillivray | 0.30 → 0.60 → 0.70 | 0.34 → 0.74 |

Photo-transect total hard coral is about 0.2–0.3. Space limitation therefore caps fast-coral recovery. In the inherited-composition runs this was masked, because first waves of 100–190 adults/ha also consumed massives (per-capita demand spills over to slow prey) and cleared space. CoralBlox group 5 contributes: its large-massive linear-extension scale is 1.16, against 0.71 in group 1.

**Decision.** Nothing is promoted.
- **Grace switch.** `COTS_STARVATION_GRACE_YEARS` stays opt-in. One year exceeds the evidence (about 6 months) and does not fix the gap.
- **Food scale.** On Lizard and MacGillivray the "food-scale" gap is a coral **space and massive-growth** problem, not a starvation-timing one.
- **Working values.** The best configuration remains `obs_m100_q30` / `acrpoc_g0` (loss 12.18).

**Next gate (proposed).** A coral-side calibration before further COTS mechanisms:
- **Target.** Fit slow-coral (massive) growth and mortality, and check fast-coral recruitment, against LTMP photo-transect group cover (fast, slow, total; 1993–2024) and manta total cover. Work per calibration reef, under group 5.
- **COTS forcing.** Run the coral fit with COTS predation **forced from observed CPUE** (adults = CPUE / q), so it is independent of the COTS model.
- **Then.** Rerun the leading COTS configuration on the refitted coral.


## Preregistered coral-side calibration with observed COTS forcing (2026-10-12)

Status: **approved by the project lead (2026-10-12); preregistered before implementation.** No default, domain-package or objective change for COTS.

### Motivation

[Run 15](runs/20261012T_15_food_grace_1991_v1/metadata.toml) showed that on Lizard and MacGillivray modelled slow corals (massives) grow from about 0.2–0.3 to 0.6–0.7 of habitable area, while LTMP photo-transect total hard coral is about 0.2–0.3. Space limitation then caps fast coral at about 0.04, below the H1 threshold, so no second wave can build. CoralBlox group 5 (the project lead's choice) gives large massives a linear-extension scale of 1.16, against 0.71 in group 1. Coral dynamics must be constrained by coral data before further COTS mechanisms are judged, and **without depending on the COTS model**.

### Forcing (new derived input and opt-in path)

**Forcing table.** `sandbox/calibration/build_observed_cots_forcing.py` (checked in) builds `sandbox/data/derived/observed_cots_forcing_1991_2024/observed_adults_ha.csv`, a site × year table:
- **Method.** For each year 1991–2024, LTMP manta COTS/tow (sum(COTS)/sum(tows) per reef-year) is interpolated with the same effort-aware IDW, crosswalk, 111 km radius and 3-donor rule as the 1991 surface.
- **Gaps.** Site-years without enough donors (the biennial gaps 2012, 2014, 2016 and 2018) are filled by linear interpolation in time.
- **Units.** Adults = CPUE / 0.015.
- **Check.** It reproduces the observed calibration-reef CPUE closely (Lizard 1998: 1.73 vs 1.92; 2013: 0.71 vs 0.74).

**Forced grazing.** New COTSMod `forced_grazing!` applies the per-capita, preferred-first grazing of the working model (10 m² per adult per year ÷ site `k`) for prescribed adults, with no population dynamics. A unit test shows it equals `cots_mortality!` when the dynamic adults equal the forcing.

**ADRIA switch.** `COTS_FORCED_ADULTS_CSV` applies forced grazing each year at the point where COTS predation occurs, and requires `ADRIA_COTS_ENABLED=false`.

### Calibrated quantities

The linear-extension scale multipliers apply to every RME calibration group present (5, 1, 4):
- **`s`, slow corals** (small and large massives): `s ∈ {0.40, 0.55, 0.70, 0.85, 1.00}`;
- **`f`, fast corals** (groups 1–3): `f ∈ {1.0, 1.2, 1.4}`.

Mortality scales, growth acceleration, recruitment, DHW, cyclones and the initial state are unchanged. The initial state is the observed 1993–95 photo-transect composition (*Acropora* + Pocilloporidae share) with inherited 2023 totals; this is the leading-candidate composition.

### Frozen coral objective

Observations:
- LTMP photo transects (9 m, `COMPOSITION`): fast = *Acropora* + Pocilloporidae, slow = all other hard coral, as fractions of benthos. They are compared with model fast and slow cover as fractions of habitable area on Lizard, MacGillivray and North Direction.
- LTMP manta 9 m hard coral, compared with model total site-area cover on all four reefs.

All observed years 1992–2024 are used.

`J = mean_reefs RMSE_fast + mean_reefs RMSE_slow + mean_reefs RMSE_manta_total`

- **Selection.** `(s*, f*)` is the forced arm with minimum `J`. A selected value on the grid boundary is reported as unidentified within the bound. The plateau is the set of arms within 10% of the minimum.
- **Secondary, not used for selection:**
  - fast-cover recovery slopes 2017–24 and 2001–13 against observed;
  - peak slow cover;
  - fast cover in 2008–13.

### Arms (16)

- **`dynamic_ref`:** the leading COTS configuration (COTS dynamics on). It must replay `acrpoc_g0` of run 15 exactly, and is shown for contrast.
- **`forced_s{040…100}_f{100,120,140}`:** forced observed COTS, COTS dynamics off.

### Expected consequences

- A lower `s` reduces slow overgrowth towards photo-transect totals and frees space, raising fast cover. A higher `f` speeds fast recovery.
- `s*` is expected well below 1.0, and `f*` near 1.0–1.2 if space was the main limit.
- *Risk:* fast-cover recovery stays too slow at any `s`. That would implicate fast-coral recruitment or larval supply, which is outside this grid.

### Next step after selection

Rerun the leading COTS configuration on the refitted coral (`s*`, `f*`), as a separate preregistered screen.

### Failures and rollback

Failures:
- an inexact replay of `dynamic_ref`;
- a forcing table that does not cover every site and year;
- non-finite cover;
- a forced arm whose grazing exceeds the available cover. The function caps it, but this would signal a units error.

Rollback: unset `COTS_FORCED_ADULTS_CSV` (default) and multipliers 1.0.

### Coral calibration decision (2026-10-13)

[v2 run](runs/20261013T_01_coral_calibration_1991_v2/metadata.toml) and [audit](runs/20261013T_01_coral_calibration_1991_v2/coral_calibration_audit.json).
- **Operational failure.** [v1](runs/20261012T_16_coral_calibration_1991_v1/metadata.toml) failed in all forced arms: `DimVector` year columns in the new forcing loader. It was fixed, a test was added, and the run was repeated as v2.
- **Checks.** `dynamic_ref` replays run 15 `acrpoc_g0` exactly (`0.0`). COTS dynamics are zero in all forced arms. Cover stays finite and in range.
- **Tests.** COTSMod 363/363 pass, including `forced_grazing!`'s equality with the dynamic per-capita path.

| `s` \ `f` | 1.0 | 1.2 | 1.4 |
| --- | --- | --- | --- |
| 0.40 | **0.381** | 0.440 | 0.499 |
| 0.55 | 0.406 | 0.459 | 0.516 |
| 0.70 | 0.442 | 0.488 | 0.536 |
| 0.85 | 0.486 | 0.524 | 0.563 |
| 1.00 | 0.537 | 0.565 | 0.598 |

The table shows `J`. For comparison, `dynamic_ref` has `J = 0.599`: fast RMSE 0.033, slow 0.364, manta 0.202.

**Selection.** `s* = 0.40`, `f* = 1.0`. **Both are on the grid boundary, so `s` and `f` are not identified within the bounds**, as the frozen rule specifies.

**Against the predictions.**
- **Slow corals.** A lower `s` reduces slow overgrowth, as predicted: mean slow cover in 2012 falls from 0.50 to 0.25, with Lizard and North Direction near observed and MacGillivray still about 0.35 against 0.2.
- **Fast corals.** The freed space lets fast coral overshoot. Even at defaults under observed COTS, fast coral 2008–13 is 0.14, against observed 0.035–0.105; at `s = 0.40` it is 0.29. Hence `f* = 1.0` sits at the lower bound.
- **The preregistered risk** (fast recovery too slow at any `s`) did **not** materialise under observed COTS: the 2017–24 slope reaches the observed range at `s ≤ 0.55`.
- **Disturbance response (not preregistered).** Observed manta cover stays at 0.07–0.10 from 2015 to 2019 after Nathan and the 2016–17 bleaching, while model cover rebounds within 1–2 years at every setting.

**Decision.** Nothing is promoted, and coral parameters are unchanged.

**Interpretation.** With observed COTS, the group-5 coral model on these reefs grows **both** coral types too fast, and loses too little to the 2015–17 disturbances. The low fast cover in the dynamic-COTS runs, which matched transects, came from massive crowding and trough grazing, not from realistic growth. A growth-only fit would let reduced growth stand in for missing disturbance mortality.

**Next gate (project lead to choose):**
- (a) extend the growth grid to `f` 0.6–1.0 and `s` 0.2–0.4;
- (b) add a disturbance-mortality sensitivity (bleaching via the existing scalar, and cyclone), scored on the 2015–19 decline and recovery;
- (c) a joint bounded fit of growth and disturbance mortality, holding out years, for example scoring 1992–2012 and validating 2013–24.


## Calibration-group provenance check (2026-10-13)

Requested by the project lead, because the domain has been updated in recent months.

- **Source of the groups in use.** `sandbox/data/rme_ml_2025_06_05/data_files/region/reefmod_gbr.gpkg`, column `CB_CALIB_GROUPS` (runner `RME_GBR_LAYER`).
- **Calibration-time canonical file.** CoralBlox calibration used `CoralBlox-params-calibration/datasets/spatial_data/rrap_canonical_2025-07-15-T10-48-29.gpkg`. Its `CB_CALIB_GROUPS` and `GBRMPA_BIOREGION` **agree with the RME layer on all 3,806 reefs** (100%).
  - Domain reefs: 58 in group 5 (GBRMPA bioregion 26; all four calibration reefs), 44 in group 1 (bioregion 3) and 11 in group 4 (bioregion 23).
  - There are 12 groups.
- **Method.** `gbrrestoration/canonical-reefs` `src/14_add_cb_calib_groups.jl`:
  - each reef takes the GBRMPA bioregion with the largest overlap, or the nearest;
  - bioregions with fewer than 4 LTMP manta locations having ≥ 4 observations spanning ≥ 10 years since 2008 are merged into the nearest data-rich bioregion (10% quantile of reef distances);
  - group 27 is split at Feather Reef's latitude;
  - groups are renumbered 1..N.
- **Upstream changes.** None to this script since 2025-07-15. Upstream `main` (`a258d69`, 2026-02-09) changed only the README. The local clone is at `56bd7e7` (2025-08-26); Zenodo 10.5281/zenodo.16945486 is a source snapshot from 2025-08-26.
- **Locally regenerated outputs.** The files in `canonical-reefs/output/` lack `CB_CALIB_GROUPS`; they are earlier or partial pipeline runs and are not used.
- **ADRIA parameter tables.** `linear_extension_group_scale_factors`, `mb_rate_group_scale_factors` and the growth-acceleration vectors were set in ADRIA commit `f931cfa42` (2025-10-08, P. Ribeiro de Almeida, "Update calibrated parameter base values"). The group factors were introduced in `b00df9ceb` (2025-07-01).
- **Calibration context.** CoralBlox calibration was run over **2008–2022 with historical DHW and observed LTMP storm/COTS mortality forced** (`historical_disturbance_mortality_rates`, published 2025-11-05). "Unknown", disease and bleaching-only events were excluded.

**Conclusion.** The group assignment is consistent and current. The group-5 parameters, however, were fitted with observed acute disturbance losses imposed, which are much larger than the RRN-hours dose-response used here. That is consistent with the coral overgrowth found in the [coral calibration](runs/20261013T_01_coral_calibration_1991_v2/metadata.toml).

## Preregistered joint coral fit: growth × cyclone strength × background mortality (2026-10-13)

Status: **directed by the project lead (2026-10-13); preregistered before implementation.** No default, domain-package or COTS change.

### Project lead's observations

- Cyclone mortality is not strong enough: the post-2015 trough is not deep enough, while post-2016 recovery looks about right.
- LTMP records other disturbances (disease, unknown and others). Mellin et al. (2019) used their incidence to derive a background mortality needed to contain coral over 2000–2010.

### Evidence

**Observed storm losses.** From the CoralBlox disturbance dataset (LTMP, per functional group), survival in the disturbance year was:

| Reef, year | Survival by group |
| --- | --- |
| North Direction 2015 (Ita/Nathan) | 0.05 *Acropora*, 0.25 corymbose non-*Acropora*, 0.24 small massives, 0.50 large massives |
| Lizard 2015 | non-*Acropora* 0.35–0.51 |

The current dose-response gives total mortality of 0.20–0.25 for Nathan and 0.06 for Ita on these reefs. Scaling the dose-response coefficient by κ, so that mortality = `1 − (1 − m)^κ` (equivalent to multiplying the −0.00444 per hour coefficient), the observed losses imply **κ ≈ 2.5–5**.

**Other-disturbance incidence.** `CoralBlox-params-calibration/datasets/ltmp_data/disturbances.nc` covers 204 LTMP reefs over 2008–2022:

| Disturbance type | Reef-years flagged |
| --- | --- |
| Unknown | 51 (1.7%/yr) |
| Disease | 1 |
| Bleaching | 51 |
| Storm | 98 |
| COTS | 61 |
| Multiple | 28 |

North Direction has "Unknown" in 2008, 2010 and 2022, and MacGillivray in 2008. Multiplying an "Unknown/Disease" incidence of about 2–20% per reef-year by a plausible impact of 10–30% cover loss gives a chronic background loss of about **0.3–6% per year**.

**Mellin et al. 2019** is the lead's reference for this approach. If it reports a specific rate, it should replace this bound.

### Mechanism (runner only; no ADRIA code change)

The domain's cyclone-mortality cube (site × year, all groups) is replaced per arm by

`m'(site, year) = 1 − (1 − m)^κ · (1 − μ_bg)`

- `κ` scales cyclone (RRN-hours) mortality.
- `μ_bg` is a constant annual background mortality (fraction of cover), applied to all groups and sizes every year.

Both enter ADRIA through the existing cyclone mortality step, before COTS grazing. With κ = 1 and μ_bg = 0 the original cube object is used, so the reference replays exactly.

### Grid (54 forced arms)

All arms force COTS grazing from observed CPUE / q (`COTS_FORCED_ADULTS_CSV`, COTS dynamics off) and use the observed 1993–95 composition with RME group 5.

| Factor | Levels |
| --- | --- |
| Slow (massive) growth `s` | 0.55, 0.75, 1.00 |
| Fast growth `f` | 0.8, 1.0 |
| Cyclone strength κ | 1, 2.5, 5 |
| Background mortality μ_bg | 0, 0.03, 0.06 per year |

**Reference.** `forced_s100_f100_k100_b000` must replay `forced_s100_f100` of [coral calibration v2](runs/20261013T_01_coral_calibration_1991_v2/metadata.toml) exactly.

### Frozen objective and checks

- **Objective.** Unchanged `J` = mean RMSE fast + mean RMSE slow (photo transects, habitable fraction) + mean RMSE total (manta, site area), over 1992–2024.
- **Selection.** Minimum `J`. Each parameter at a grid boundary is reported as unidentified within the bound. The plateau is the set within 10% of the minimum.
- **Held-out reefs.** Re-select on each set of 3 of the 4 reefs and report the held-out reef's error.
- **Secondary reports:**
  - `J` split by period (1992–2012 and 2013–2024);
  - 2017–19 trough depth (manta) against observed 0.07–0.10;
  - post-2016 recovery slope;
  - 2008–13 fast cover;
  - peak slow cover.

### Expected consequences

- κ ≈ 2.5–5 deepens the 2015–17 trough to observed levels without changing pre-2014 years.
- μ_bg > 0 lowers both coral types through 2000–2012, so slow overgrowth shrinks at a higher `s` than in the growth-only fit.
- `f < 1` controls the fast overshoot.

If the optimum still sits at `s` 0.55 or `f` 0.8, growth remains biased high even with disturbances added.

### Failures and rollback

Failures:
- an inexact reference replay;
- a combined mortality outside [0, 1];
- non-finite cover;
- COTS dynamics active in a forced arm.

Rollback: κ = 1, μ_bg = 0 (the domain cube), and multipliers 1.

**Clarifications by the project lead (2026-10-13, before results):**
- **Cyclones are calibrated on modelled outputs, not forced with observed losses.** κ scales the modelled RRN-hours dose-response and is selected by the fit of modelled coral to LTMP. The observed LTMP storm survival values only motivated the preregistered κ range; they are never imposed.
- **The Mellin et al. (2019) background mortality is a separate component from the CoralBlox storm/COTS disturbance dataset**, and is to be included as its own term (μ_bg). The grid value μ_bg is a placeholder bound from LTMP Unknown/Disease incidence; it should be replaced by, or checked against, the Mellin et al. derivation (disturbance incidence × impact) when that is available.

## Preregistered Mellin background-mortality follow-up (2026-10-13)

Status: **directed by the project lead (2026-10-13); preregistered before the joint-fit results were seen.**

### Source

The lead's COTS_Model, `github.com/sammatthews990/COTS_Model`:
- `COTSMod.R`, `COTSModel_COTSFunctions.R::doCoralDisturbance`;
- parameters in `RData/COTSMod_bckp.RData`.

Copies are kept in git-ignored `sandbox/data/external/COTS_Model/`. Hashes: RData `14a0628e…0631`, `COTSMod.R` `221bdad3…403c`.

How disease and unknown events act:
- **Effect.** Each event acts on log cover, `log(C) += B`, with B_disease −0.1681 (sd 0.0303) and B_unknown −0.2681 (sd 0.0331). That is 15.5% and 23.5% cover loss per event.
- **Incidence.** Events are simulated at 1% per cell-year each (`sampling::srswor(round(n*0.01), n)`) and replaced by LTMP records where available. Observed LTMP incidence in the RData is lower: 0.32% (disease) and 0.53% (unknown) of 38,971–39,515 observed cell-years.

### Derived input

The table is built in two checked-in steps:
1. `export_mellin_disturbance.R` exports the coefficients and LTMP reef-year records (1996–2017).
2. `build_mellin_background_mortality.jl` writes `sandbox/data/derived/mellin_background_mortality_1991_2024/background_mortality.csv` (site × year).

How each reef-year gets its value:
- **Recorded reef-years.** `m = 1 − exp(B_d·disease + B_u·unknown)`. Sources are the CoralBlox `disturbances.nc` flags (2008–2022, preferred) and the COTS_Model LTMP tables (1996–2017).
- **Unrecorded reef-years.** The expected loss at 1% incidence each: `m = 0.0039`/yr.
- **Recorded events in the domain:** 17 reef-years. They include North Direction 2008, 2010 and 2022 (unknown) and MacGillivray 2008 (unknown). The COTS_Model table dates MacGillivray's event to 2009; the newer CoralBlox source is preferred.

### Arms (18)

| Factor | Levels |
| --- | --- |
| `s` | 0.55, 0.75, 1.00 |
| `f` | 0.8, 1.0 |
| κ | 1, 2.5, 5 |
| Background | Mellin table |

Combined mortality is `1 − (1 − m_cyc)^κ · (1 − m_bg[site, year])`. COTS are forced from observed CPUE / q.

**Reading.** The arms are compared with the joint run's matching `μ = 0` arms (same s, f, κ) using the frozen `J`, and the full model is selected across both runs by `J`. There is no replay arm, because every arm adds the Mellin term. The `μ = 0` arms of the joint run serve as the paired controls.

**Failures:**
- a background table that does not cover every site and year;
- mortality outside [0, 1);
- COTS dynamics active.

### Joint coral and Mellin decision (2026-10-13)

**Runs:**
- **Joint fit:** split across [part 1](runs/20261013T_02_joint_coral_1991_v1), [part 2](runs/20261013T_04_joint_coral_1991_v1_part2) and [part 3](runs/20261013T_05_joint_coral_1991_v1_part3), audited together; the audit lives in part 3.
- **Mellin follow-up:** [run](runs/20261013T_03_mellin_coral_1991_v1).

**Operational disclosures:**
- **Slowdown.** Per-arm time grew within a long-running process. Parts 1 and 2 were stopped and the remaining arms completed in fresh processes via the new `OWEN_1991_ARM_FILTER`. The audit checks that every arm appears exactly once across directories.
- **Metadata.** Parts 1 and 2 lack `metadata.toml`. Their arm definitions and code paths are identical to part 3's metadata. Part 1 ran a runner version without the `MELLIN`/filter additions, but those are gated and do not touch `joint_coral`.
- **Audit bug.** The first audit attempt failed on a relative-path bug after writing its CSV and plots. Those partial outputs were deleted and the audit rerun.
- **Replay.** The reference arm replays exactly.

**Joint fit.** The best `J` is 0.240 at `s = 0.55`, `f = 0.8`, κ = 2.5, μ = 0.06, against 0.537 for the reference.

| Profile | Best `J` by level |
| --- | --- |
| μ = 0 / 0.03 / 0.06 | 0.317 / 0.265 / **0.240** |
| `f` = 0.8 / 1.0 | **0.240** / 0.281 |
| κ = 1 / 2.5 / 5 | 0.246 / **0.240** / 0.243 |
| `s` = 0.55 / 0.75 / 1.00 | **0.240** / 0.241 / 0.261 |

- **Bounds.** μ and `f` are at their bounds; `s` is nearly flat; κ is weakly identified.
- **Held-out reefs.** Leave-one-reef-out selections are nearby arms, and every held-out error is about half the reference's.
- **Trough depth.** At κ = 5 the 2017–19 total is 0.10, against observed 0.071.

**Mellin follow-up.** The Mellin background improves `J` against μ = 0 in all 18 pairs (ΔJ −0.02 to −0.05). The best is 0.298 (κ = 5), against 0.240 with the constant μ = 0.06.

**Against the predictions:**
- κ ≈ 2.5–5 deepens the post-2015 trough, as predicted, but κ is weakly identified.
- μ > 0 shrinks slow overgrowth, so `s` can stay higher, as predicted.
- `f < 1` is needed, as predicted.
- **Not predicted:** a constant μ removes cover from all groups, and fast corals refill the space. Fast cover then overshoots photo transects 2–4× in 2005–13 even at `f = 0.8`.

**Decision.** No coral parameters are promoted.

**Interpretation.** On these reefs the model needs more coral loss than storms, COTS and the recorded disease/unknown events supply. That missing chronic loss is about 6% per year against the Mellin evidence of 0.4% per year, and it falls disproportionately on fast corals.

**Next gate (project lead to choose):**
- (a) DHW bleaching-response sensitivity: ADRIA's bleaching mortality, scored on the 1998, 2002, 2016, 2017 and 2020–24 years, with the group-specific susceptibility already in ADRIA;
- (b) a fast-coral-specific chronic loss or recruitment limit, with evidence needed;
- (c) an extended grid (`f` 0.5–0.8; μ fixed at the Mellin value) to test whether growth bias alone explains the fast excess.

Constant μ = 0.06 should not be adopted without a mechanism, because it exceeds the incidence evidence.


## Preregistered bleaching-response sensitivity (2026-10-13)

Status: **directed by the project lead (2026-10-13); preregistered before implementation.**

### Lead's principles

- **Leave CoralBlox's calibrated growth untouched as far as possible.** Those rates were calibrated to LTMP data (2008–2022, with observed storm/COTS losses imposed). Growth multipliers stay `s = f = 1`.
- **Background loss.** The lead considers background loss and cyclone strength between the joint-fit best and the Mellin evidence "pretty good".

### Evidence and expectation

- **ADRIA bleaching.** Mortality is the fraction of each group/size class whose heat tolerance (truncated normal; lower bound `HEAT_LB = 1` DHW, or the previous bleaching DHW) falls below the depth-adjusted DHW. `ADRIA_DEBUG_BLEACHING_SCALAR` scales the resulting loss linearly (β).
- **DHW history on the calibration reefs** (v0.3 historical; mean of sites, depth_med 5 m):
  - 2016: about 8.6;
  - 2017: about 8.9;
  - 2020: 5.1–6.6;
  - 2024: about 6;
  - 2021/2022: about 4;
  - 1994: about 3.6–4;
  - 2002/2004/2010: 2–2.7;
  - **1998 below 2.**
- **Expectation.** Bleaching can deepen the 2016–17 trough and slow the 2020s recovery. It **cannot** explain the 2005–13 fast-coral overshoot, because there is no material heat stress there. If the fit still needs a large chronic loss, that loss is not bleaching.

### Arms (12)

COTS are forced from observed CPUE / q, with the observed composition and RME group 5. Growth is fixed at `s = f = 1`.

| Factor | Levels |
| --- | --- |
| Background | `mel` (Mellin table), `mel3` (Mellin plus a constant 3% per year chronic loss, midway between Mellin's 0.4% and the joint best of 6%) |
| Cyclone strength κ | 2.5, 5 |
| Bleaching scalar β | 1.0, 1.5, 2.5 |

Combined non-bleaching mortality is `1 − (1 − m_cyc)^κ · (1 − m_mellin) · (1 − μ)`.

**Reference.** `bleach_mel_k250_b100` must replay `mellin_s100_f100_k250` of [run 03](runs/20261013T_03_mellin_coral_1991_v1/metadata.toml) exactly.

**Objective and reports.**
- **Objective.** The frozen `J`, with selection by minimum `J`; bounds are flagged.
- **Secondary reports:**
  - `J` by period (1992–2012, 2013–2024);
  - 2017–19 total trough against manta (0.071);
  - 2020–24 fast-recovery slope against photo transects;
  - 2005–13 fast overshoot.

**Failures.**
- inexact reference replay;
- mortality outside [0, 1];
- COTS dynamics active.

### Bleaching decision and runtime fix (2026-10-13)

**Runtime fix.** The per-arm slowdown in long runner processes came from `set_initial_cover!`. Its use of `collect(cube.species)` nested DimensionalData lookups each arm, so DataCube rebuild time grew geometrically (diagnostic runs 06–08; standalone reproduction). It is fixed with plain label vectors captured once. Results were unaffected and replays are exact. Runner arm timing now prints `run_scenario`, setup (with checkpoints) and post-processing times.

**Bleaching screen** ([run](runs/20261013T_09_bleaching_sensitivity_1991_v1/metadata.toml), [audit](runs/20261013T_09_bleaching_sensitivity_1991_v1/bleaching_audit.json)). The reference replays exactly. **No arm is identified.**

| Arm | `J` | 2017–19 total (obs 0.071) | Fast 2005–13 (obs 0.065) | Fast slope 2020–24 (obs 0.015) |
| --- | --- | --- | --- | --- |
| Mellin, κ 2.5, β 1 (ref) | 0.456 | 0.284 | 0.140 | 0.007 |
| Mellin, κ 5, β 1 | 0.427 | 0.219 | 0.141 | 0.013 |
| Mellin + 3%, κ 5, β 1 | 0.341 | 0.157 | 0.155 | 0.018 |
| Mellin + 3%, κ 5, β 1.5 | 0.316 | 0.126 | 0.152 | 0.007 |
| **Mellin + 3%, κ 2.5, β 2.5 (min `J`)** | **0.284** | 0.130 | 0.143 | −0.001 |
| Mellin + 3%, κ 5, β 2.5 | 0.286 | 0.086 | 0.144 | 0.000 |

Against the predictions:
- **Trough.** Bleaching deepens the 2016–17 trough, as expected.
- **Pre-2015 overshoot.** It cannot fix the 2005–13 fast overshoot, as expected.
- **Recovery.** At β 2.5 it removes fast coral after 2020 and in 1994, which contradicts observed recovery.

**Decision.** Nothing is promoted.

**Recommended working coral configuration** for returning to COTS dynamics (**project lead to confirm**), consistent with the lead's principle of keeping CoralBlox growth untouched:

| Component | Setting |
| --- | --- |
| Growth | `s = f = 1` |
| Background | Mellin + 3% per year (between Mellin's 0.4% and the joint best of 6%) |
| Cyclone strength | κ 5 |
| Bleaching | β 1–1.5 |

**Known residual.** Fast coral overshoots observations about 2× in 2005–13 under observed COTS. Possible reasons: COTS undercounted by manta in troughs, fast-coral recruitment, or site representativeness. This matters for second-wave food supply.


## Preregistered COTS re-baseline on the working coral configuration (2026-10-13)

Status: **directed by the project lead (2026-10-13); preregistered before implementation.** Runner-only screen; no ADRIA, COTSMod, default, domain-package, objective or observation-model change.

### Project lead's decisions

- **Working coral configuration confirmed:**

| Component | Setting |
| --- | --- |
| Growth | CoralBlox calibrated growth untouched (`s = f = 1`; RME groups, group 5 on the calibration reefs) |
| Background | Mellin table + 3%/yr chronic loss |
| Cyclone strength | κ = 5 |
| Bleaching scalar | β = 1–1.5 (a bracket, not a point value) |

- **Candidate to re-baseline:** run 15 `acrpoc_g0` (= run 14 `obs_m100_q30`). This is the working model ("condition, no switching"), Ext×5 and predator pit `Q = 30` (both labelled assumptions), gradual Allee, and the observed 1993–95 composition.
- **`Q` must be re-checked**, because the coral configuration changes coral supply.
- **Manta detectability (recorded for a later screen, not used here).** The lead states that manta tows under-detect COTS at low density. Evidence they supplied: manta tows detect about 20% of the COTS counted on SCUBA searches, and under 5% of the total population (Fernandes et al. 1990). Diver methods such as culling and SALAD reach about 50–82% for adults, falling markedly for small adults and juveniles (MacNeil et al. 2016a; Pratchett & Caballes 2025; Pratchett et al. 2026). The observation model stays frozen at `q = 0.015` in this screen, so that the re-baseline is comparable with run 15.

### Hypothesis and mechanism

**Hypothesis.** The run-15 candidate's small second waves on Lizard and MacGillivray come partly from slow-coral (massive) overgrowth. That overgrowth exists because the coral model loses too little coral to disturbance. With the working coral configuration:
- slow cover is held down;
- fast coral can recover above the starvation threshold θ = 0.05 between waves;
- larger second waves become possible.

**Change.** No equation changes. Each arm replaces the domain cyclone-mortality cube (site × year, all groups and sizes; applied before COTS grazing) with

`m' = 1 − (1 − m_cyclone)^κ · (1 − m_mellin[site, year]) · (1 − μ)`

and sets `ADRIA_DEBUG_BLEACHING_SCALAR = β`. This is the bleaching screen's construction, now with COTS dynamics on instead of forced COTS.

### Arms (13)

All arms use the run-15 candidate, `q = 0.015`, `k_unit = 19.61`, Allee scale 3, `ADRIA_RNG_SEED` 90972102 and Owen connectivity seed 20260930.

| Arm | κ | Background | β | `Q` | Role |
| --- | --- | --- | --- | --- | --- |
| `ref_q30` | 1 | none | 1 (unset) | 30 | must replay run 15 `acrpoc_g0` exactly |
| `k500_q30` | 5 | none | 1 | 30 | component attribution (cyclones only) |
| `mel3_q30` | 1 | Mellin + 3% | 1 | 30 | component attribution (background only) |
| `wc_b100_q{10,20,30,45,70}` | 5 | Mellin + 3% | 1.0 | grid | working coral, `Q` re-check |
| `wc_b150_q{10,20,30,45,70}` | 5 | Mellin + 3% | 1.5 | grid | working coral, `Q` re-check |

The `Q` grid extends one step below the previous lower bound of 20, because run 14 selected `Q = 20` at the bound for one growth setting.

### Expected consequences

Forced-COTS runs (run 09, working configuration) give directional expectations:
- **First wave.** Pre-outbreak fast cover is about 10–20% lower at β 1 (3%/yr chronic loss). At β 1.5 it is about 30–40% lower, because the 1994 heat stress (DHW about 4) removes fast coral. First peaks on Lizard and MacGillivray therefore fall further below observed (risk: losing first peaks at β 1.5). Timing may move up to 1 year later.
- **Between waves.** Slow cover stays at about 0.3–0.4 instead of 0.6–0.7, so fast coral on Lizard and MacGillivray recovers above θ = 0.05 in 2003–12. Second waves should be larger and more frequent.
- **Troughs.** More inter-wave food raises trough adults (risk: back above 3–4 adults/ha). A higher `Q` may then be selected.
- **κ 5.** Heavier 2014–15 cyclone losses truncate food after 2014. Second waves that peak before 2014 should be little affected, and post-2015 COTS should starve sooner.
- **Components.** `mel3_q30` should carry most of the second-wave effect, and `k500_q30` should change only post-2014 years.

### Selection rules (frozen before results)

- **`Q` per β.** The minimum summed `smoothed_3y` loss (frozen `cots_cycle_score`, `q = 0.015`) over the 4 reefs, requiring a first peak on every reef; ties go to the smaller `Q`. A selected `Q` of 10 or 70 is reported as unidentified within the bounds.
- **β for later screens.** β is a coral parameter, so it is **not** chosen by the COTS loss. The new baseline for screens 2–4 uses the β with the lower frozen coral `J` (fast + slow photo-transect RMSE + manta total RMSE, 1992–2024), computed in these dynamic-COTS runs at each β's selected `Q`. COTS metrics are reported at both β.
- **Qualitative gate (unchanged).** Two peaks on each of Lizard, MacGillivray and North Direction; 7/7 observed peaks matched; worst trough-to-smaller-peak ratio ≤ 0.10.

### Report

- the gate, summed loss and matched peaks;
- per reef: first-peak year and height (Lizard and MacGillivray against observed 1.92 / 0.86 COTS/tow), second-peak height, and 2010–17 maximum against observed;
- trough minimum adults/ha;
- fast cover 2000–13 against θ, and slow cover in 2003 and 2012;
- coral `J` with its components;
- juvenile predation removed;
- component-arm attribution.

### Failures and rollback

**Failures:**
- an inexact reference replay;
- external supply other than Ext×5;
- combined mortality outside [0, 1];
- non-finite states;
- failed arms (logged with their error).

**Interpretation commitments.**
- If no `Q` gives second peaks on Lizard and MacGillivray, the second-wave gap is not explained by coral disturbance alone. Screen 3 (second-wave food) then starts from this baseline.
- If first peaks are lost on any calibration reef at a β, that β is reported as incompatible with the current first-wave mechanism. It is not hidden by retuning `Q`.

**Rollback.** κ 1, no background, β unset (1.0), `Q = 30`: the run-15 configuration.

### Re-baseline decision (2026-10-13)

[Run](runs/20261013T_11_cots_rebaseline_1991_v1/metadata.toml), [audit](runs/20261013T_11_cots_rebaseline_1991_v1/cots_rebaseline_audit.json) and [plot](runs/20261013T_11_cots_rebaseline_1991_v1/cots_rebaseline.png).

**Operational failure.** The first attempt, [`20261013T_10_cots_rebaseline_1991_v1`](runs/20261013T_10_cots_rebaseline_1991_v1/FAILED.txt), stopped before any arm ran: the runner's treatment list mixed `,` and `;` in one array literal. This is a Julia *lowering* error, so `test_calibration_scripts_parse.jl` (which only calls `Meta.parseall`) passed. The list was fixed and checked with `Meta.lower` on every top-level expression. The run was repeated as run 11; the empty directory is kept with a `FAILED.txt` note.

**Checks.**
- `ref_q30` replays run 15 `acrpoc_g0` exactly (`0.0` in trajectories, flows and mechanisms).
- External immigration is identical to Ext×5 in all 13 arms (`0.0`).
- States are finite and mortality stays in [0, 1].

**Result.** **No arm passes the gate.** Every working-coral arm has two peaks on all three two-wave reefs, but matches 6/7 peaks and has a worst trough ratio of 0.24–0.30 (gate ≤ 0.10, failing on Lizard).

| Arm | Summed loss | Matched | Two-wave reefs | Mean first peak (yr / COTS/tow) | Lizard / MacG first peak (obs 1.92 / 0.86) | Lizard / MacG 2010–17 max (obs 0.74 / 0.67) | Trough min (adults/ha) | Coral `J` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ref_q30 (run 15) | 12.18 | 5/7 | 2 | 1998.3 / 0.66 | 0.44 / 0.45 | 0.17 / 0.14 | 2.9 | 0.599 |
| k500_q30 (κ 5 only) | 13.90 | 4/7 | 1 | 1998.3 / 0.66 | 0.43 / 0.45 | 0.17 / 0.15 | 2.9 | 0.510 |
| mel3_q30 (Mellin + 3% only) | **8.75** | 6/7 | 3 | 1998.3 / 0.64 | 0.42 / 0.44 | 0.20 / 0.41 | 2.9 | 0.415 |
| wc_b100_q10 / q20 / q30 / q45 / q70 | 9.23 / **9.12** / 9.53 / 9.89 / 10.42 | 6/7 | 3 | 1997.0–2000.0 / 0.60–0.64 | 0.41 / 0.37 (q20) | 0.16 / 0.38 (q20) | 2.8–3.0 | 0.365–0.371 |
| wc_b150_q10 / q20 / q30 / q45 / q70 | 9.16 / **9.06** / 9.46 / 9.54 / 10.04 | 6/7 | 3 | 1997.0–1999.8 / 0.46–0.56 | 0.35 / 0.35 (q20) | 0.16 / 0.40 (q20) | 2.6–2.8 | 0.354–0.369 |

**Selection (frozen rules).**
- `Q` = 20 at both β, an interior value. The differences between `Q` values are small (about 0.4 in loss between 20 and 30).
- β = 1.5 by coral `J` at `Q = 20`: 0.364 against 0.369. **The margin (0.005) is negligible.**
- **New baseline for screens 2–4: `wc_b150_q20`** (summed loss 9.06, against 12.18 for run 15).

**Against the predictions:**
- **First wave.** Lower, as predicted, and more so at β 1.5 (the 1994 heat stress): Lizard 0.35 and MacGillivray 0.35 COTS/tow against observed 1.92 / 0.86. North Direction (0.62 vs 0.56) and Eyrie (0.75, against an observed peak of about 0.13) remain high. Timing is unchanged or earlier (1997; observed 1996–98).
- **Slow coral.** Lower, as predicted: 2012 slow cover 0.43–0.48 against 0.62–0.70 in run 15. It remains 2–3× the photo transects (0.15–0.24).
- **Second waves. Partly as predicted.**
  - North Direction gains its 2012 wave (0.83 vs observed 0.71), worth −2.6 in loss.
  - MacGillivray's 2013 wave rises from 0.14 to 0.40 (observed 0.67).
  - **Lizard is unchanged** (0.16 vs 0.74).
- **The predicted mechanism did not occur.** Freeing space did **not** raise inter-wave fast coral. Lizard and MacGillivray fast cover in 2000–13 stays at 0.026–0.029, against 0.026–0.030 in run 15 (observed 0.047 / 0.027); it never exceeds θ on Lizard. The larger MacGillivray and North Direction waves build on 0.04–0.05 fast cover, where the food-survival factor already reaches 0.9–1.0.
- **Troughs.** Unchanged at 2.6–3.0 adults/ha, inside the observed range. COTS do not fall away between waves: adults stay at 2–12 adults/ha and subadults at 11–130/ha through 1999–2012.
- **Components.** The background term (`mel3`) carries all the gain (−3.4). κ 5 alone is *worse* (+1.7), because the 2014–15 cyclones remove MacGillivray's late food and its second peak (MacGillivray loss 1.8 → 3.8). In combination the background effect dominates.

**Diagnosis.** Under the same working coral configuration, inter-wave fast coral is:
- about **2× too high** with COTS forced from manta CPUE (run 09: 0.14–0.16 in 2005–13);
- **too low** with COTS dynamics (0.026–0.029 on Lizard and MacGillivray; North Direction 0.02 vs observed 0.075).

The difference is the trough COTS density. The model keeps 2–12 adults/ha plus a large subadult pool, removing 0.3–1.0% of habitable cover per year from preferred prey, while manta tows read about zero. This places the manta-detectability question (item 5) at the centre of both the coral residual and the second-wave food. It does not change the frozen observation model: that needs the project lead's approval.

**Decision.**
- Nothing is promoted.
- The working coral configuration is adopted for COTS screens as the lead directed: Mellin + 3%, κ 5, growth untouched, with β 1.5 and `Q = 20` as the new baseline.
- `Q = 20` and the Ext×5 boundary remain labelled assumptions.

**Next (order set by the lead):** (2) first-wave height across the plausible composition range on this baseline; (3) second-wave food; (4) the 1999 crash; (5) manta detectability, pending the lead's approval of an observation-model change.


## Preregistered first-wave composition range on the re-baselined model (2026-10-13)

Status: **item 2 of the order set by the project lead (2026-10-13); preregistered before implementation.** Runner-only; no ADRIA, COTSMod, default or objective change.

### Question

Can the plausible range of 1991 fast-coral shares raise first-wave peaks on Lizard and MacGillivray (0.35 / 0.35 COTS/tow at the new baseline, against observed 1.92 / 0.86) without losing the troughs?

### Inputs

The composition levels are those of run 15 (1993–95 photo-transect fast share = (*Acropora* + Pocilloporidae + w × "Other") / all hard coral; totals and within-block structure preserved):

| Level | w | Lizard | MacGillivray | North Direction | Other reefs (mean) |
| --- | --- | --- | --- | --- | --- |
| `acrpoc` | 0 | 0.209 | 0.125 | 0.250 | 0.195 |
| `halfother` | 0.5 | 0.249 | 0.201 | 0.383 | 0.278 |
| `other` | 1 | 0.290 | 0.277 | 0.516 | 0.361 |

### Arms (6)

Composition {`acrpoc`, `halfother`, `other`} × β {1.0, 1.5}, all on the re-baselined model:
- working coral (Mellin + 3%, κ 5, growth untouched);
- `Q = 20`, Ext×5, gradual Allee;
- the working COTS model.

β is crossed because the 1994 heat stress directly lowers pre-outbreak food, which is this screen's subject, and the β choice in run 11 was decided by a negligible coral-`J` margin.

**References.** `wc_b150_q20_acrpoc` and `wc_b100_q20_acrpoc` must replay `wc_b150_q20` and `wc_b100_q20` of [run 11](runs/20261013T_11_cots_rebaseline_1991_v1/metadata.toml) exactly.

### Expected consequences

- **MacGillivray.** The fast share can rise 2.2×, so its first peak should rise towards 0.86.
- **Lizard.** The share rises only 1.4× (0.21 → 0.29). Composition alone is expected **not** to bring Lizard near 1.92. If confirmed, Lizard's first-wave deficit lies elsewhere: supply (connectivity or external), or the food scale.
- **North Direction and Eyrie** (shares up to 0.52 and 0.36). Their first peaks are already at or above observed, and are expected to overshoot further. Troughs are expected to rise, as in run 15 (4.5–6.5 adults/ha at `halfother` / `other`).

### Adoption rule (frozen)

Composition is an uncertain **input**, so it is not chosen by the COTS loss alone. A level is recommended for adoption only if all of the following hold:
1. it raises the Lizard **and** MacGillivray first peaks;
2. the mean trough minimum stays ≤ 4 adults/ha;
3. the summed loss does not worsen by more than 0.5;
4. the mean fast-cover RMSE against the photo transects worsens by ≤ 0.01.

Otherwise `acrpoc` (the direct observation) stays. The project lead decides adoption.

### Report

- per reef: first-peak year and height, second-wave height, and trough minimum;
- summed loss and the gate;
- coral `J` and fast RMSE;
- pre-outbreak (1992–96) fast cover per reef.

### Failures and rollback

**Failures:**
- an inexact replay;
- 1991 totals not preserved;
- external supply other than Ext×5;
- non-finite states.

**Rollback:** the `acrpoc` composition.

### First-wave composition decision (2026-10-13)

[Run](runs/20261013T_12_composition_rebaseline_1991_v1/metadata.toml), [audit](runs/20261013T_12_composition_rebaseline_1991_v1/composition_rebaseline_audit.json) and [plot](runs/20261013T_12_composition_rebaseline_1991_v1/composition_rebaseline.png).

**Checks.**
- Both `acrpoc` arms replay run 11 exactly (`0.0`).
- 1991 totals are identical across compositions, and external supply is identical in all arms.

**Result.** **No arm passes the gate, and no level meets the adoption rule**: every higher share raises the mean trough above 4 adults/ha.

| Arm | Summed loss | Matched | First peak: Lizard / MacG / North Dir. / Eyrie (obs 1.92 / 0.86 / 0.56 / ~0.13) | Trough min (adults/ha) | Fast RMSE | Coral `J` |
| --- | --- | --- | --- | --- | --- | --- |
| b100 `acrpoc` (ref) | 9.12 | 6/7 | 0.41 / 0.37 / 0.76 / 0.91 | 3.0 | 0.028 | 0.369 |
| b100 `halfother` | 9.24 | 6/7 | 0.50 / 0.58 / 2.23 / 1.12 | 4.6 | 0.034 | 0.351 |
| b100 `other` | 16.34 | 4/7 | 0.62 / 0.87 / 3.24 / 1.42 | 5.2 | 0.043 | 0.333 |
| **b150 `acrpoc` (baseline)** | **9.06** | 6/7 | 0.35 / 0.35 / 0.62 / 0.75 | 2.7 | 0.032 | 0.364 |
| b150 `halfother` | 9.32 | 6/7 | 0.42 / 0.56 / 1.59 / 1.01 | 4.3 | 0.035 | 0.331 |
| b150 `other` | 9.58 | 6/7 | 0.51 / 0.76 / 2.39 / 1.17 | 4.2 | 0.040 | 0.331 |

**Against the predictions.** All three were confirmed:
- MacGillivray reaches its observed first peak (0.76–0.87 at `other`).
- Lizard does not: at most 0.62 against 1.92.
- North Direction and Eyrie overshoot strongly (up to 3.2 and 1.4), and troughs rise to 4.2–5.2 adults/ha.

**Diagnosis: the between-reef ranking of first waves is inverted.**
- **Observed ranking:** Lizard 1.92 > MacGillivray 0.86 > North Direction 0.56 > Eyrie about 0.13. Early CPUE (1993–94) is already highest on MacGillivray (0.43) and Lizard (0.23), and low on North Direction (0.06).
- **Model ranking** (baseline): Eyrie 0.75 > North Direction 0.62 > Lizard 0.35 ≈ MacGillivray 0.35.

Two inputs drive this:
1. **The 1991 initial state.** The 1989–91 IDW CPUE surface seeds North Direction at 3.3 adults/ha (0.049 COTS/tow), against 0.5 on Lizard and 0.2 on MacGillivray.
2. **Near-uniform external supply.** First-wave recruitment on all three reefs is 70–95% Ext×5 boundary supply (about 65–100 recruits/ha/yr in 1994–98). Local retention is ≤ 1.4/ha/yr, and internal immigration is 6–50/ha/yr.

No reef-level food or composition change can therefore make Lizard's wave about 3× North Direction's. Lizard lies in the Lizard–Cairns initiation region of the third outbreak wave, which suggests self-seeding and local retention there, but that is a hypothesis to test, not an established mechanism.

**Decision.**
- Nothing is promoted, and the baseline is unchanged: `wc_b150_q20`, `acrpoc`.
- Composition is closed as a lever for Lizard's first wave. MacGillivray's observed peak lies inside the plausible composition range.

**Candidate next tests for first-wave height (project lead to choose):**
- (a) a reef-differentiated 1991 initial state (an IDW surface check against 1992–94 observations);
- (b) local retention / self-recruitment at initiation-region reefs, preregistered with evidence;
- (c) a spatially structured external supply (outbreak-source connectivity) instead of a uniform ×5.


## Preregistered manta detectability: density-dependent observation sensitivity (2026-10-13)

Status: **approved by the project lead (2026-10-13), who brought item 5 forward and chose the detection curve; preregistered before any re-scoring.** This is a *labelled observation-model sensitivity*. The frozen observation model (`q = 0.015`, linear) is **not** replaced, and nothing is promoted. No ADRIA, COTSMod or simulation change: existing trajectories are re-scored.

### Evidence (supplied by the project lead)

- **Manta tows** detect about 20% of the COTS counted on SCUBA searches, and under 5% of the total population present (Fernandes et al. 1990), because COTS are cryptic.
- **Diver methods** (culling, SALAD) detect adults at about 50–82%, with detectability falling markedly for smaller adults and juveniles (MacNeil et al. 2016a; Pratchett & Caballes 2025; Pratchett et al. 2026).
- **Interpretation.** At low density, the COTS present are disproportionately small, cryptic and dispersed, so the manta detection fraction falls. At outbreak density, large, exposed and aggregated animals are seen at close to the conventional rate.

### Observation model (treatment)

`CPUE = q · A · A / (A + A½)`

- `A` is the model adult density (COTS ha⁻¹), and `q = 0.015` COTS/tow per adult/ha (unchanged).
- The relative detection `d(A) = A / (A + A½)` rises from 0 towards 1, so `q` is kept at outbreak densities. At `A = A½`, detection is half the conventional rate.
- **Control.** `A½ = 0` is the frozen linear model and must reproduce the stored scores exactly.
- **Grid.** `A½ ∈ {1, 3, 10, 30}` adults/ha.
  - 1–3 brackets the model trough densities (2–3 adults/ha);
  - 10 sits below the 15 adults/ha outbreak threshold (0.22 COTS/tow);
  - 30 is a deliberately strong upper bracket, which also lowers CPUE at the outbreak threshold to a third.
- **Scale.** This is an observation property, so it applies to the reef-level mean trajectories already written by the runner. Sites within a reef follow near-identical trajectories, so applying it at site level would change little.

### What is re-scored

The frozen `cots_cycle_score` (`smoothed_3y`) and the gate metrics (matched peaks, two-wave reefs, trough-to-smaller-peak ratio) are recomputed on `CPUE(A)`. Arms:
- all 13 arms of [run 11](runs/20261013T_11_cots_rebaseline_1991_v1/metadata.toml) (baseline `wc_b150_q20`, the `Q` grid and the component arms);
- the run-15 candidate `acrpoc_g0` as the legacy reference.

**Not in scope.** Under this observation model, a consistent treatment would also invert the 1991 initial state (`A` from observed CPUE) and the forced-COTS coral tables. Re-scoring holds the simulations fixed, so it tests the **observation side only**. Inverting the initial state is taken up in the reef-differentiated starting-state screen.

### Expected consequences

- **Troughs.** Simulated trough CPUE falls sharply. At `A = 3`, `A½ = 10` gives 0.010 COTS/tow against 0.045 linear, so the trough-to-peak ratio falls (towards the ≤ 0.10 gate) and normalised peak prominence rises.
- **The 1999 drop looks more abrupt** (item 4). A decline from about 25 to 10 adults/ha becomes 0.27 → 0.08 COTS/tow at `A½ = 10`.
- **Peak heights fall modestly** at `A½ ≤ 10` (Lizard 23 adults/ha: −30% at `A½ = 10`). This **worsens** the first-wave height deficit, so the re-score cannot by itself fix item 2.
- **`Q`.** The selected `Q` may move, because trough misfit is penalised less.

### Report

Per `A½`: summed loss, matched peaks, two-wave reefs, worst trough ratio, the gate, first and second peak heights per reef, and the `Q` that would be selected under the run-11 rule.

### Failures and rollback

**Failures:**
- the `A½ = 0` control does not reproduce stored losses exactly;
- non-finite CPUE.

**Rollback:** none needed, because the frozen observation model is untouched. Adopting `A½ > 0` as the observation model would be a separate decision for the lead, requiring:
- independent evidence for `A½` (for example, paired manta–SCUBA or culling data);
- a re-derived initial state and forcing tables;
- re-validation of earlier conclusions.

### Manta detectability decision (2026-10-13)

[Re-score v2](runs/20261013T_14_detectability_rescore_1991_v2/metadata.toml), [audit](runs/20261013T_14_detectability_rescore_1991_v2/detectability_audit.json), [detection curve](runs/20261013T_14_detectability_rescore_1991_v2/manta_detection_curve.png) and [trajectories](runs/20261013T_14_detectability_rescore_1991_v2/detectability_trajectories.png). New checked-in tools:
- `rescore_detectability.jl` (stored trajectories → both observation treatments);
- `audit_1991_detectability.py`;
- `plot_manta_detectability.py`.

**v1 is superseded.** [v1](runs/20261013T_13_detectability_rescore_1991_v1/metadata.toml) scored the `smoothed_3y` treatment only. It was kept, and v2 (raw and smoothed) was run after the audit finding below.

**Checks.** The `A½ = 0` control reproduces every stored loss exactly (0.0) for 14 arms × 4 reefs × 2 observation treatments. CPUE is finite.

**Audit finding (pre-existing; affects every gate report since the gate was set).**
- **Where the 7 comes from.** The gate requires "all seven observed peaks matched". Those 7 peaks exist on the **raw** LTMP series: Lizard 1998/2013, MacGillivray 1996/2015, North Direction 1995/2013 and Eyrie 2013.
- **What the audits counted.** They count matched peaks from the **`smoothed_3y`** score. Smoothing removes North Direction's 1995 peak (0.56, between 0.06 and 0.12), leaving 6, so "7/7" could never be reached as audited.
- **Effect on run 11.** Against the raw series, the run-11 baseline `wc_b150_q20` already matches **7/7**: since run 11 it has failed the gate **only** on the trough ratio (0.25). Run 15 `acrpoc_g0` matches 6/7 raw.
- **Decision.** The gate definition is **not** changed here: the lead decides. Both counts are reported from now on.

**Result.**

| `wc_b150_q20` under | Summed loss | Matched raw / smoothed | Worst trough ratio (Lizard) | Mean first peak (COTS/tow) | Model mean 2000–10 CPUE (obs 0.08) |
| --- | --- | --- | --- | --- | --- |
| `A½ = 0` (frozen) | 9.06 | 7/7 · 6/6 | 0.252 | 0.52 | 0.118 |
| `A½ = 1` | 9.14 | 7/7 · 6/6 | 0.198 | 0.50 | 0.105 |
| `A½ = 3` | 9.32 | 7/7 · 6/6 | 0.150 | 0.48 | 0.089 |
| `A½ = 10` | 9.81 | 7/7 · 6/6 | 0.102 | 0.40 | 0.062 |
| `A½ = 30` | 13.00 | 6/7 · 5/6 (Lizard second peak lost) | — | 0.29 | 0.037 |

**Gate passes with the raw count.** Five arms pass, all β 1.5 and all at a `Q` above the loss-selected 20:

| Arm | `A½` | Worst trough ratio | Summed loss |
| --- | --- | --- | --- |
| `wc_b150_q30` | 10 | 0.092 | 10.20 |
| `wc_b150_q45` | 10 | 0.096 | 10.34 |
| `wc_b150_q30` | 30 | 0.069 | 10.98 |
| `wc_b150_q45` | 30 | 0.071 | 11.13 |
| `wc_b150_q70` | 30 | 0.082 | 11.62 |

No arm passes with the gate as previously audited (smoothed count against 7).

**Against the predictions:**
- **Troughs.** They fall as predicted: the worst trough ratio goes from 0.25 to 0.10 at `A½ = 10`, and modelled 2000–10 CPUE reaches the observed level (0.06 vs 0.08).
- **Peak heights.** They fall, as predicted (−22% mean first peak at `A½ = 10`), so the **first-wave deficit grows** (Lizard 0.24 vs 1.92). The frozen loss rises monotonically with `A½`, because its height terms dominate and its trough terms are weak.
- **The 1999 drop does *not* become abrupt.** Lizard's 1998→1999 ratio is 0.61 → 0.50 at `A½ = 10`, against 0.03 observed. The model's decline is spread over 1998–2003, so detectability alone does not explain item 4.
- **`Q`.** It stays at 20 under the run-11 rule for `A½ ≤ 3`, and moves to 30 at higher `A½` for one β.

**Interpretation.**
- **What this shows.** Within the range the evidence allows, density-dependent manta detection can reconcile the modelled troughs (2–3 adults/ha) with near-zero manta counts. It is the first route found to the trough criterion.
- **What it does not show.** `A½` is **not identified by independent evidence**. The passing arms need `A½ ≥ 10` and `Q ≥ 30`, so choosing them because they pass would be tuning the observation model to the gate. The first-wave height deficit and the gradual 1999 decline are untouched.
- **Consistency gap.** The simulations still start from a 1991 state inverted with the linear model.

**Decision.**
- Nothing is promoted, and the frozen observation model is unchanged.
- The detectability curve is carried as a labelled sensitivity with `A½ = 10` as the working comparator.
- **For adoption, the lead needs:** independent evidence for `A½` (paired manta–SCUBA or culling data, for example from the Fernandes et al. 1990 or MacNeil et al. 2016a datasets), a decision on the gate's peak count, and a detection-consistent initial state (next screen).


## Preregistered reef-differentiated 1991 starting state (2026-10-13)

Status: **directed by the project lead (2026-10-13), item 2 follow-up (a); preregistered before implementation.** Runner plus one new derived input; no ADRIA, COTSMod, default or objective change.

### Diagnosis motivating the screen

- **Each reef's seed comes from its own survey.** The 1989–91 IDW surface (tows/distance² weights, 1 km floor) is dominated by each calibration reef's own survey:

| Reef | 1989–91 COTS / tows | CPUE |
| --- | --- | --- |
| Lizard | 2 / 255 | 0.0078 |
| MacGillivray | 0 / 15 | 0 |
| North Direction | 1 / 16 | 0.0625 |
| Eyrie | 1 / 133 | 0.0075 |

- **North Direction's seed is a one-animal survey.** Its 3.3 adults/ha (8× Lizard) rests on one COTS in 16 tows. This is sampling noise, and it inverts the first-wave ranking (run 12).
- **The seed was inverted with the linear observation model.** The detectability sensitivity shows that true low-density populations may be several times larger than CPUE/q implies.

### Treatments

**Initial CPUE surface:**
- **`idw`.** The current surface (reference).
- **`eb`.** Each donor reef's 1989–91 rate is replaced by its empirical-Bayes (Poisson–gamma) posterior mean `(c + α) / (t + β)`, and the IDW is then re-run with unchanged weights, radius and minimum donors.
  - The prior is fitted by negative-binomial marginal maximum likelihood to the **31 regional donors** (those within 111 km of any domain site), not to the GBR-wide set, which includes 1989–91 outbreak reefs: α 0.386, β 80.9, prior mean 0.0048 COTS/tow.
  - Resulting donor rates: North Direction 0.0625 → 0.0143, Lizard 0.0078 → 0.0071, Eyrie 0.0075 → 0.0065, MacGillivray 0 → 0.0040.
  - New checked-in builder: `build_1991_eb_surface.py`, with hashes and the prior recorded.

**Inversion from CPUE to adults:**
- **`lin`.** `A = C / q` (reference).
- **`h10`.** The detection-consistent inverse of `C = q·A²/(A + A½)` at `A½ = 10`: `A = (C/q + sqrt((C/q)² + 4·A½·C/q)) / 2`.
  - Example: CPUE 0.0071 → 2.5 adults/ha against 0.47 linear.
  - `A½ = 10` is the detectability screen's working comparator. It is a labelled assumption, not identified.

Juvenile and subadult preloads keep the 6:3:1 ratio to adults.

### Arms (8)

{`idw`, `eb`} × {`lin`, `h10`} × `Q` {20, 30}, on the baseline `wc_b150_q20` model (working coral, β 1.5, `acrpoc` composition, Ext×5, gradual Allee). `Q` 30 is included because the gate-passing detectability arms needed `Q ≥ 30`.

**Reference.** `wc_b150_q20_idw_lin` must replay run 11 `wc_b150_q20` exactly.

**Scoring.** The runner's frozen linear score is kept. All arms are also re-scored with `rescore_detectability.jl` (both observation treatments, all `A½`). The `h10` arms are read primarily under the matching `A½ = 10` observation model; this is the consistent pairing.

### Expected consequences

- **`eb`.** North Direction's seed falls from about 3.3 to under 1 adult/ha, so its first peak falls and comes later (currently 0.62 vs 0.56). Lizard is almost unchanged and MacGillivray rises slightly, so the first-wave ranking becomes less inverted. The summed loss may rise if North Direction loses its first-peak match.
- **`h10`.** Every seed rises 3–5× (Lizard about 2.5, North Direction about 9 adults/ha), crossing the Allee scale (3/ha) sooner, so first waves come earlier and larger.
  - Risk: peaks before the observed 1996–98.
  - Under the matching observation model, the larger true peaks are partly offset by lower detection at the start.
- **Lizard's first-wave deficit.** It is not expected to close, because the seed is not Lizard-specific enough. If confirmed, the deficit lies in supply (local retention or structured external supply), not in the starting state.

### Report

- per reef: 1991 seed adults/ha, first-peak year and height, second wave and trough ratio;
- summed loss under both scorings, matched peaks (raw and smoothed) and the gate;
- the observed and model first-wave rankings.

### Failures and rollback

**Failures:**
- an inexact reference replay;
- an EB surface that does not cover every site, or a non-finite surface;
- non-finite states;
- external supply other than Ext×5.

**Rollback:** the `idw` surface with linear inversion.

### Starting-state decision (2026-10-13)

[Run](runs/20261013T_17_initial_state_rebaseline_1991_v1/metadata.toml), [EB surface](runs/20261013T_15_1991_initial_surfaces_eb/metadata.json), [re-score](runs/20261013T_18_detectability_rescore_initial_state_1991_v1/metadata.toml), [audit](runs/20261013T_17_initial_state_rebaseline_1991_v1/initial_state_audit.json) and [plot](runs/20261013T_17_initial_state_rebaseline_1991_v1/initial_state_rebaseline.png).

New checked-in tools:
- `build_1991_eb_surface.py`. It reproduces the stored IDW surface to 1e-12 before shrinking.
- `audit_1991_initial_state_rebaseline.py`.

**Operational failure.** The first attempt, [run 16](runs/20261013T_16_initial_state_rebaseline_1991_v1/FAILED.txt), stopped before any arm. The runner's EB-surface check demanded exact equality on a float column copied through pandas, which can differ by 1 ulp. The check was changed to `atol = 1e-12` and the run repeated as run 17. An empty re-score directory created by a call against the missing run was removed; it never held output.

**Checks.**
- `wc_b150_q20_idw_lin` replays run 11 `wc_b150_q20` exactly (0.0).
- External supply is identical in all arms.
- The re-score control reproduces all stored losses (0.0).

**1991 seeds (adults/ha):**

| Surface / inversion | Lizard | MacGillivray | North Direction | Eyrie |
| --- | --- | --- | --- | --- |
| `idw` / `lin` | 0.52 | 0.24 | 3.26 | 0.46 |
| `eb` / `lin` | 0.47 | 0.35 | 0.82 | 0.41 |
| `idw` / `h10` | 2.55 | 1.66 | 7.57 | 2.38 |
| `eb` / `h10` | 2.41 | 2.06 | 3.30 | 2.23 |

**Result.**

| Arm (`Q` 20) | Read with | Summed loss | Matched raw / smoothed | Worst trough ratio | Gate (raw count) | First wave 1995–2000: Lizard / MacG / North Dir. / Eyrie | Rank corr. vs observed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `idw_lin` (ref) | `A½ = 0` | 9.06 | 7/7 · 6/6 | 0.252 | no | 0.43 / 0.54 / 0.93 / 0.84 | −0.8 |
| `eb_lin` | `A½ = 0` | 9.10 | 7/7 · 6/6 | 0.253 | no | 0.43 / 0.54 / 0.92 / 0.84 | −0.8 |
| `idw_h10` | `A½ = 10` | 8.56 | 6/7 · 6/6 | 0.087 (two-wave reefs 2) | no | 0.24 / 0.42 / 0.26 / 0.68 | −0.8 |
| **`eb_h10`** | **`A½ = 10`** | 9.92 | 7/7 · 6/6 | **0.088** | **yes** | 0.25 / 0.42 / 0.84 / 0.69 | −0.8 |

At `Q = 30`, `eb_h10` (0.081) and `idw_h10` (0.081) also pass under `A½ = 10`.

**Against the predictions:**
- **`eb` and North Direction.** `eb` cut North Direction's seed 4×, but **its first wave did not change** (0.93 → 0.92). This was not predicted. Lizard and MacGillivray did not change either.
- **The starting state hardly shapes the first-wave size.** First waves are set by supply and local production during 1992–96, not by the 1991 seed. The predicted easing of the inverted ranking **did not occur**: the rank correlation stays −0.8 in every arm.
- **Large seeds hurt.** As the risk stated, `idw_h10` (North Direction 7.6 adults/ha) grazes North Direction's food early: its first wave collapses to 0.26 and comes late (mean first-peak year 2001), losing a matched peak.
- **Lizard's first-wave deficit is not closed, as predicted.** It is 0.25–0.43 against 1.92.

**Not predicted: abrupt declines under the consistent pairing.** With the detection-consistent start read under the matching observation model (`eb_h10` at `A½ = 10`), the first-wave decline becomes abrupt on MacGillivray and North Direction:

| Reef | Linear one-year drop after the 1997 peak | `eb_h10` one-year drop |
| --- | --- | --- |
| MacGillivray | 0.30 | **0.09** |
| North Direction | 0.35 | **0.19** |
| Lizard | 0.70 | 0.46 |

Observed 1998→1999: Lizard 0.03, MacGillivray 0.00, North Direction 0.00. The modelled peaks fall in 1997, a year before the observed crash year.

**First internally consistent gate pass.** `eb_h10` at the loss-selected `Q = 20` is the first internally consistent configuration to pass the qualitative gate (raw peak count):
- the empirical-Bayes 1991 state;
- inverted with the same detection curve used to observe it.

It does so **only** under assumptions the lead has not adopted:
- `A½ = 10`, which is not identified by independent evidence;
- the raw-series reading of "7 observed peaks".

It also leaves the first-wave heights far too low (Lizard 0.25 vs 1.92) and the between-reef ranking inverted. It is a **candidate for evidence gathering, not for promotion**.

**Decision.**
- Nothing is promoted, and the baseline stays `wc_b150_q20` (`idw`, `lin`) under the frozen observation model.
- The starting state is closed as a lever for first-wave height and ranking.
- The `eb` surface is recommended as the better-founded 1991 input (it treats single-animal surveys as noise), pending the lead's decision; it is neutral for the loss (9.06 → 9.10).
- **`eb_h10` with `A½ = 10`** is recorded as the leading *detection-consistent* candidate.

**Next (lead to choose):**
- (i) Evidence for `A½`, and a decision on the gate's peak count. Together these decide whether `eb_h10` becomes the working configuration.
- (ii) First-wave height and ranking: test supply, either local retention or self-recruitment in the Lizard–Cairns initiation region, or spatially structured external supply.
- (iii) Items 3–4 on whichever configuration (i) settles.


## Amended qualitative gate G2′ (2026-10-13)

Status: **directed by the project lead (2026-10-13): "use both [peak counts] but amend the gate so we're actually comparing what we need … maybe the gates are too strict".** Written before any arm was evaluated against it. It is a **post-hoc amendment**: the old gate's results for runs 11–18 had already been seen, and both the old and the new gate are reported for every arm.

### Why the old gate compared the wrong things

1. **Peak count.** "All 7 observed peaks matched" was audited on the smoothed score, which has only 6 observed peaks (smoothing removes North Direction 1995). It was unattainable as audited.
2. **Trough ratio ≤ 0.10** divides the model trough by the **model's own** smaller peak. When peaks are too low (a peak-height failure, which belongs to promotion gate 3), a trough sitting at the observed level still fails. On Lizard, the observed smoothed peaks are 0.96 / 0.74 with an inter-wave mean of 0.040 COTS/tow. The model's peaks are 0.35 / 0.15 with a comparable trough, giving a ratio of 0.25.
3. **No observational test of the trough.** The ratio never compared trough level with what manta tows actually recorded, or with their sampling uncertainty.

### G2′: a mechanism produces the missing qualitative behaviour when all of the following hold

1. **Two waves.** Two detected model peaks (frozen detector) on each of Lizard, MacGillivray and North Direction. *Unchanged.*
2. **Peaks matched, like with like.**
   - All observed peaks are matched under the **raw** observation treatment (7/7);
   - **and** all are matched under the **`smoothed_3y`** treatment (6/6).
3. **Distinct waves.** On each two-wave reef, the smoothed model minimum between its first two peaks is **≤ 0.5 × the smaller of those peaks**. This only excludes plateaus and shoulders, which the 0.12 normalised-prominence rule of the detector can still admit.
4. **Trough level consistent with the observations.**
   - **Trough years.** On each two-wave reef, these are the years strictly between the two observed raw peaks in which observed CPUE < 0.22 COTS/tow (the outbreak level). Surveyed years: Lizard 1999–2011, MacGillivray 1999–2011 (excluding the 2004 survey at 0.38, which is ≥ 0.22), North Direction 1996–2011. There are no 2012 surveys, which is why each window ends in 2011. *(Wording corrected after the evaluation output was read. The rule itself was not changed; only its description was wrong. The earlier text said the 2004 MacGillivray survey was included and that the windows ran to 2012.)*
   - **Test.** The mean model CPUE over those years, under the observation model in use, must be **≤ the upper 95% Poisson bound** of the observed mean (Σ counts / Σ tows, Garwood interval).
   - **Strictness.** Clustered COTS counts are over-dispersed, so the Poisson bound is narrower than the true sampling interval. The test therefore errs towards failing.

The old trough ratio is still reported as a diagnostic. Peak heights, peak prominence and timing errors stay with **promotion gate 3** and are reported alongside every G2′ result.

### Evaluation

A checked-in script, `evaluate_gate_g2prime.py`, evaluates G2′ on the re-score outputs (raw and smoothed scores, plus detected-CPUE series at every `A½`). It covers runs 11 and 15 (re-score v2) and run 17 (re-score run 18). For each arm it reports the old gate, G2′ and each G2′ component.

**Observation-model pairing.**
- Under the frozen linear model, G2′ is read at `A½ = 0`.
- The detection-consistent arms (`h10`) are read at `A½ = 10`.
- A linear-seeded arm read at `A½ > 0` is a labelled sensitivity, not a consistent configuration.

### G2′ evaluation and evidence for `A½` (2026-10-13)

**G2′ evaluation** ([evaluation](runs/20261013T_19_gate_g2prime_evaluation_v1/gate_g2prime.json), `evaluate_gate_g2prime.py`, run on re-scores 14 and 18).

**Observed inter-wave windows** (mean COTS/tow; upper 95% Poisson bound):

| Reef | Years | Mean | Upper bound | Sample |
| --- | --- | --- | --- | --- |
| Lizard | 1999–2011 | 0.041 | 0.054 | 41 COTS in 1,029 tows |
| MacGillivray | 1999–2011 | 0.021 | 0.054 | 2004 excluded |
| North Direction | 1996–2011 | 0.054 | 0.089 | |

| Arm | Read at `A½` | Peaks matched raw · smoothed | Separation (worst) | Model trough CPUE: Lizard / MacG / North Dir. | G2′ | Old gate (raw count) |
| --- | --- | --- | --- | --- | --- | --- |
| run 15 `acrpoc_g0` | 0 | 6/7 · 5/6 | — (North Dir. one wave) | 0.133 / 0.063 / 0.312 | no | no |
| baseline `wc_b150_q20` | 0 | 7/7 · 6/6 | 0.25 | 0.105 / 0.092 / 0.206 | no (trough) | no |
| baseline `wc_b150_q20` | 10 | 7/7 · 6/6 | 0.10 | 0.048 / 0.040 / 0.136 | no (North Dir. trough) | no |
| `wc_b150_q20_eb_h10` | 10 | 7/7 · 6/6 | 0.088 | 0.041 / 0.040 / 0.134 | no (North Dir. trough) | yes |
| `wc_b150_q30_idw_h10` | 10 | 7/7 · 6/6 | 0.081 | 0.045 / 0.038 / 0.065 | **yes** | yes |

**How G2′ changes the picture.**
- **Strictness removed.** The unattainable peak count and the peak-height confounding in the trough ratio no longer block a pass. Every working-coral arm now meets components 1–3.
- **Strictness added.** The absolute trough test fails on every arm read with the linear observation model, on all three reefs.
- **Under `A½ = 10`.** Lizard and MacGillivray pass the trough test, but **North Direction's modelled inter-wave population is 2.5× the observed level** (0.13 vs 0.054 COTS/tow).
- **The only G2′ pass.** `wc_b150_q30_idw_h10` passes, because its large North Direction seed depletes North Direction's food early. The same arm loses peak timing elsewhere: mean first-peak year 1997.5.

G2′ therefore turns the remaining qualitative failure into one specific, observable target: **the inter-wave population on North Direction** (item 4: the juvenile pipeline keeps refilling adults).

**Evidence for `A½`** (lead's request; searched 2026-10-13):
- **Fernandes et al. (1990), *Coral Reefs* 9.**
  - **Sightability.** Manta tows saw 22.7% of the COTS counted on SCUBA on the same 10 m strips, falling with strip width (under 5% for routine tows).
  - **Drivers.** Sightability was driven chiefly by the **proportion of cryptic starfish** and **reef complexity**.
  - **Correction model.** A correction equation explained only 39% of the variation.
  - **Density.** No density term is reported in the abstract. The companion paper (Fernandes 1990, *Coral Reefs* 9:161–165, on the distribution and density of targets) and De'ath (1992, *Coral Reefs* 11:119–121, comment) are paywalled and were not read.
- **MacNeil et al. (2016), *PeerJ* 4:e2310.**
  - **Diver detection** was 0.82 [0.77, 0.87]. The strongest covariate was **disc width** (above 0.8 for COTS over 30 cm; below 0.5 under 15 cm). There is no density term, and no manta-tow detectability.
  - **Cull CPUE** is **hyperstable** (CPUE ∝ N^0.33). That is the opposite curvature to `A½ > 0`, but it concerns targeted culling dives, not manta tows.
- **Pratchett & Caballes (2025), *Biology* 14:1391.**
  - **Detectability.** About 0.74–0.78 by depletion sampling at Lizard Island (40 transects of 200 m², 0.2–7.2 COTS per 200 m², i.e. 10–360 per ha). Mark–recapture at Rib Reef gave 0.65 by day and 0.78 by night.
  - **Covariates.** **Density did not explain variation among transects.** Detectability rose with size (0.50 at ≤ 20 cm; 0.89–0.91 above 30 cm) and was lower on the slope than the crest (0.70 vs 0.90).
  - **Relevance.** These densities lie mostly above the range where `A½ ≈ 3–10` matters.
- **Pratchett et al. (2026), *Coral Reefs*, doi 10.1007/s00338-026-02841-y (SALAD; Lizard, Cairns and Far Northern regions, 2021–25).**
  - **Detectability.** Mean 0.51, varying by region with **COTS size and feeding**. Feeding starfish were 83% exposed from above, non-feeding 24%.
  - **Implication for manta.** The authors stress that methods relying on detection from above (manta, towed camera) are significantly limited.
  - **Regional means against inferred density:**

| Region | Inferred density (COTS/ha) | Detectability |
| --- | --- | --- |
| Cairns | 4.1 | 0.30 |
| Princess Charlotte Bay | 8.7 | 0.43 |
| Cape Grenville | 17.1 | 0.59 |
| Lizard | 20.3 | 0.55 |

  - **Fitting our curve.** Fitting `d = d_max·A/(A + A½)` to these four points gives **`A½` = 6.0 ± 2.1 COTS/ha** (`d_max` 0.74).

**Assessment.**
- **No study measures manta-tow detectability as a function of density.**
- The direct within-site test (Pratchett & Caballes 2025) found no density effect at 10–360 COTS/ha.
- Across studies, detectability is governed by **size, crypsis and feeding (exposure)**, all of which shift towards cryptic, small and non-feeding animals at low density. This is the mechanism behind `A½`.
- The only quantitative handle is the regional SALAD fit (`A½` ≈ 6 ± 2). It is indirect: four regional means, confounded with size and region, from diver surveys, not manta tows.
- **Conclusion.** The evidence supports `A½` **in the range ~3–10 COTS/ha as plausible**, with the mechanism being size- and exposure-dependent detectability. It does not identify `A½`.
- **Recommended route to identification:** paired manta–SALAD surveys on the same reefs and years (the Lizard region has both, 2021–25), or the paywalled Fernandes (1990) density analysis.


## Preregistered Lizard Island local-retention screen (2026-10-13)

Status: **directed by the project lead (2026-10-13); preregistered before implementation.** Runner-level input treatment plus a checked-in, unit-tested transform; no ADRIA, COTSMod, default or objective change.

**Decisions recorded with this screen:**
- **Empirical-Bayes 1991 surface adopted** as the working starting state (`eb`, linear inversion; run 17 `wc_b150_q20_eb_lin`, loss 9.10). The lead's reason: it is better founded, and it is not the cause of Lizard's suppressed peaks.
- **Gate:** G2′ as amended above, with both peak counts.

### Hypothesis (project lead)

Lizard Island holds a **remnant population** between waves; the observed trough CPUE is about 0.04, non-zero. It also has an **internal lagoon** that may retain larvae. The Owen reef-scale connectivity under-represents this retention, so Lizard's own production does not feed back into its population. That would keep both of Lizard's waves small, and the updated (finer) connectivity should eventually resolve it.

### Evidence

- **Owen matrices (2018–22).** Lizard is three domain reefs: 14-116b (North East), 14-116c (Coconut Bay) and 14-116d (**Lagoon**).
  - **Self-retention.** The diagonal entries are 0.0001–0.0008 for 14-116b and 14-116c, and 0.0007–0.0049 for the Lagoon.
  - **Within-complex supply.** All three sub-reefs send most of their within-complex supply to the Lagoon.
  - **Survival.** Each row sums to about 0.06–0.08, so about 93% of production dies in the plankton.
  - **Retained share.** The Lizard complex keeps 1–8% of its surviving larvae.
- **Flushing.** The Lizard lagoon flushes in hours to 10 days, longer in low wind (Philipps & Bellwood 2024, *Coral Reefs*, doi 10.1007/s00338-024-02511-x). COTS larval duration is weeks, so lagoon retention alone is partial.
- **Sticky water.** Dense reef clusters are poorly flushed and retain more larvae ("sticky water"; Andutta et al. 2012). Reef-scale residence times are 0.5–5.6 days.
- **Self-seeding rates.** Large-scale GBR models give self-seeding under 2% for passive coral larvae and over 20% for fast-developing fish larvae.
- **Fine-scale model.** A 30 m Lizard dispersal model exists (Gouezo et al. 2025, *Ecological Applications*).
- **What is missing.** No study gives a COTS retention fraction for Lizard. The levels below therefore bracket "how much retention would be needed"; they are not estimates.

### Mechanism (source-conserving redistribution)

For each annual Owen matrix `M` (rows are sources, columns are sinks; entries are the fraction of source production arriving), and each source reef `s` in a retention group `G`:
- `S = Σ_j M[s, j]` is the source's surviving fraction, and `c = Σ_{j∈G} M[s, j]` is the part already retained in `G`.
- **Within the group:** `M'[s, j] = M[s, j] · ρS / c` for `j ∈ G`. This keeps Owen's within-group pattern, in which the Lagoon dominates.
- **Exports:** `M'[s, j] = M[s, j] · (1 − ρ)S / (S − c)` for `j ∉ G`.
- **Conservation.** Row sums, and so each source's total larval survival, are **preserved exactly**: more larvae stay home, but no extra larvae are created. The maximum possible retention is therefore `S` (about 7% of production).
- **Implementation.** `connectivity_retention.jl` holds the transform, and `test_connectivity_retention.jl` tests row-sum preservation, orientation, the target share and the identity at the Owen share.
- **Unchanged inputs:** external supply and source larval survival.

### Arms (5)

All arms run on the adopted baseline: working coral, β 1.5, `Q = 20`, `eb` seeds with linear inversion, Ext×5 and gradual Allee.

| Arm | Treatment |
| --- | --- |
| `…_owen` | Owen matrices unchanged. **Must replay run 17 `wc_b150_q20_eb_lin` exactly.** |
| `…_liz25` / `…_liz50` / `…_liz90` | Lizard complex `G` = {14-116b, 14-116c, 14-116d} retains ρ = 0.25 / 0.5 / 0.9 of its surviving larvae (Owen: 0.01–0.08). |
| `…_all50` | Discriminating control. Every domain reef retains ρ = 0.5 of its survivors on itself (Lizard as one complex). This tests whether any benefit is Lizard-specific or a generic retention effect that would also inflate North Direction and Eyrie. |

### Expected consequences

**Lizard's own production.** At about 250 larvae/ha/yr at the peak, it supplies 0.06–0.25 × 250 × 0.07 ≈ 1–4 recruits/ha/yr even at ρ = 0.9. That is small against Ext×5 external supply (65–100 recruits/ha/yr). The feedback loop (more adults → more production → more retention) is weak, because only about 10% of settlers reach adulthood.

**Prediction:**
- **First wave.** A **modest** rise for Lizard (≤ +30% at ρ = 0.9) with no timing change.
- **Remnant and second wave.** A larger inter-wave remnant population, and an earlier or larger second wave, because the trough population is then partly self-sustaining.
- **Risk.** Higher troughs, which can fail the G2′ trough test.
- **Control.** `all50` should raise every reef's troughs.

**Discrimination.** If even ρ = 0.9 leaves Lizard's first wave far below observed, then retention within source-conserving bounds cannot explain it. Lizard's deficit would then need larval supply beyond the current external boundary (an initiation-region source), or higher larval survival near Lizard, not retention.

### Report

- **Gates.** G2′ (linear observation), with a sensitivity re-score at `A½ = 10`, and the old gate.
- **Per reef:** first- and second-peak heights and years, inter-wave mean CPUE, and the trough population.
- **Lizard fluxes:** local retention, internal immigration and external immigration, 1992–2015.
- **Mass check:** maximum absolute row-sum change over all annual matrices.

### Failures and rollback

**Failures:**
- an inexact reference replay;
- a row-sum change above 1e-12;
- a negative matrix entry;
- a retention group not found in the matrix;
- non-finite states.

**Rollback:** Owen matrices unchanged.

### Local-retention decision (2026-10-13)

[Run](runs/20261013T_20_local_retention_1991_v1/metadata.toml), [audit](runs/20261013T_20_local_retention_1991_v1/local_retention_audit.json), [plot](runs/20261013T_20_local_retention_1991_v1/local_retention.png), [re-score](runs/20261013T_21_detectability_rescore_local_retention_1991_v1/metadata.toml) and [G2′](runs/20261013T_22_gate_g2prime_local_retention_v1/gate_g2prime.json).

New checked-in tools: `connectivity_retention.jl`, `test_connectivity_retention.jl` (14/14 pass), and `audit_1991_local_retention.py`. The parse test now also covers `rescore_detectability.jl` and the two retention files (35/35 pass).

**Checks.**
- `…_owen` replays run 17 `wc_b150_q20_eb_lin` exactly (0.0).
- The largest row-sum change over all annual matrices is ≤ 1e-12, recorded in metadata.
- External immigration is identical in every arm.
- The re-score control is exact.

**Result.** Retention **does not change Lizard materially**, and no arm passes G2′.

| Arm | Summed loss | Lizard first / second wave (obs 1.92 / 0.74) | Lizard inter-wave CPUE | North Dir. inter-wave CPUE (obs 0.054) | Lizard settlers from domain sources 1992–2015 |
| --- | --- | --- | --- | --- | --- |
| `owen` (ref) | 9.10 | 0.432 / 0.160 | 0.105 | 0.128 | 15% (local 4, internal 285, external 1,670 /ha) |
| `liz25` | 9.10 | 0.435 / 0.158 | 0.105 | 0.127 | 16% |
| `liz50` | 9.10 | 0.439 / 0.158 | 0.105 | 0.127 | 18% |
| `liz90` | 9.09 | 0.446 / 0.158 | 0.105 | 0.126 | 20% (local 60) |
| `all50` | 9.01 | 0.428 / 0.179 | 0.104 | 0.136 | 12% |

**Against the predictions.** These were confirmed, and the effect was even smaller than predicted:
- **Lizard.** The first wave rises by +3% at ρ = 0.9 against a predicted ceiling of +30%. The second wave and the inter-wave population are unchanged.
- **Why so small.** Source-conserving retention can return at most the surviving ~7% of production, and Lizard's own production over 1992–2015 (1,217 /ha) is smaller than its external supply (1,670 /ha).
- **Control.** `all50` hardly changes anything either (North Direction's trough up slightly).

**Diagnosis: Lizard's first wave is limited by food, not supply.**
- **Supply is high.** Lizard receives 65–140 settlers/ha/yr in 1994–99.
- **Maturation is low.** Only 3–15 subadults/ha/yr mature, because Lizard's modelled fast coral falls from 0.07–0.09 (1992–93) to 0.04 after the 1994 heat stress and to 0.02 by 1997. Food survival then falls to 0.48, and body condition to 0.5, by 1998.
- **North Direction contrast.** It keeps fast cover at 0.05–0.11 until 1996 and food survival at 1.0 through 1997, so the same settlement produces 61 adults/ha against Lizard's 29.
- **The data puzzle.** Observed Lizard reached 1.92 COTS/tow (about 128 adults/ha) on photo-transect fast cover of about 0.02–0.04. This is the food-scale question first raised after run 14.
- **Interpretation.** The threshold θ = 0.05 and the condition saturation of 0.10 are defined on 9 m *Acropora* + Pocilloporidae cover. They do not let Lizard's starfish feed at the scale observed. This is consistent with the run-12 finding that a higher starting composition helps Lizard only weakly, because the 1994 heat stress resets the fast cover.

**Decision.**
- Nothing is promoted, and the Owen matrices are unchanged.
- Within source-conserving bounds, local retention is **closed** as an explanation for Lizard's suppressed peaks.
- `connectivity_retention.jl` is kept for the updated connectivity, which should resolve sub-reef retention directly.
- **Two ways retention could still matter:**
  - higher survival for retained larvae (a non-conserving treatment, not tested);
  - a Lizard-specific supply that the current Ext×5 boundary misses.
  Neither changes the food limit identified here.

**Next (lead to choose):**
- (a) **Food scale for Lizard.** Gather evidence on what Lizard's COTS ate in 1995–98: the full prey breadth (*Montipora*, *Isopora*, crest *Acropora* above 9 m) and the reef-wide fast cover. Then decide whether θ and the condition saturation should be defined on a broader prey pool. This is a mechanism change on an independently set threshold, so it needs evidence, not a fit.
- (b) **North Direction's inter-wave population** (the remaining G2′ failure; item 4): why 80–130 subadults/ha persist through 1999–2002 and refill adults.
- (c) **Bleaching β.** β = 1.0 instead of 1.5 removes less of Lizard's pre-outbreak food in 1994 (run 11: Lizard first peak 0.41 vs 0.35). A small lever, already measured.
