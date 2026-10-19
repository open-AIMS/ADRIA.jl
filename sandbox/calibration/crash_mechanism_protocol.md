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


## Lead's decisions and new evidence (2026-10-14)

**North Direction (G2′ component 4).** The lead is not concerned by North Direction's higher inter-wave density. Unlike most reefs, it continued to record > 0.05 COTS/tow during the 2000s, but it was surveyed only intermittently, so the remnant population was probably unobserved in some years. **The G2′ trough test therefore gates Lizard and MacGillivray only.** North Direction's inter-wave level is still reported. `evaluate_gate_g2prime.py` gains a `--trough-reefs` option; the default keeps all three reefs, for replay.

**Pratchett (2005), *Coral Reefs* 24:453–462** (the lead supplied it): the Lizard outbreak of 1994–99, from 10 locations × 2 depths, 200 m² belt transects searched twice (under 0.8% missed).
- **Densities.** About 12.5 COTS/ha in October 1994 (from Sweatman et al. 1998). The build-up was gradual, reaching 40 COTS/ha in October 1996, then 40–55 COTS/ha until June 1998.
- **Crash.** 0.92 → 0.06 per 200 m² between June 1998 and January 1999 (46 → 3 per ha), then under 0.5 per ha by September 1999.
- **Cohorts.** At least four, with recruitment from about 1990 onwards: a gradual build-up, a probable **primary outbreak** from which outbreaks spread south, and **self-seeding** suggested. Small starfish were found only at sheltered back-reef sites (North Reef, Casuarina, Corner Beach). Starfish were **almost absent from the lagoon**.
- **Cause of the crash.** Coral cover was **still 22.9% at the end, against 32.2% at the start**, so the crash was **not food-driven**. All size classes declined, so it was **not senescence**. An infectious disease was demonstrated in starfish from Lizard in January 1999 (Pratchett 1999).
- **Implications for this work:**
  - **Food.** Lizard's real outbreak was not food-limited (reef-wide cover about 23–32%), whereas the model's Lizard fast cover is 0.02–0.04.
  - **Crash mechanism.** The crash is a disease-like, density-associated mortality, not starvation.
  - **Retention.** The lagoon is not where Lizard's starfish were.

**Manta detectability from paired data**
- **Run:** [run 02 (v2)](runs/20261014T_02_manta_detectability_evidence_v2/manta_detectability_fit.json), via the checked-in `analyse_manta_detectability.py`; derived tables are in the git-ignored `sandbox/data/derived/manta_detectability_pairs_lizard/`.
- **Superseded:** [v1](runs/20261014T_01_manta_detectability_evidence_v1/SUPERSEDED.txt), whose pairing reused manta surveys.
- **Data and pairing.** SALAD 2021–26 (the lead's master workbook; "Eagle" is Eyrie) and the Pratchett (2005) belt densities are paired with LTMP and Control Program manta surveys on the same reef within 6 months: one pair per manta survey, 47 pairs.
- **Effective `q` (manta COTS/tow per diver COTS/ha):**
  - **2021–26:** 0.002–0.006 at diver densities of about 2–35 COTS/ha (SALAD hotspots up to 80/ha). That is **3–7× below the frozen `q = 0.015`**, with no detectable trend within that range.
  - **Lizard 1996–97 outbreak:** 0.022–0.045 (40–43 COTS/ha against 0.90 and 1.92 COTS/tow), **above 0.015**.
  - **During the crash (November 1998):** 13.5 COTS/ha on belt transects against 0.05 COTS/tow on manta in October 1998, an effective `q` of 0.004.
- **Model fit.** Across all pairs the power model gives an exponent of 1.4–1.7 (hyperdepletion), and the linear model is no better. `A½` is **not identified**: the profile 95% interval is unbounded, because the two density clusters come from different eras, methods and spatial scales (SALAD follows feeding scars into hotspots; belt transects and manta tows are reef-wide).
- **Conclusion.** The evidence contradicts a constant `q` at both ends, and supports strong density-dependent manta detection. It does not fix `A½`.
- **Most material finding.** At Lizard's real outbreak density (diver 40–55/ha), manta read 0.9–1.9 COTS/tow. The frozen `q = 0.015` maps 1.92 COTS/tow to 128 adults/ha, about 2.5–3× the diver density. The model's Lizard first-wave deficit against **diver** density is about 1.6× (29 vs 46 adults/ha), not the 5× implied by `q = 0.015`.
- **Who decides.** `q` is a lead-adopted convention; whether to revisit it is the lead's decision.

**Time-convention finding (for the lead; nothing changed).**
- **What `year` means.** In `reef_cots.csv` (and `reef_manta.csv`, `report_year`), `year` is the **AIMS reporting year** (July–June, named for its ending year): every survey from July to December carries the next calendar year. All four calibration reefs are surveyed in September–November, so each observation is labelled **about 1 year after it was made**. For example, Lizard's "1998" peak (1.92) was surveyed in October 1997, and the "1999" crash value (0.05) in October 1998.
- **What the model compares.** Model year Y includes year-Y disturbances (DHW is indexed by the bleaching summer), and the scorer compares it with observations labelled Y, so the model is effectively scored against surveys made about a year earlier.
- **What follows:**
  - modelled peaks are probably not "a year early" as previously read;
  - the 2016–17 coral trough is offset by a year;
  - the observed-COTS forcing table (built by `year`) grazes a year late.
- **Proposed migration** (lead to decide; it changes a time convention, so it needs a compatibility test): score against `floor(survey_date)` (the survey calendar year) for COTS and coral, and rebuild the forcing table on the same basis.

## Preregistered bleaching threshold floor × starting fast-coral share (2026-10-14)

Status: **directed by the project lead (2026-10-14), item (a); preregistered before implementation.**

### Lead's reasoning

- **Bleaching.** The 1994 loss of fast coral to bleaching is suspect: there should be no mortality until heat stress exceeds about 4 DHW, and then only some. A strong response at high levels is acceptable, provided it does not assume mortality at low heat stress.
- **Starting fast cover.** The photo transects sample a small part of each reef, while the COTS densities come from manta tows over the whole reef. The starting fast-coral cover should therefore be allowed to increase where that helps match the data.

### Evidence

- **Effective DHW on the calibration reefs** (5 m, ADRIA's depth attenuation):

| Year | Effective DHW |
| --- | --- |
| 1994 | 3.3–3.6 |
| 1996 | about 1.6 |
| 1998 | about 1.1 |
| 2002, 2004, 2010 | 1.5–2.5 |
| 2016–17 | 7.5–8.8 |
| 2020 | 4.7–6.1 |
| 2022 | 3.5–4.0 |
| 2024 | about 5.7 |

- **Current ADRIA behaviour.** Mortality is the fraction of each coral group's tolerance distribution (truncated normal) lying below effective DHW. The truncation's lower bound is `HEAT_LB` = 1.0 DHW, so every year above 1 DHW kills the low-tolerance tail. ADRIA's own comment in `bleaching_mortality!` reads "with 4.0 DHW-weeks as the minimum susceptibility threshold", which disagrees with the constant.
- **Literature.** GBR coral loss begins at about 3–4 DHW and is severe above about 6 (Hughes et al. 2018, *Nature* 556:492, which ADRIA already cites).
- **Pratchett 2005.** Reef-wide coral cover at Lizard was 32% in October 1996.

### Mechanism

**Opt-in ADRIA switch `ADRIA_BLEACHING_DHW_FLOOR`** (DHW-weeks; default `HEAT_LB` = 1.0, giving identical results).
- `bleaching_mortality!` takes the floor as its minimum susceptibility threshold. It is used for the active-location test, as the lower truncation bound of the tolerance distribution (still `max(floor, previous bleaching DHW)`), and in the post-mortality mean update.
- **Below the floor** there is no mortality.
- **Above the floor** the tolerance mass that lay below it is renormalised above it. *(Corrected 2026-10-14, before any run: the original wording said this makes the response at high DHW stronger. The unit test showed the opposite: renormalising upwards slightly **lowers** mortality above the floor, to 0.992 against 0.996 at 8 DHW for mean 4 ± 1.5. The lead's "strong response at high levels" therefore comes from β, which the factorial crosses with the floor.)*
- `ADRIA_DEBUG_BLEACHING_SCALAR` (β) still applies.
- **Test:** `test_bleaching_floor.jl` checks that the default equals the explicit `HEAT_LB`, that there is no mortality below the floor, and that mortality at 8 DHW is at least the default value.

**Starting fast share.** The observed 1993–95 photo-transect share (`acrpoc`) is multiplied by `m` ∈ {1, 1.5, 2, 3}, capped at 0.9, on every reef, with total cover and within-block structure preserved (run-12 method). The Lizard shares are 0.21 / 0.31 / 0.42 / 0.63.

### Arms (12)

Floor {1 (default), 4} × β {1.0, 1.5} × `m` {1, 2, 3}, on the adopted baseline (working coral, `Q = 20`, `eb` seeds with linear inversion, Ext×5, gradual Allee). The `m = 1.5` level is dropped to keep the grid small.

**Reference.** `floor100_b150_m100` must replay run 17 `wc_b150_q20_eb_lin` exactly.

### Expected consequences

- **Floor 4.** No bleaching loss in 1994, 1996, 1998, 2002, 2004 or 2010. Lizard's pre-outbreak fast cover stays at about 0.07–0.09 rather than falling to 0.04, so the first wave rises and becomes less food-limited (Lizard first peak up perhaps 20–40%). Second waves, which depend on inter-wave fast cover, rise. The 2016–17 and 2020–24 losses are similar or stronger.
- **`m` > 1.** First waves rise on every reef. The risk (run 12) is North Direction and Eyrie overshooting, and Lizard and MacGillivray troughs rising.
- **β.** Matters less under floor 4, because only years above 4 DHW are affected.

### Selection rule (frozen)

The `m` scale is an uncertain **input** that the lead allows to be fitted.
- **Eligible arms** must have a G2′ trough test on Lizard and MacGillivray no worse than the reference, and a reef-wide manta total-cover RMSE (the reef-wide coral observation) no more than 10% worse than the reference. The photo-transect fast RMSE is reported but does not gate, because the lead's premise is that the transects under-represent reef-wide fast cover.
- **Choice among eligible arms:** the minimum summed frozen loss.
- **Floor and β.** These are coral parameters. The floor is set by evidence (4) if the coral `J` (manta component) is not worse; β is then chosen by coral `J`, as before.
- **Readings.** All arms are re-scored under `A½` (sensitivity) and evaluated with G2′ (trough gate on Lizard and MacGillivray).

### Report

- summed loss, G2′, and peaks per reef;
- Lizard's pre-outbreak fast cover (1992–96) and maturation;
- coral `J` components;
- bleaching loss by year.

### Failures and rollback

**Failures:**
- an inexact replay;
- the floor test failing;
- non-finite states.

**Rollback:** `ADRIA_BLEACHING_DHW_FLOOR` unset, and `m = 1`.

### Bleaching-floor decision (2026-10-14)

[Run](runs/20261014T_03_bleaching_floor_1991_v1/metadata.toml), [audit](runs/20261014T_03_bleaching_floor_1991_v1/bleaching_floor_audit.json), [plot](runs/20261014T_03_bleaching_floor_1991_v1/bleaching_floor.png), [re-score](runs/20261014T_05_detectability_rescore_bleaching_floor_1991_v1/metadata.toml) and [G2′, trough gated on Lizard and MacGillivray](runs/20261014T_06_gate_g2prime_bleaching_floor_v1/gate_g2prime.json).

**ADRIA change (opt-in).**
- `bleaching_mortality!(…; heat_lb = HEAT_LB)` takes the floor, and `run_model` reads `ADRIA_BLEACHING_DHW_FLOOR`; it must be finite and ≥ `HEAT_LB`, and unset means `HEAT_LB`.
- `test_bleaching_floor.jl` passes 8/8: the default is identical to the explicit `HEAT_LB`, there is no mortality at or below 4 DHW, and mortality above the floor is ≥ 95% of the default.
- The G2′ evaluator gains `--trough-reefs`. Re-evaluating runs 11, 15, 17 and 20 with the trough gated on Lizard and MacGillivray ([evaluation](runs/20261014T_04_gate_g2prime_lizmacg_trough_v1/gate_g2prime.json)) shows the **adopted baseline `wc_b150_q20_eb_lin` passes G2′ at `A½ = 10`**, but not under the linear observation model.

**Checks.**
- `f100_b150_m100` replays run 17 `wc_b150_q20_eb_lin` exactly (0.0).
- 1991 totals are identical in all arms, and the re-score control is exact.

**Result.** No arm passes G2′ under the linear observation model. Under `A½ = 10` only the reference passes.

| Arm (β 1.5) | Summed loss | Lizard first / second (obs 1.92 / 0.74) | Lizard peak adults /ha (diver obs 40–55) | Lizard fast 1993 → 1994 | MacG / North Dir. first (obs 0.86 / 0.56) | Lizard trough CPUE, linear / `A½` 10 (obs ≤ 0.054) | Coral `J` (manta RMSE) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| floor 1, m 1 (ref) | 9.10 | 0.43 / 0.16 | 29 | 0.091 → 0.042 | 0.54 / 0.92 | 0.105 / 0.048 | 0.363 (0.131) |
| **floor 4, m 1** | 9.79 | **0.59** / 0.17 | **39** | **0.091 → 0.099** | 0.67 / 1.42 | 0.144 / 0.081 | 0.392 (0.144) |
| floor 1, m 2 | 12.49 | 1.61 / 0.50 | 107 | 0.170 → 0.083 | 0.71 / 3.35 | 0.516 / 0.417 | 0.286 (0.121) |
| floor 4, m 2 | 11.12 | 1.85 / 0.90 | 123 | 0.170 → 0.179 | 1.06 / 3.54 | 0.762 / 0.653 | 0.332 (0.137) |
| floor 4, m 3 | 17.19 | 2.22 / 0.87 | 148 | 0.247 → 0.257 | 2.21 / 3.83 | 1.05 / 0.94 | 0.325 (0.141) |

At floor 4, β 1.0 and 1.5 give identical COTS dynamics, because no bleaching occurs before 2016.

**Against the predictions:**
- **Floor 4 and Lizard: as predicted.** The 1994 loss is removed, and Lizard's first wave rises by 36%; its peak adult density of 39/ha now matches Pratchett's diver density of 40–55/ha. Second waves barely change.
- **Floor 4 overall: not predicted.**
  - North Direction's first wave overshoots (1.42).
  - Troughs rise, and Lizard fails the trough test even at `A½ = 10`.
  - Coral `J` worsens, mainly through manta total cover, because coral now grows more.
  - **Reason:** the floor also removes the 1.5–2.5 DHW bleaching losses of 1996, 1998, 2002, 2004 and 2010. Those losses had been supplying part of the chronic coral loss that the coral fit needed (the "Mellin + 3%" term was fitted with them in place).
- **`m` > 1.** First waves rise everywhere. Lizard reaches its *manta-scale* target only at m 2–3 (adults 107–148/ha, which is what `q = 0.015` implies), and Lizard's second wave reaches 0.9 at floor 4, m 2. But North Direction and Eyrie overshoot heavily, and troughs explode. A uniform multiplier cannot help Lizard without breaking the other reefs.

**Frozen selection.** Only the reference meets the trough and manta-cover conditions, so `m = 1`. The floor-4 condition (manta coral `J` not worse) fails. **The baseline is unchanged.**

**Synthesis for the lead.**
- (1) **Lizard's first-wave height is largely a units question.** On the diver scale (Pratchett 2005), Lizard's outbreak was 40–55 COTS/ha. With floor 4, the model reaches 39 adults/ha. The apparent 5× deficit came from converting 1.92 COTS/tow with `q = 0.015` (128/ha). The paired data give `q` ≈ 0.02–0.045 at outbreak densities.
- (2) **The floor is evidence-supported, but it exposes a coral-side gap.** With mortality below 4 DHW removed, the coral model needs its chronic loss re-fitted. That is a coral-side step: a forced-COTS fit of the background μ under floor 4 (and possibly the time-convention correction) before COTS screens use floor 4.
- (3) **A uniform starting-share multiplier is rejected.** A Lizard-specific (reef-wide) share would need its own evidence, such as Pratchett's unpublished reef-wide coral data (32% total cover, 1996) or the RHIS benthic records in the lead's CSIRO workbook. It should not be fitted.


## Lead's decisions (2026-10-14, second round) and Pratchett-anchored observation models

**Time convention.** No migration: observations stay scored at their LTMP reporting year. Survey timing is now **tracked**.
- **Tracker.** The new `build_survey_timing.py` writes [run 07](runs/20261014T_07_survey_timing_v1/survey_timing_summary.json), with tables in the git-ignored `sandbox/data/derived/ltmp_survey_timing/`.
- **Lag.** 99% of calibration-reef surveys (COTS manta, coral manta, photo transects) are made in September–December, a median 0.21 years before the end of the model year they are scored with. They therefore predate that model year's bleaching summer.
- **Surveys predating ≥ 3 DHW bleaching the model has already applied:**
  - **1994 (3.3–3.5 DHW):** surveyed October 1993.
  - **2017 (about 8 DHW):** surveyed September–December 2016, so they show the 2016 bleaching but not 2017's.
  - **2020** (MacGillivray, 5.5 DHW).
  - **2021–22** (3.5–4.0 DHW).
  - **2024** (5.5–5.8 DHW): surveyed November 2023, and January 2024 on Eyrie.
- **Use.** Audits should read these flags whenever a coral or COTS residual falls in a flagged year.

**Pratchett-anchored observation models (re-score only; seeds still inverted at `q = 0.015`).**
- **Fits to Pratchett's four Lizard pairs:**
  - Poisson-fit linear `q` = **0.030**;
  - saturating fit runs to its `A½` bound (`q_max` 0.055 at `A½` 30), so it is not identified;
  - power model exponent b = 3.5. The November 1998 crash point drives the steepness.
- **Re-score.** `rescore_detectability.jl` gains `--q=`. [`q` 0.030](runs/20261014T_08_rescore_pratchett_q0030_v1/metadata.toml) and [`q` 0.055](runs/20261014T_09_rescore_pratchett_q00552_v1/metadata.toml) were evaluated with [G2′ (0.030)](runs/20261014T_10_gate_g2prime_rescore_pratchett_q0030_v1/gate_g2prime.json) and [G2′ (0.055)](runs/20261014T_11_gate_g2prime_rescore_pratchett_q00552_v1/gate_g2prime.json). Both directories were created without a sequence number and renamed (`RUN_RENAMES.csv`).
- **Linear `q` = 0.030** lowers the baseline's frozen loss **9.10 → 7.29**, the best to date. Peaks double (Lizard 0.69, MacGillivray 0.70, North Direction 1.42 COTS/tow smoothed). But trough CPUE doubles too (Lizard 0.21, MacGillivray 0.19, against an observed upper bound of 0.054), so G2′ fails.
- **Saturating (0.055, `A½` 30)** loses Lizard's second peak (loss 10.0).
- **Conclusion.** A linear `q` cannot match both the outbreak peaks and the troughs. Pratchett's pairs need **strong curvature**, high detection at outbreak density and very low detection in troughs, which is consistent with the SALAD pairs.

## Preregistered coral re-fit: bleaching floor × chronic background × forcing `q` (2026-10-14)

Status: **directed by the project lead (2026-10-14), item 3; preregistered before implementation.** Runner-only (forced COTS); uses the opt-in ADRIA floor switch.

**Lead's reasoning.** With the bleaching floor raised, coral loss previously attributed to low-DHW bleaching (1–3.6 DHW, every few years) is removed. The background should be re-fitted, so that this loss is attributed to other chronic processes. A floor of 3 should also be tried.

### Setup

- **Model.** COTS forced from observed CPUE ÷ `q` (COTS dynamics off), observed `acrpoc` composition, RME group 5, CoralBlox growth untouched, κ 5, β 1.5.
- **Background.** Mellin table + μ: combined mortality `1 − (1 − m_cyc)^κ · (1 − m_mellin) · (1 − μ)`.
- **Forcing table.** For `q = 0.030`, the forcing table is scaled by 0.5 at run time and written into the run directory with its hash.

### Arms (24)

Floor {1 (default), 3, 4} × μ {0.03, 0.045, 0.06, 0.08} × forcing `q` {0.015, 0.030}.

**Reference.** `refit_f100_mu030_q015` must replay run 09 `bleach_mel3_k500_b150` exactly.

### Frozen objective and selection

- **Objective.** Coral `J`, unchanged: mean RMSE fast + slow (photo transects) + total (manta), 1992–2024.
- **Selection.** Per floor and `q`: μ at minimum `J`. A μ at a grid bound is reported as unidentified.
- **Floor choice.** The floor is evidence-set. Floor 4 (or 3) is preferred if its best `J` is within 5% of floor 1's best at the same `q`; otherwise this is reported.
- **Secondary reports:**
  - Lizard pre-outbreak fast cover (1992–96) and fast cover 2005–13, against photo transects;
  - 2017–19 total trough, against manta (survey-timing flags apply);
  - 2020–24 recovery.

### Expected consequences

- **Floor 3–4** removes the 1994–2010 low-DHW losses, so the best μ rises (about 0.045–0.06) to restore chronic loss. `J` should be within a few percent of floor 1.
- **Lizard's pre-outbreak fast cover** rises at floor 3–4, because the 1994 loss is gone.
- **`q = 0.030`** halves forced grazing, so the best μ rises further.
- **Risk.** A higher uniform μ removes fast and slow coral alike. The 2005–13 fast overshoot (about 2×) persists.

### Failures and rollback

**Failures:** an inexact replay; mortality outside [0, 1]; COTS dynamics active.

**Rollback:** floor unset, μ = 0.03, `q` = 0.015.

### Coral re-fit decision (2026-10-14)

[Run 12](runs/20261014T_12_coral_refit_1991_v1/metadata.toml), [audit](runs/20261014T_12_coral_refit_1991_v1/coral_refit_audit.json) and [plot](runs/20261014T_12_coral_refit_1991_v1/coral_refit.png). New audit script: `audit_1991_coral_refit.py`.

**Checks.**
- `refit_f100_mu030_q015` replays run 09 `bleach_mel3_k500_b150` exactly (0.0).
- COTS dynamics are zero in all forced arms, and cover stays in range.

**Result.** Best `J` per floor and forcing `q`:

| Floor | `q` 0.015: best μ / `J` | `q` 0.030: best μ / `J` |
| --- | --- | --- |
| 1 | 0.08 (bound) / **0.254** | 0.08 (bound) / 0.282 |
| 3 | 0.08 (bound) / 0.287 | 0.08 (bound) / 0.325 |
| 4 | 0.08 (bound) / 0.295 | 0.08 (bound) / 0.338 |

For reference, floor 1 with μ 0.03 at `q` 0.015 (the current configuration) has `J` 0.315.

- **μ is at the grid bound in every cell**, so it is unidentified. The coral model wants more chronic loss than 8%/yr, as in the earlier joint fit (6% at its bound). Raising μ improves slow cover and manta total, but worsens fast RMSE (0.065 → 0.077).
- **Floors 3 and 4 are 13–16% worse than floor 1, even with μ re-fitted.** That fails the preregistered 5% condition. The loss is in 1992–2012 (`J` 0.246 at floor 1 against 0.288 and 0.296).
- **Forcing at `q = 0.030`** (half the grazing) is worse at every floor.

**Why floors 3 and 4 lose (survey dates from the timing tracker).** The 9 m photo transects record a fast-coral **loss across the 1994 heat stress**:
- **Lizard:** 0.039 in October 1993 → 0.023 in October 1994 (−41%).
- **North Direction:** 0.064 → 0.036 (−44%), with almost no COTS (0.06 COTS/tow), so not grazing.
- **MacGillivray:** 0.031 → 0.036 (+16%).
- **Model, floor 1:** Lizard −73% (bleaching plus forced COTS at `q` 0.015), North Direction −50%, MacGillivray −80%.
- **Model, floor 4:** Lizard −11%, North Direction **+22%**.

So the data contradict "no mortality below 4 DHW" at these reefs in 1994 (3.3–3.5 effective DHW): a loss occurred. Floor 1 overshoots it where COTS are also forced, and floor 4 removes it.

**Decision.**
- Nothing is promoted, and the floor stays at 1 (the default).
- The floor switch stays opt-in. The evidence supports *some* mortality at about 3.3–3.5 DHW on these reefs. A weaker low-DHW response (rather than a hard floor) would need its own evidence.
- The best μ is unidentified (bound at 0.08). It is carried as a factor below, not adopted.

## Preregistered Lizard hard-coral initialisation × observation `q` × chronic background (2026-10-14)

Status: **directed by the project lead (2026-10-14), item 4; preregistered before implementation.** Runner-only.

### Evidence

- **Pratchett (2005)** measured reef-wide hard coral at Lizard (10 locations, 3 and 7 m belt transects): **32.2% in October 1996** and 22.9% in January 1999. At the same time, LTMP at 9 m recorded 12.7% (manta) and 10.7% (photo) in October 1996, and 10.6% and 12.2% in October 1998. LTMP under-represents Lizard's coral **2.2–2.5×**, where the starfish actually lived (sheltered back reef and crest at 3 m, exposed slopes at 7 m).
- **The model's Lizard cover** tracks LTMP manta: about 0.16 of site area in 1991, against manta 0.19 in October 1993.
- **RHIS (snorkel, 2–3 m, 2018–26)** shows Lizard's shallow hard and fast cover at about half its neighbours' (hard 0.10–0.18 against 0.27–0.37). That is after severe bleaching in 2016–17, 2020 and 2024, so it does not inform 1991, and it is **not** used to scale. It is reported as a caveat.
- **The other reefs** have no comparable reef-wide 1990s evidence, so the scale applies to Lizard's three sub-reefs only.

### Treatment

- **Lizard scale `c_L`.** The 1991 initial cover (all groups and sizes) at Lizard sites is multiplied by `c_L` ∈ {1, 1.5, 2}, keeping the observed `acrpoc` composition share. 2 is roughly Pratchett's ratio. Total habitable cover must stay ≤ 0.9.
- **Observation `q`** ∈ {0.015 (frozen), 0.030 (Pratchett pairs)}, applied consistently: seeds are `eb` CPUE ÷ `q`, and the trajectories are scored at the same `q` (the `q = 0.030` arms via `rescore_detectability.jl --q=0.030`).
- **Chronic background μ** ∈ {0.03 (current), 0.08 (re-fit bound)}.
- **Everything else** is the adopted baseline: working coral, floor 1, β 1.5, `Q = 20`, Ext×5, gradual Allee.

**Arms (12).** `wc_b150_q20_eb_lin_obs{015,030}_mu{030,080}_c{100,150,200}`.

**Reference.** `…_obs015_mu030_c100` must replay run 17 `wc_b150_q20_eb_lin` exactly.

### Expected consequences

- **`c_L` 2** gives Lizard about 2× the food, so its first wave rises (perhaps towards 40–55 adults/ha on the diver scale), its second wave rises, and troughs rise.
- **`q` 0.030** halves the seeds; first waves change little (run 17) and read 2× higher. Troughs also read 2× higher, so they risk the G2′ trough test.
- **μ 0.08** lowers coral everywhere, so waves and troughs both shrink.

### Report

- summed loss under the arm's own `q`;
- G2′ (trough gated on Lizard and MacGillivray), read at `A½` 0 and 10;
- Lizard peak adults against Pratchett's 40–55/ha;
- Lizard first and second waves;
- coral `J` components, with Lizard total cover against manta and Pratchett's 32% (1996).

### Failures and rollback

**Failures:**
- an inexact replay;
- Lizard total cover > 0.9 of habitable area;
- non-finite states.

**Rollback:** `c_L = 1`, `q` 0.015, μ 0.03.

### Lizard hard-coral decision (2026-10-14)

[Run 13](runs/20261014T_13_lizard_coral_1991_v1/metadata.toml), [audit](runs/20261014T_13_lizard_coral_1991_v1/lizard_coral_audit.json) and [plot](runs/20261014T_13_lizard_coral_1991_v1/lizard_coral.png). New audit script: `audit_1991_lizard_coral.py`.

Supporting runs:
- re-scores at [`q` 0.015](runs/20261014T_14_rescore_lizard_coral_q015_v1/metadata.toml) and [`q` 0.030](runs/20261014T_15_rescore_lizard_coral_q030_v1/metadata.toml);
- G2′ evaluations at [`q` 0.015](runs/20261014T_16_gate_g2prime_lizard_coral_q015_v1/gate_g2prime.json) and [`q` 0.030](runs/20261014T_17_gate_g2prime_lizard_coral_q030_v1/gate_g2prime.json) (trough gated on Lizard and MacGillivray).

**Checks.**
- `obs015_mu030_c100` replays run 17 `wc_b150_q20_eb_lin` exactly (0.0).
- External supply is identical in all arms, and the `q` 0.015 re-score control is exact.

**Result** (each arm read at its own `q`):

| Arm | Loss | Lizard first / second (COTS/tow) | Lizard peak adults (Pratchett 40–55) | Lizard total HC 1996 (Pratchett 0.32) | MacG / North Dir. second (obs 0.67 / 0.71) | G2′ at `A½` 10 (Lizard / MacG trough; bound 0.054) | Coral `J` |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `q` .015, μ .03, `c_L` 1 (ref) | 9.10 | 0.43 / 0.16 | 29 | 0.32 | 0.40 / 0.83 | pass (0.048 / 0.040) | 0.363 |
| `q` .015, μ .03, `c_L` 2 | 9.19 | 0.31 / 0.17 | 21 | 0.41 | 0.40 / 0.83 | pass | 0.400 |
| **`q` .015, μ .08, `c_L` 1** | **8.86** | 0.40 / 0.21 | 27 | 0.29 | **0.54** / 0.90 | fail by a hair (0.051 / 0.054) | **0.267** |
| `q` .030, μ .03, `c_L` 1 | **7.30** | 0.87 / 0.32 | 29 | 0.32 | 0.80 / 1.67 | fail (0.097 / 0.080) | 0.363 |
| `q` .030, μ .08, `c_L` 1 | 7.58 | 0.82 / 0.42 | 27 | 0.29 | 1.07 / 1.77 | fail | 0.267 |

**Against the predictions.**
- **Lizard hard coral: the opposite of what was predicted.** Raising Lizard's starting coral **lowers** its first wave (peak adults 29 → 24 → 21 at `c_L` 1 / 1.5 / 2). The model total in 1996 rises only 0.32 → 0.41, because the reef is near space-full.
- **Why.** Fast coral cannot grow in 1992–93. The 1994 heat stress removes a larger stock (0.111 → 0.045). The doubled slow coral then crowds fast-coral recovery (0.025 in 1996, against 0.042), so Lizard's starfish have **less** food. In this model, more *total* coral at a fixed composition does not feed COTS; only fast coral that survives into 1995–97 would.
- **`q` 0.030 (Pratchett, with consistent seeds).** As in the re-score, peaks double, the loss falls to its lowest value yet (7.30), and troughs double, failing G2′.
- **μ 0.08** (the coral re-fit's choice, made independently with forced COTS): coral `J` 0.363 → 0.267, the COTS loss improves (9.10 → 8.86), and second waves rise. This is a **cross-validated improvement** on both sides, but μ sits at its grid bound.

**Exploratory, post-hoc re-score (not preregistered; used for no decision).**
- [Runs 18–19](runs/20261014T_19_gate_g2prime_lizard_coral_q040_exploratory_v1/gate_g2prime.json): density-dependent detection with `q_max` anchored to Pratchett's outbreak pairs (0.04) and `A½` 10.
- **Result.** Loss 7.4–8.0, distinct waves (worst separation 0.06–0.10), and MacGillivray's second wave 0.85. But inter-wave CPUE is 0.13–0.15, about 7 adults/ha over 1999–2011, against the observed level of about 4.6.
- **The remaining excess** is the model's slow decline after the first peak: 12–20 adults/ha persist through 1999–2002, whereas Pratchett saw 46 → 3 per ha in seven months.

**Decision.**
- Nothing is promoted.
- A Lizard hard-coral scale is **rejected**: it lowers Lizard's wave.
- **For the lead:**
  - adopting μ 0.08 (better coral fit and COTS loss) would need an extended μ grid first, to identify it;
  - adopting a detection model with Pratchett's outbreak `q`.
- **The binding failure** is now the post-peak crash (item 4). Pratchett (2005) shows it was disease-associated, affected all size classes, and happened with coral still at 23%.


### Correction (2026-10-14): the 1994 fast-coral loss and North Direction's COTS

The coral re-fit decision above wrongly stated that North Direction's 1993 → 1994 fast-coral drop (−44%) happened "with almost no COTS (0.06 COTS/tow)". The figure 0.06 is the **report-year** 1994 value, which was surveyed in **October 1993**.

At the **October 1994** survey, the one showing the fast-coral drop, COTS were at outbreak level:
- North Direction 0.56 COTS/tow;
- Lizard 0.53;
- MacGillivray 0.29, where fast coral rose.

The 1993 → 1994 transect losses therefore coincide with the **early COTS outbreak**, as the lead judged. They are not evidence for bleaching mortality at 3.3–3.5 DHW.

The forced-COTS table is keyed by reporting year, so the forced runs apply the October-1994 densities in model year 1995, one year late. In model year 1994 they graze at October-1993 levels. That penalised the floor-3 and floor-4 arms in 1994: with low-DHW bleaching removed, nothing removed fast coral in model 1994.

The "floor rejected by the coral data" conclusion is **withdrawn**. It is re-tested in the re-fit v2 screen below, with forcing aligned to survey dates as a factor.


## Lead's decisions (2026-10-14, third round)

- **No crash mechanism.** The lead does not think one is needed; item 4 is closed.
- **Chronic background.** Extend the μ grid and identify the level.
- **Bleaching below 4 DHW.**
  - **Lead's position:** there is no evidence for bleaching mortality below 4 DHW, and the 1993 → 1994 transect fast-coral loss is much more likely to have come from the early COTS outbreak (confirmed by the correction above).
  - **Change wanted:** not necessarily a floor, but at least a strong reduction in low-DHW mortality, combined with a higher background.
- **Starting fast cover.** With no 1991 data, use the entire range of fast shares observed at each reef in all photo-transect data, limit it to plausible values by comparison with the manta tows, and allow the initialisation to be partially calibrated.

## Preregistered coral re-fit v2: low-DHW bleaching reduction × extended μ × forcing alignment (2026-10-14)

Status: **directed by the project lead; preregistered before implementation.**

### Mechanism (opt-in ADRIA switch, default identical)

- **Switches.** `ADRIA_BLEACHING_LOW_DHW_SCALE` (`s`, default 1) and `ADRIA_BLEACHING_LOW_DHW_THRESHOLD` (`D`, default 4 DHW).
- **Rule.** In `bleaching_mortality!`, if effective DHW < `D` at a location, that year's bleaching mortality fraction is multiplied by `s`. Mortality at or above `D` is unchanged.
- **How it differs from the floor.** It does not renormalise the tolerance distribution, so the response at high DHW is exactly the default. It is the lead's "strong reduction under 4 DHW".
- **Known simplification.** The stored bleaching DHW (the next year's lower truncation bound) still updates when mortality is scaled.
- **Test:** `test_bleaching_floor.jl` is extended to check that the default is identical, that mortality is scaled below `D`, and that it is unchanged above.

### Forcing alignment (new derived input)

- **Builder change.** `build_observed_cots_forcing.py` gains `--year-basis {report, survey}`. The default `report` reproduces the existing table (checked by hash). `survey` keys each LTMP survey by `floor(survey_date)`, so model year Y is grazed by the COTS observed during calendar year Y.
- **Scope.** This is a coral-forcing input only. Scoring stays at the reporting year (lead).

### Arms (31)

- **Low-DHW treatment `L`:** default, `s` = 0.25 (below 4), `s` = 0 (below 4).
- **μ:** {0.06, 0.08, 0.10, 0.12, 0.14}.
- **Forcing alignment:** {report, survey}.
- Plus `ref` = default, μ 0.03, report, which **must replay run 09 `bleach_mel3_k500_b150` exactly**.
- **Fixed:** κ 5, β 1.5, growth untouched, `q` 0.015 forcing, observed composition.

### Frozen objective and selection

- **Objective.** Coral `J`, unchanged.
- **μ.** Per (`L`, alignment), μ at minimum `J`; a value at the grid edge is reported as unidentified.
- **`L`.** Set by the lead's evidence position: `s` = 0.25 or 0 is preferred if its best `J` is within 5% of the default's at the same alignment.
- **Alignment.** Reported. It changes no scoring.

**Expected.**
- μ is identified in 0.08–0.12.
- Under survey alignment, the 1994 fast-coral loss is reproduced by COTS rather than bleaching, so `s` = 0–0.25 is no worse than the default.

## Preregistered partial calibration of the 1991 starting fast share (2026-10-14), stage 2 after re-fit v2

- **Plausible range per reef.** The photo-transect fast share (*Acropora* + Pocilloporidae over hard coral, 9 m), taken from all years in which that reef's LTMP manta total cover was **at or above its own median**, i.e. healthy, pre-outbreak-like states consistent with 1991. Years with degraded cover are excluded.
- **Levels.** The minimum, 25th, 50th and 75th percentiles and maximum of that range for Lizard, MacGillivray and North Direction. Other reefs (including Eyrie, which has no transects) take the mean of the three reefs' values at the same level.
- **Application.** Applied as in run 12: totals and within-block structure preserved.

**Stage 2a (5 arms + reference).**
- Arms at each level, all reefs together, on the coral configuration selected in v2 and the adopted COTS baseline.
- The reference arm `acrpoc` replays the matching v2-configuration arm.

**Stage 2b (partial calibration).**
- Each calibration reef's level is chosen by **its own** smoothed frozen loss, which is approximately separable because seeds and coral are reef-specific.
- One confirmation arm then runs the per-reef choices together.
- **Adoption condition:** the confirmation arm must not worsen G2′ (trough gated on Lizard and MacGillivray) relative to stage 2a's `acrpoc` arm.

**Report.**
- per-reef chosen levels;
- first and second waves;
- G2′;
- coral `J` (expected to worsen for the transects in 1993–95 by construction, since they are not the 1991 state);
- Lizard peak adults against Pratchett's 40–55/ha.

### Coral re-fit v2 decision (2026-10-14)

[Run 20](runs/20261014T_20_coral_refit_v2_1991_v1/metadata.toml), [audit](runs/20261014T_20_coral_refit_v2_1991_v1/coral_refit_v2_audit.json) and [plot](runs/20261014T_20_coral_refit_v2_1991_v1/coral_refit_v2.png).

**New tools.**
- **ADRIA** (opt-in, default identical): `ADRIA_BLEACHING_LOW_DHW_SCALE` and `ADRIA_BLEACHING_LOW_DHW_THRESHOLD` in `bleaching_mortality!` and `run_model`. `test_bleaching_floor.jl` now passes 13/13.
- **Forcing builder:** `build_observed_cots_forcing.py --year-basis survey`. The default `report` rebuild reproduces the existing table hash exactly. The survey-year table is in `sandbox/data/derived/observed_cots_forcing_1991_2024_survey_year/`.
- **Audit:** `audit_1991_coral_refit_v2.py`.

**Operational notes.**
- The first audit attempts failed before any result was read:
  - `grid.align` resolved to pandas' `DataFrame.align` method, not the column, so the column was renamed `alignment`;
  - a plot path key was wrong.
- Partial audit outputs written by those attempts were deleted before the rerun.

**Checks.**
- The reference replays run 09 exactly (0.0).
- COTS are off in all forced arms, and cover stays in range.

**Result** (best μ per low-DHW treatment and forcing alignment; frozen `J`; the current configuration is 0.315):

| Low-DHW bleaching | Reporting-year forcing: μ / `J` | Survey-year forcing: μ / `J` |
| --- | --- | --- |
| default (`s` = 1) | **0.10** / 0.249 | 0.10 / 0.253 |
| **`s` = 0.25 below 4 DHW** | **0.12** / 0.259 | 0.12 / 0.259 |
| `s` = 0 below 4 DHW | 0.14 (bound) / 0.289 | 0.12 / 0.294 |

- **μ is identified.** It is now interior: 0.10 with default bleaching, 0.12 with low-DHW mortality cut to a quarter. As the lead expected, a stronger low-DHW reduction is compensated by a higher chronic background.
- **`s` = 0.25** is within 5% of the default (`J` 0.259 against 0.249), so it is **selected** under the frozen rule, consistent with the lead's evidence position.
- **`s` = 0** is 16% worse at every μ. Its manta fit worsens (2017–19 total 0.05 against 0.071 observed), so a complete removal is not supported, but a strong reduction is.
- **Survey-year forcing** leaves `J` unchanged. But it drives Lizard's 1992–96 fast coral to about 0.004 (observed 0.035), because the forced densities (`q = 0.015`) are applied in the year they were observed, which is the 1993–94 outbreak build-up. This is the forcing-`q` question again: Pratchett's diver densities imply the forced adults are about 2× too high.

**Decision.**
- Nothing is promoted.
- **Selected coral configuration ("cv2")** for the next COTS screens:
  - CoralBlox growth untouched;
  - Mellin + **12%/yr** chronic;
  - κ 5, β 1.5;
  - bleaching mortality ×**0.25 below 4 DHW**.
- Reporting-year forcing is retained for coral fits, consistent with scoring.

### Starting fast-share partial calibration decision (2026-10-14)

**Inputs.** [Range (run 21)](runs/20261014T_21_initial_share_range_v1/initial_share_range.json), built by `build_initial_share_range.py` (table in the git-ignored `sandbox/data/derived/initial_fast_share_range/`). Healthy-year photo-transect fast share (years with manta ≥ the reef's median):

| Reef | Range | Median | Healthy years kept |
| --- | --- | --- | --- |
| Lizard | 0.147–0.453 | 0.308 | 11 of 24 |
| MacGillivray | 0.050–0.387 | 0.129 | 12 of 25 |
| North Direction | 0.191–0.627 | 0.326 | 12 of 24 |
| Others (mean of the three) | 0.129–0.489 | | |

**Runs.**
- Stage 2a: [run 22](runs/20261014T_22_initial_share_1991_v1/initial_share_audit.json).
- Stage 2b confirmation: [run 25](runs/20261014T_25_initial_share_1991_v2/initial_share_audit.json).
- Re-scores 23 and 26; G2′ evaluations 24 and 27 (trough gated on Lizard and MacGillivray).
- New tools: `audit_1991_initial_share.py`. The runner gains the `initial_share` experiment and `OWEN_INITIAL_SHARE_PICKS`.

**Checks.**
- `…_ref` (old coral, observed composition) replays run 17 `wc_b150_q20_eb_lin` exactly (0.0) in both runs.
- 1991 totals and external supply are identical across arms, and the re-score controls are exact.

**Result** (cv2 coral configuration unless stated):

| Arm | Summed loss | Lizard / MacG / North Dir. 1991 share | Lizard first / second (obs 1.92 / 0.74) | Lizard peak adults (Pratchett 40–55) | MacG second (obs 0.67) | North Dir. first (obs 0.56) | G2′ `A½` 10 trough, Lizard / MacG (≤ 0.054) | Coral `J` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ref (old coral, `acrpoc`) | 9.10 | 0.21 / 0.13 / 0.25 | 0.43 / 0.16 | 29 | 0.40 | 0.92 | 0.048 / 0.040 (pass) | 0.363 |
| cv2 `acrpoc` | 8.58 | 0.21 / 0.13 / 0.25 | 0.48 / 0.24 | 32 | 0.60 | 1.36 | 0.078 / 0.069 | 0.239 |
| cv2 min | 8.50 | 0.15 / 0.05 / 0.19 | 0.42 / 0.26 | 28 | 0.62 | 0.92 | 0.054 / 0.061 | 0.253 |
| cv2 q50 | 8.21 | 0.31 / 0.13 / 0.33 | **0.79 / 0.35** | **53** | 0.58 | 2.33 | 0.186 / 0.070 | 0.244 |
| cv2 q75 / max | 11.87 / 16.62 | up to 0.45 / 0.39 / 0.63 | 1.16–1.79 / 0.43–0.57 | 77–119 | 0.52–0.58 | 3.3–3.6 | 0.26–0.48 / 0.08–0.62 | 0.255–0.282 |
| **cv2 pick** (Lizard min, MacG q50, North Dir. q25) | **7.97** | 0.15 / 0.13 / 0.31 | 0.44 / 0.24 | 29 | 0.58 | 1.96 | 0.056 / 0.070 | 0.260 |

**Against the expectations.**
- **The cv2 coral configuration helps both models.** Coral `J` improves (0.363 → 0.239) and the COTS loss improves (9.10 → 8.58), with larger second waves, MacGillivray's near observed. Troughs rise (Lizard and MacGillivray about 0.07–0.08 at `A½` 10).
- **Within the plausible range**, the median share puts Lizard's first wave at 53 adults/ha, inside Pratchett's diver-observed 40–55.
- **But Lizard's own frozen loss picks the minimum share.** The frozen loss (linear `q` 0.015) penalises the higher troughs that a larger first wave brings more than it rewards the closer peak height.
- **Per-reef separability is only approximate.** In the pick arm, Lizard's own loss (1.77) is higher than in the uniform-min arm (1.47): North Direction's larger first wave at q25 changes Lizard's supply.
- **Confirmation.** Loss 7.97, the lowest yet under the frozen observation model. G2′ is not worse than cv2 `acrpoc`; both fail the trough test at `A½` 10 (Lizard 0.056, MacGillivray 0.070), so the adoption condition holds.

**Caveats for the lead.**
- (1) **North Direction's first wave (1.96 against 0.56 observed) is barely penalised**, because the smoothed observed series has no 1995 peak. The score partly ignores this overshoot.
- (2) **Lizard's chosen level conflicts with Pratchett's diver density** (29 against 40–55 adults/ha). The higher level that matches Pratchett is rejected by a loss built on `q = 0.015`. The observation-model choice (`q`) and Lizard's initial share are coupled.

**Decision.**
- Nothing is promoted.
- **Recommended working configuration** (lead to confirm):
  - **Coral:** cv2 (low-DHW bleaching ×0.25 below 4 DHW, Mellin + 12%/yr, κ 5, β 1.5, growth untouched).
  - **1991 fast share:** partially calibrated (Lizard 0.15, MacGillivray 0.13, North Direction 0.31, others 0.25).
  - **Everything else:** the adopted baseline (`eb` seeds, `Q` 20, Ext×5, gradual Allee).
  - **Score:** loss 7.97, coral `J` 0.260.
- G2′ fails on Lizard and MacGillivray troughs at `A½` 10 (0.056 and 0.070 against 0.054), and passes at neither `A½` 0 nor 10.


## Lead's decisions (2026-10-14, fourth round) and close-out plan

**Decisions.**
1. **Working configuration confirmed:** cv2 coral plus the partially calibrated 1991 fast shares (run 25 `wc_b150_q20_eb_lin_cv2_pick`; loss 7.97, coral `J` 0.260).
2. **Re-score it under Pratchett's `q` (0.030).**
3. **Scoring amendments:**
   - **North Direction's** first peak is scored against the **raw** observed series (its 1995 peak counts).
   - **Eyrie is not penalised** for a missing observed first peak: sparse surveys most likely missed it, and the IDW / outbreak-prediction surface would show it.
4. **Close-out.**
   - **Goal:** a short list of the most plausible candidate parameter sets to carry into simulations.
   - **Problem to fix:** peak timing, which has drifted later.
   - **Then:** promote the best-performing options and start on COTS removal from the culling data.

**Timing diagnosis (from stored raw scores).** Modelled first peaks in reporting-year terms:

| Arm | Lizard | MacGillivray | North Direction | Eyrie |
| --- | --- | --- | --- | --- |
| run 15 | 1998 | 1997 | 1998 | 2000 |
| run 17 baseline | 1997 | 1997 | 1997 | 2000 |
| cv2 `acrpoc` | 1998 | 1997 | 1998 | 2000 |
| cv2 pick | 1997 | 1997 | **1999** | 2000 |
| cv2 q50 | 1999 | 1997 | 1999 | 2000 |

- **Observed** (raw, reporting year): Lizard 1998, MacGillivray 1996, North Direction 1995. In survey terms that is October 1997, October 1995 and October 1994.
- **Gap.** North Direction and MacGillivray first waves come 2–4 years late. The lag grows as coral and starting fast share rise.

The plan to work through is recorded in `CURRENT_STATE.md` ("Close-out plan").


## Preregistered frozen score v2 (2026-10-14), close-out step A1

**Why.** This implements lead decision 3 of the fourth round.
- **North Direction.** The `smoothed_3y` observed series has no 1995 peak, because 3-year calendar smoothing of 0.06 / 0.56 / 0.12 falls below the peak criteria. A North Direction first wave of any height or timing is therefore nearly unpenalised. The raw series has the 1995 (0.56) and 2013 (0.71) peaks.
- **Eyrie.** Surveys ran every three years (1995, 1998 at 0.13, 2001 at 0.09). A 1999–2000 first wave would fall between surveys, so Eyrie's observed series has only the 2013 peak. Under v1, every arm that produces an Eyrie first wave pays the full peak-count penalty (1.25).

**Definition.** Score v2 is a wrapper. `cots_cycle_metrics.jl` and its weights, peak rules and observation model are unchanged.
1. **Per-reef observation treatment:**
   - Lizard, MacGillivray and Eyrie: `smoothed_3y` (as v1);
   - North Direction: `raw`.
2. **Eyrie, unobserved first wave.** Model peaks earlier than (earliest observed peak − peak tolerance), i.e. before 2007, are dropped from the peak count. Then `peak_count_penalty_v2 = min(|n_sim_v2 − n_obs| / n_obs, 1)` and `loss_v2 = loss_v1 − 1.25 × (peak_count_penalty_v1 − peak_count_penalty_v2)`.
   - Such peaks can never be matched under the 6-year tolerance, so no other term changes.
   - The series terms (RMSE, bias and lag correlation against the 1998 and 2001 surveys) are kept.
   - An early peak cannot stand in for a missing 2013 peak: with only an early peak, v2 is worse than v1.
3. **Summed loss v2** = Lizard + MacGillivray + Eyrie (`smoothed_3y`, Eyrie adjusted) + North Direction (`raw`).
4. **Reporting.** v1 (summed `smoothed_3y`) and v2 are reported side by side for every candidate. G2′ is unchanged; its v2 peak count (7 observed peaks: 2 + 2 + 2 + 1) is reported alongside.

**Implementation.**
- `score_v2.jl` holds the wrapper as pure functions; it is tested by `test_score_v2.jl`.
- `rescore_detectability.jl` adds the columns `peak_count_penalty`, `loss_v2` and `score_v2_primary` to every row. Its `loss` column is unchanged, so the `A½ = 0` replay check still applies.
- `evaluate_gate_g2prime.py` adds `summed_loss_v2` and the v2 matched-peak counts when those columns are present.

**Control.** For Lizard and MacGillivray rows, and for Eyrie rows with no model peak before 2007, `loss_v2 == loss` exactly. Re-scoring run 25 must reproduce its stored v1 losses (`A½ = 0`, `q` 0.015).

**Expected consequences**, computed from the stored run-26 per-reef scores before implementation:

| Arm | v1 summed | v2 summed (expected) |
| --- | --- | --- |
| `ref` (run-17 baseline) | 9.10 | about 6.85 |
| cv2 `acrpoc` | 8.58 | about 7.33 |
| cv2 pick | 7.97 | about 8.30 |
| cv2 q50 | 8.21 | about 8.61 |

- **Eyrie.** All four arms lose the 1.25 penalty.
- **North Direction.** The working configuration's overshoot (1.72 against 0.56) and lateness (1999) are now scored: 2.30 → 3.88.
- **So v2 is expected to re-rank** the working configuration below `acrpoc` and the run-17 baseline. The North Direction share pick (q25), made under v1, should be revisited in B2.

**Failure criteria.**
- The implementation fails if the control does not replay exactly, or if any v2 value differs from the hand computation above beyond rounding.
- A v2 ranking that differs from v1 is a finding, not a failure.

**Rollback.** v2 lives in separate columns and files; v1 is untouched.


## Score v2 results on stored runs (2026-10-14, close-out A1 / A2 re-score)

**Runs.**
- Run `20261014T_28_rescore_score_v2_initial_share` re-scores run 25 at `q` 0.015. Its control replays with a difference of exactly 0. The v2 sums match the preregistered hand computation: `ref` 6.847, cv2 `acrpoc` 7.326, cv2 pick 8.302, cv2 q50 8.608.
- Run 29 is G2′ at `q` 0.015.
- Run 30 is a score-only re-score at `q` 0.030. Its seeds were still inverted at 0.015, so it is not consistent. Run 31 is G2′ at `q` 0.030.

**Summed loss** (`A½ = 0`; v1 = summed `smoothed_3y`):

| Arm | v1 @ 0.015 | **v2 @ 0.015** | v2 @ 0.030 (score only) | North Direction v2 @ 0.015 [first, second peak] |
| --- | --- | --- | --- | --- |
| `ref` (run-17 baseline) | 9.10 | 6.85 | 6.48 | 1.85 [1997, 2012] |
| cv2 `acrpoc` | 8.58 | 7.33 | 8.83 | 2.93 [1998, 2012] |
| **cv2 min** (all reefs at minimum share) | 8.50 | **6.45** | 6.90 | 1.99 [1997, 2012] |
| cv2 pick (working configuration) | **7.97** | 8.30 | 10.90 | 3.88 [1999, 2013] |
| cv2 q50 | 8.21 | 8.61 | 10.04 | 4.32 [1999, 2013] |

- **Re-ranking.** v2 re-ranks the arms, as expected. North Direction's first-wave overshoot and lateness under the q25 share dominate the working configuration's v2 loss.
- **Per-reef v2 best** among the run-25 share levels:
  - Lizard: min (1.47);
  - MacGillivray: q50 (1.25);
  - North Direction: min (1.99; first peak 1997).
- **G2′.** At `A½ = 10`, only `ref` passes at `q` 0.015. cv2 min fails only on MacGillivray's trough (0.061 against 0.054; Lizard 0.054).

## Late first waves: flux diagnosis (2026-10-14, close-out B1)

Source: run 25 `flows.csv` and `mechanisms.csv`, North Direction and MacGillivray, 1991–2001.

**Settlement is set externally.**
- It is dominated by external supply (Ext×5), identical in every arm: North Direction 38, 47, 95, 104, 116, 131 per ha per year over 1992–97.
- Local fecundity matters only from 1995. The 1991 seeds (North Direction 0.7 adults/ha) are negligible against it.

**So the rise is set by settlement plus maturation**, and is the same in every share arm.
- North Direction adults reach 15 per ha in 1995 and 37 in 1996.
- The observed peak of 0.56/tow equals 37/ha at `q` 0.015. The model reaches it about a year late (observed report year 1995, i.e. the October 1994 survey).
- MacGillivray reaches about 21–25 per ha in 1996, against its observed 1996 peak. Its model peak is 1997.

**The peak year is set by food exhaustion.** Adults keep rising until fast coral runs out.

| North Direction arm | Fast coral 1995 | Fast coral gone (≤ 0.002) | Peak adults/ha (year) |
| --- | --- | --- | --- |
| cv2 min | 0.086 | 1998 | 61 (1998; scored peak 1997) |
| `ref` | 0.072 | 1998 | 61 (1997) |
| cv2 pick (q25 share) | 0.119 | 1998–99 | **131 (1999)** |

- The q25 share therefore gives both the overshoot (1.72/tow against 0.56) and the late peak.
- The observed population instead fell within a year of its peak (0.56 → 0.12). With no crash mechanism (lead's decision), the remaining timing levers are:
  - **less food**, i.e. a lower share (earlier exhaustion and a lower peak);
  - **a faster rise:** a lower predator-pit `Q` (run 11: `Q` 10 brought the mean first peak about one year earlier than `Q` 30);
  - **the observation `q`:** at 0.030, 0.56/tow is 19 adults/ha, reached in 1995.

## Preregistered timing factorial (2026-10-14), close-out B2 with the A2 consistent-`q` arms

**Hypothesis.** The late and high North Direction first wave comes from extra food (the q25 share) plus a slow rise. A lower `Q` and a lower North Direction share bring the first peaks earlier and lower. A consistent `q` 0.030 makes the rise visible earlier in CPUE terms (but, from run 30, doubles troughs).

**Change.** No model change. New runner experiment `timing`: cv2 coral, `eb` seeds, Ext×5, gradual Allee. The factorial is 3 × 2 × 2 = 12 arms, named `wc_b150_q<Q>_eb_lin_obs<q>_cv2_<share>`:
- **predator-pit `Q`** {10, 15, 20} COTS/ha/yr;
- **share set**, with other reefs at q50:
  - `pick`: Lizard min, MacGillivray q50, North Direction q25 (the working configuration);
  - `pnd`: the same with North Direction at min;
- **observation `q`** {0.015, 0.030}. The `eb` seeds are inverted with the arm's `q` (linear, 6:3:1 preload), as in `lizard_coral` `obs030`. The runner scores at 0.015, and `obs030` arms are read via `rescore_detectability.jl --q=0.030`.

**Control.** `wc_b150_q20_eb_lin_obs015_cv2_pick` must replay run 25's `wc_b150_q20_eb_lin_cv2_pick` exactly (trajectories, flows and mechanisms). External immigration and 1991 total cover must be identical across arms.

**Expected consequences.**
- **Share `pnd`:**
  - North Direction first peak in 1997 (against 1999), with a lower height (about 0.9 against 1.7/tow at `q` 0.015);
  - North Direction v2 loss about 2.0 (against 3.9);
  - other reefs nearly unchanged (they are connected, so not identical).
- **`Q` 10:** first peaks on all reefs about 0–1 year earlier and higher, with higher troughs. A G2′ trough failure is the expected cost.
- **`q` 0.030:**
  - first-wave CPUE about double, so peak-height penalties rise where the model already overshoots;
  - Lizard's height fits better (Pratchett 40–55/ha corresponds to 1.2–1.6/tow against the observed 1.92);
  - trough failures.

**Read-out.** Each arm at its own `q`:
- summed v2 loss (primary) and v1, both at `A½ = 0`;
- G2′ at `A½` 0 and 10 (troughs gated on Lizard and MacGillivray);
- first-peak year per reef against the observed year, in both reporting-year and survey-year terms;
- first- and second-wave heights;
- Lizard peak adults/ha against Pratchett's 40–55;
- coral `J`.

**Selection rule for C4** (the timing fix):
1. Take the arm with the lowest summed v2 loss among arms that match 7/7 v2 peaks and pass G2′ criteria 1–3.
2. Arms within 0.10 of the minimum are ranked by the mean absolute first-peak year error over Lizard, MacGillivray and North Direction; then by the smaller trough.
3. If no `q` 0.030 arm is selected, the best `q` 0.030 arm is still carried as C2 (Pratchett `q`).

**Failure criteria.** The reference does not replay, or an arm fails to run. A screen in which no arm improves first-peak timing over `cv2_pick` is a finding: the timing then needs a mechanism (excluded by the lead for the crash), and is reported rather than tuned.

**Bounds and evidence.**
- `Q` 10–20 is within the run-11 grid. The Allee threshold is held at its evidence value.
- The shares stay within each reef's healthy-year photo-transect range.
- `q` 0.030 is the Poisson fit to the Pratchett 2005 pairs.

**Rollback.** Runner-only experiment; nothing changes in the defaults.


## Timing factorial results (2026-10-14), close-out B2 + A2

**Runs.**
- Run `20261014T_32_timing_1991_v1`: 12 arms, no failures.
- Re-scores: run 33 (`q` 0.015) and run 34 (`--q=0.030`).
- G2′: runs 35 and 36.
- Audit: `runs/20261014T_32_timing_1991_v1/timing_audit.json` (`audit_1991_timing.py`).

**Controls.**
- `wc_b150_q20_eb_lin_obs015_cv2_pick` replays run 25's `cv2_pick` exactly (trajectories, flows and mechanisms: 0.0).
- External immigration and 1991 total cover are identical across arms.
- The run-33 `A½ = 0` control replays with a difference of 0.0.

**Results** (each arm at its own `q`, `A½ = 0`; first peaks in reporting years; observed: Lizard 1997 smoothed, MacGillivray 1997 smoothed, North Direction 1995 raw):

| Arm (`Q`, `q`, share set) | v1 | **v2** | coral `J` | First peaks (Lizard, MacGillivray, North Direction) | North Direction first wave (/tow) | Troughs at `A½ 10` (Lizard, MacGillivray) |
| --- | --- | --- | --- | --- | --- | --- |
| **20, 0.015, `pnd`** | 8.50 | **6.44** | **0.251** | 1997, 1997, **1997** | 0.93 | 0.056, 0.070 |
| 15, 0.015, `pnd` | 8.97 | 6.89 | 0.253 | 1997, 1996, 1997 | 1.03 | 0.057, 0.074 |
| 10, 0.015, `pnd` | 9.51 | 7.36 | 0.254 | 1996, 1996, 1997 | 1.16 | 0.056, 0.083 |
| 20, 0.015, `pick` (reference) | 7.97 | 8.30 | 0.260 | 1997, 1997, 1999 | 1.96 | 0.056, 0.070 |
| 20, 0.030, `pnd` | 9.26 | 8.57 | 0.250 | 1997, 1997, 1997 | 1.98 | 0.113, 0.139 |
| 20, 0.030, `pick` | 10.73 | 10.56 | 0.261 | 1997, 1997, 1999 | 4.01 | 0.113, 0.139 |

Remaining arms are in the audit.
- Every arm matches 7/7 v2 peaks and passes G2′ criteria 1–3.
- None passes G2′ criterion 4 (troughs) at `A½ 10`.

**Findings.**
1. **North Direction share is the timing lever.** At its minimum (0.19, against 0.31 at q25), the North Direction first wave peaks in 1997 at 0.93/tow (against 1999 at 1.96).
   - North Direction v2 loss falls from 3.88 to 2.02, and the summed v2 loss from 8.30 to 6.44.
   - Lizard and MacGillivray are essentially unchanged; coral `J` improves.
   - It is still two reporting years late (observed 1995, i.e. the October 1994 survey) and 1.7× too high.
2. **Lower `Q` brings Lizard and MacGillivray one year earlier** (1996), but raises v2 loss (Lizard loss 1.77 → 2.46 at `Q` 10) and troughs. Not selected. As in run 11, `Q` trades timing against trough and height.
3. **The consistent `q` 0.030 is worse under v2** (8.57 against 6.44).
   - It helps Lizard's height (Lizard loss 1.77 → 1.30). It overshoots MacGillivray and North Direction and doubles troughs.
   - Halving the seeds barely changes adult dynamics: Lizard peak 29.7 against 29.5 adults/ha. External supply dominates, so the 1991 state is not the lever.
4. **Lizard's first wave stays at about 29 adults/ha** in every arm, below Pratchett's 40–55. This needs the Lizard share (C3), not `Q` or `q`.

**Selection (preregistered rule).**
- **C4 (timing fix):** `wc_b150_q20_eb_lin_obs015_cv2_pnd`.
- **C2 (best Pratchett-`q` arm):** `wc_b150_q20_eb_lin_obs030_cv2_pnd`.
- Four arms improve first-peak timing over the reference: `pnd` at `Q` 20 and 15, at both `q` values.

**Decision.** Nothing is promoted. C4 is recommended to replace C1 as the working configuration (lead to decide). It is C1 with North Direction's 1991 fast share at the minimum of its healthy-year range.


## Lead's decisions (2026-10-14, fifth round)

1. **Working configuration changed.** North Direction's 1991 starting fast-coral share is lowered to the minimum of its healthy-year range (0.191, from 0.307). Run 32 `wc_b150_q20_eb_lin_obs015_cv2_pnd` (loss v2 6.44, coral `J` 0.251) becomes the working configuration. The previous one (run 25 `cv2_pick`) is kept as a comparison candidate.
2. **Larger Lizard first outbreak.** Test Lizard's starting fast share at the median of its healthy-year range (0.308, against 0.147) with both manta detectabilities: the current `q` 0.015, and Pratchett's `q` 0.030 with seeds inverted at the same `q`.
3. **Culling.**
   - (a) **Stages:** culled COTS under 15 cm (cohort 1) come off model juveniles; 15 cm and over (cohorts 2–4) come off adults.
   - (b) **Location:** spatially explicit. Cull-site centroids are fuzzy-matched to model sites on the same reef, not spread across the reef in proportion to modelled COTS.
   - (c) **Year:** LTMP reporting-year convention (September–December culls count toward the following year).
   - (d) **Timing:** culls before January (September–December) are removed before spawning; then spawning occurs; then January–August culls are removed.
4. **Replicates for promotion:** vary both the connectivity-sequence seed and ADRIA's RNG seed (at least 5 pairs).
5. **Process:** decision requests to the lead are written in plain ecological terms, without internal labels.

## Preregistered candidate set and replicate-seed run (2026-10-14), close-out B3 + B4 (replicates) + C3 of the plan

**Candidates.** All use `eb` seeds, Ext×5 and gradual Allee. Other reefs keep the q50 share.

| ID | Arm suffix | Coral | 1991 fast share (Lizard / MacGillivray / North Direction) | `Q` | Observation `q` (seeds inverted at it) |
| --- | --- | --- | --- | --- | --- |
| W (working, was C4) | `cv2_pnd` | cv2 | min / q50 / min | 20 | 0.015 |
| C1 (previous working) | `cv2_pick` | cv2 | min / q50 / q25 | 20 | 0.015 |
| C2 (Pratchett `q`) | `cv2_pnd` | cv2 | min / q50 / min | 20 | 0.030 |
| C3a (Lizard median) | `cv2_pndl50` | cv2 | **q50** / q50 / min | 20 | 0.015 |
| C3b (Pratchett-consistent Lizard) | `cv2_pndl50` | cv2 | **q50** / q50 / min | 20 | 0.030 |
| C5 (structural control) | `ref` | run-17 coral (Mellin + 3%/yr, default bleaching) | `acrpoc` | 20 | 0.015 |

**Seeds.**
- **Pairs:** five (connectivity-sequence seed, ADRIA RNG seed) pairs. Pair `s1` is the archived pair (20260930, 90972102). Pairs `s2`–`s5` are (20261101, 11101), (20261102, 22202), (20261103, 33303) and (20261104, 44404).
- **Per arm:** both seeds are set via `COTS_CONNECTIVITY_SEED` and `ADRIA_RNG_SEED`. The runner's Julia `Random.seed!` stays at 20260930 and is recorded.
- **Arm names:** `wc_b150_q20_eb_lin_obs<q>_s<k>_<suffix>`. That is 6 × 5 = 30 arms; the runner experiment is `candidates`.

**Controls (`s1`, exact replay).**
- `obs015_s1_cv2_pick` replays run 25 `cv2_pick`.
- `obs015_s1_cv2_pnd` replays run 32 `obs015_cv2_pnd`.
- `obs030_s1_cv2_pnd` replays run 32 `obs030_cv2_pnd`.
- `obs015_s1_ref` replays run 25 `ref`, i.e. run 17.

**Expected.**
- **C3a:** Lizard's first wave at about 50 adults/ha (Pratchett 40–55), peaking about 1999 with a higher loss at `q` 0.015. Its Lizard height penalty rises (0.75/tow against the observed 1.92 raw).
- **C3b:** Lizard's first wave of about 1.5/tow against the observed 1.92, i.e. a better height. MacGillivray and North Direction overshoot, and troughs double, as in run 32.
- **Replicates:** peak years shift by at most about 1 year between seeds. Single-seed rankings that differ by less than the between-seed sd are not distinguished.

**Read-out.**
- Per candidate at `s1`: v1 and v2 summed loss, G2′ at `A½` 0 and 10, first-peak years, first- and second-wave heights, Lizard peak adults/ha, and coral `J`.
- Over the five seeds, per reef: mean ± sd of v2 loss; peak years, heights and prominences of the first two peaks; the inter-peak period; and the trough mean.
- AGENTS gate 3 is assessed as "peak timing, height, prominence and period better than C1 in at least 4 of 5 seeds without worse coral `J`".

**Failure criteria.** A control fails to replay, or any arm errors.

**Not decided here.** Which candidates go forward is the lead's decision. The held-out tests and ±25% sensitivity (AGENTS gates 4–5) follow on the chosen set.

## Preregistered observed-culling switch (2026-10-14), close-out C2

**Hypothesis.** Culling (171,948 COTS on 72 domain reefs, mostly 2012–15 and 2021–26) lowered the second wave and post-2015 populations where it was intense, both directly and through reduced larval supply from heavily culled source reefs. The calibration so far has ignored it.

**Change.** COTSMod gains a pure, opt-in removal step: `apply_removals!(model, juveniles_ha, adults_ha) -> (removed_juveniles_ha, removed_adults_ha)`.
- Removals are capped at the COTS present: `min(request, N[2])` and `min(request, N[3])`.
- Adult age classes and `large_adults` are scaled by the same surviving fraction.
- Recruits (`N[1]`) are never culled.
- Nothing runs unless it is called.

**ADRIA adapter.**
- A new env `COTS_CULL_REMOVALS_CSV` points to a table with columns `site_id, report_year, phase, juveniles_ha, adults_ha` (COTS per ha of model site area).
- **Unset** gives the established path, unchanged.
- **When set:**
  - for each annual step whose calendar year Y equals `report_year`: `phase = pre` rows (September–December dives of Y − 1) are applied immediately before the COTS step (grazing, survival, spawning, maturation);
  - `phase = post` rows (January–August dives of Y) are applied immediately after the COTS step and before larval dispersal;
  - removed amounts per site-year are returned as `cots_culled` (juveniles and adults; additive output);
  - unknown site IDs or negative values are errors.
- Removal before the step reduces that year's spawning, grazing and maturation. Removal after it does not affect this year's spawning, but does lower the end-of-year (logged) population.

**Derived table.** `build_cull_removals.py` v2 gains:
- **fuzzy matching:** each cull-site centroid's counts are spread over model sites on the same reef with Gaussian weights `area × exp(−d² / 2σ²)` for sites within 3σ, falling back to the nearest site;
- `σ` = 210 m, half the median nearest-neighbour spacing of cull sites (422 m; model-site spacing 274 m). Sensitivities are 0 m (nearest site only) and 420 m;
- **adapter tables** (`removals_<σ>.csv`) with phase from the dive month, cohort 1 as juveniles and cohorts 2–4 as adults.
- Mass balance is checked for every σ.

**Tests.**
- **COTSMod unit tests:**
  - zero request is the identity;
  - removals are capped at available;
  - age classes and `large_adults` are scaled consistently;
  - recruits are untouched;
  - non-finite or negative requests are errors.
- **ADRIA adapter tests:**
  - table parsing and validation;
  - an unset env gives an identical run;
  - a pre / post ordering check on a two-site toy.

**Experiment** (after the candidate run):
- the working configuration (W) at `s1` × {no cull, observed cull at σ 210} for 1991–2024, plus σ 0 and 420 sensitivities;
- C3 and C5 with and without culling, if the lead keeps them.
- The no-cull arm must replay W exactly.

**Expected consequences.**
- North Direction's 2013 peak is lowered: 417 culled, about 20/ha reef-wide in 2013, concentrated at 2 model sites.
- Eyrie's 2021–24 population is lowered (1–3/ha/yr).
- Lizard and MacGillivray show only small direct effects.
- Internal immigration is reduced in 2013–16 and 2022–26 through heavily culled sources (Startle, Irene, Eye, 14-133). This is damped because external supply dominates settlement.
- Coral is slightly higher after 2013 on culled reefs.

**Read-out.**
- Per reef: second-wave height and timing; 2016–24 trough; v2 loss and G2′; coral `J`.
- **Removal diagnostics:** requested vs removed (the unmet share where the model has fewer COTS than were culled), and internal immigration changes.

**Failure criteria.**
- The replay or a test fails, or mass balance fails.
- If more than 25% of requested removals are unmet on any calibration reef, the spatial matching or the modelled density is inconsistent with the cull data. That is reported to the lead before any interpretation.

**Bounds.** The removals are observed data; no parameter is tuned.

**Rollback.** The env is unset by default.


## Candidate set and replicate-seed results (2026-10-14), close-out B3 + B4 (replicates) + C3

**Runs.**
- Run `20261014T_37_candidates_1991_v1`: 30 arms, no failures.
- Re-scores: run 38 (`q` 0.015) and run 39 (`--q=0.030`).
- G2′: runs 40 and 41.
- Audit: `candidates_audit.json` (`audit_1991_candidates.py`).

**Controls.** All four seed-pair-1 arms replay their prior arms exactly (0.0): C1 → run 25, W and C2 → run 32, C5 → run 25 `ref`.

**Over 5 seed pairs** (each candidate at its own `q`; mean ± sd):

| Candidate | Loss v2 | Loss v1 | Coral `J` | Mean first-peak error (yr) | Lizard peak adults/ha | Lizard trough at `A½ 10` | MacGillivray trough at `A½ 10` | Seeds passing G2′ at `A½ 10` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **W** (working) | **6.79 ± 0.79** | 8.98 ± 0.83 | 0.250 | **0.67 ± 0** | 29 ± 2 | 0.056 | 0.075 | 0 |
| C1 (previous working) | 9.81 ± 1.44 | 10.22 ± 1.85 | 0.256 | 1.33 ± 0 | 29 ± 2 | 0.056 | 0.075 | 0 |
| C2 (W at `q` 0.030) | 7.75 ± 0.87 | 8.58 ± 0.64 | 0.250 | 0.67 ± 0 | 29 ± 2 | 0.112 | 0.145 | 0 |
| C3a (Lizard median share) | **6.72 ± 0.70** | 8.92 ± 0.77 | **0.228** | 1.33 ± 0 | **56 ± 7** | **0.191** | 0.075 | 0 |
| C3b (C3a at `q` 0.030) | 7.78 ± 1.05 | 8.64 ± 0.86 | 0.228 | 1.33 ± 0 | 57 ± 7 | 0.391 | 0.146 | 0 |
| C5 (run-17 coral, structural control) | 6.84 ± 0.49 | 9.22 ± 0.55 | **0.363** | 0.60 ± 0.15 | 28 ± 3 | **0.047** | **0.037** | **4 of 5** |

Observed trough bound: 0.054. Per-reef replicate statistics (peak years, heights, prominence, period, trough) are in `candidates_replicates.csv`.

**Findings.**
1. **Seeds.** Seed variation in loss v2 is 0.5–1.4 (sd). W, C3a and C5 are not distinguishable on loss (6.72–6.84). Peak timing is robust: first-peak years vary by ≤ 1 year, and W's first-peak error is identical in all 5 seeds.
2. **W against C1.** W is better than C1 on loss, timing and coral in 4 of 5 seeds, so it meets the preregistered replicate criterion (AGENTS gate 3). C3a (5/5) and C3b (4/5) also meet it.
3. **The Lizard median share (C3a)**:
   - gives Pratchett's density (56 adults/ha);
   - gives the best coral fit (`J` 0.228);
   - but the Lizard first peak moves to 1999 (observed 1997/98);
   - and Lizard's inter-wave trough is 3.5× the observed bound (0.191 against 0.054).

   A larger Lizard outbreak leaves too many COTS between waves.
4. **Pratchett's `q` (C2, C3b)** doubles troughs and raises loss in every seed. It does not help.
5. **The old coral (C5)** is the only configuration that keeps the troughs within the observed bound (4 of 5 seeds pass G2′). Its coral fit is much worse (`J` 0.363 against 0.250).
   - The cv2 coral changes (more chronic background mortality, less low-DHW bleaching) leave more food after the first wave, and so more COTS between waves.
   - This trade-off between coral fit and COTS troughs is the main unresolved tension.

**Decision.** Nothing is promoted. The candidate set for the lead to choose from is W, C3a and C5. C1, C2 and C3b are dominated.


## Observed-culling screen results (2026-10-14), close-out C2

**Implementation.**
- **COTSMod:** `apply_removals!` (opt-in, capped, scales adult age classes and large adults). The COTSMod test suite passes: 380/380, including the new removal tests.
- **ADRIA:** `load_cots_cull_removals` and `apply_cots_cull!` in `ecosystem/cots.jl`; hooks in `scenario.jl` (pre-spawning before the COTS step, post-spawning after it); output `cots_culled_log`. Adapter test `test_cull_removals.jl`: 16/16.
- **Derived tables:** `build_cull_removals.py` v2 (`sandbox/data/derived/cots_cull_removals_lizard_v0_3_v2/`, git-ignored) with fuzzy matching at σ 0 / 210 / 420 m. Mass balance is preserved at every σ.

**Run.**
- Run `20261014T_42_culling_1991_v1`: the working configuration at seed pair 1 × {no cull, σ 0, 210, 420}, no failures.
- The no-cull arm replays run 32 exactly (0.0).
- Re-score: run 43. G2′: run 44. Audit: `culling_audit.json`.

**Removals requested vs applied** (COTS counts, σ 210):

| Scope | 2012–16 requested | Unmet | 2021–24 requested | Unmet |
| --- | --- | --- | --- | --- |
| Lizard | 631 | 4% | 573 | 51% |
| MacGillivray | 35 | 0% | 12 | 7% |
| North Direction | 429 | 0% | 29 | 50% |
| Eyrie | 141 | 0% | **3,770** | **88%** |
| Whole domain | 106,118 | **57%** | 16,714 | **70%** |

**Preregistered failure criterion triggered.** More than 25% of requested removals are unmet on calibration reefs, for all σ on Eyrie and Lizard. Interpretation is therefore limited to the following facts, and the matter goes to the lead.
1. **2012–16, on the calibration reefs:** removals are almost fully applied. The model's second-wave densities there are consistent with the culls.
2. **2012–16, domain-wide:** 57% unmet. Heavily culled reefs outside the calibration set (Startle, Irene, Eye, U/N 14-133) had more COTS at cull sites than the model holds there.
3. **2021–24:** most removals cannot be applied. The cull programme removed 16,700 COTS in the domain (Eyrie 3,770, about 7/ha), but the model has no outbreak in the 2020s.
   - LTMP manta recorded little at Eyrie in these years (0.00 in 2021, 0.09 in 2024).
   - So either a 2020s outbreak was under-detected by manta (consistent with the hyperdepletion found in the paired diver/manta data), or culling concentrates on local aggregations far denser than the reef average.

**Effects where removals apply** (σ 210 against no cull):
- North Direction's 2013 peak falls from 0.82 to 0.55/tow (observed 0.71), and its v2 loss rises from 2.02 to 2.16.
- Internal immigration falls 2–4% on the calibration reefs over 2012–24. It is damped because external supply dominates settlement.
- Other reefs, coral `J` and G2′ are essentially unchanged. Summed v2 goes from 6.44 to 6.57.

**Decision.** Nothing is promoted. The switch works and is opt-in. Whether to include culling in the candidate simulations, and how to read the 2020s cull numbers, is for the lead.


## Lead's decisions (2026-10-14, sixth round)

1. **Simulations:** carry three configurations into simulation:
   - the working configuration (W);
   - W with Lizard at its median starting fast share (C3a);
   - the run-17 coral structural control (C5).
2. **Coral-fit vs between-outbreak trade-off:** run a small test of intermediate coral settings before any formal sensitivity analysis.
3. **2020s cull numbers:** compare cull catch with manta tows at the same reefs and times.
   - The lead regards the 2020s culls, against near-zero manta counts, as further evidence that manta tows detect COTS poorly at low densities.
   - So more modelled COTS between outbreaks is acceptable.
   - Between-outbreak (trough) results are therefore also reported under the stronger low-density under-detection setting (`A½ = 30`). The G2′ definition is not changed.
4. **Culling is off** in calibration comparisons until item 3 is resolved.

## Preregistered intermediate coral settings (2026-10-14)

**Hypothesis.** The cv2 coral settings fit coral better than the run-17 settings, but leave more food after the first wave, and so more COTS between waves. Intermediate settings may give a better compromise.

**Change.** No model change. New runner experiment `coral_mix`:
- the working configuration (1991 shares: Lizard min, MacGillivray q50, North Direction min; Q 20; `q` 0.015; seed pair 1);
- chronic background coral mortality μ {0.03, 0.06, 0.09, 0.12} per year, with the same expression as cv2 and run 17: `1 − (1 − m_cyclone)^5 · (1 − m_Mellin) · (1 − μ)`;
- × bleaching mortality multiplier below 4 DHW {1, 0.5, 0.25};
- = 12 arms, plus the run-17 coral control arm (C5);
- arm names `wc_b150_q20_eb_lin_obs015_s1_m<μ>_l<scale>_cv2_pnd`, and `wc_b150_q20_eb_lin_obs015_s1_ref`.

**Controls.**
- (μ 0.12, ×0.25) must replay W exactly (run 37 `obs015_s1_cv2_pnd`).
- The `ref` arm must replay C5 (run 37 `obs015_s1_ref`).
- (μ 0.03, ×1) is the run-17 coral with W's starting shares. Its difference from C5 isolates the share effect.

**Expected.** Lower μ and a larger low-DHW multiplier lower coral cover after 1998 and so lower between-wave COTS. Coral `J` worsens toward 0.36.

**Read-out.** Per arm:
- loss v2;
- coral `J`;
- G2′, with between-wave levels for Lizard and MacGillivray at `A½` 10 and 30;
- first-peak timing and heights.

**Selection rule.** Report the trade-off front (coral `J` against the worse of the Lizard and MacGillivray troughs at `A½ 10`). An arm is a candidate replacement for W only if:
- coral `J` ≤ 0.30;
- loss v2 ≤ W + 0.5;
- and its troughs pass at `A½ 10`, or (as a weaker standard, given decision 3) at `A½ 30`.

The lead decides.

**Failure criteria.** A control does not replay, or an arm errors.

**Rollback.** Runner-only.

## Preregistered cull-catch vs manta comparison (2026-10-14)

**Question.** Are cull-dive catches consistent with manta tows that strongly under-detect at low density, as hyperdepletion predicts?

**Data.** CSIRO workbook, data licence, local only; only summaries go into the protocol.
- **Pairs:** `Cull` dives (count, bottom time) and `Manta` tows (count per tow) on the same domain reef within the same voyage. LTMP manta reef-years (`reef_cots.csv`) are a secondary pairing by reporting year.

**Measures.**
- cull catch per unit effort (COTS per diver-hour);
- manta COTS per tow;
- the share of reef-visits with positive cull catch but zero manta COTS;
- fits of manta CPUE against cull CPUE, linear vs hyperdepleted (`m = a · C² / (C + C½)`), compared by AIC on the visit means (Poisson likelihood on manta counts given tows).

**Expectation under hyperdepletion:**
- many reef-visits with cull catch but zero manta COTS;
- the hyperdepleted fit preferred, with `C½` well above the low-density cull rates.

**Caveats.** Cull divers target aggregations. Cull CPUE also saturates at high density (handling time), which biases toward a linear or concave shape. A hyperdepleted result is therefore conservative.

**Output.** `analyse_cull_vs_manta.py` writes derived tables to the git-ignored `sandbox/data/derived/cull_vs_manta_lizard/`. Only summary statistics are recorded here.


## Intermediate coral settings results (2026-10-14)

**Runs.**
- Run `20261014T_45_coral_mix_1991_v1`: 13 arms, no failures.
- Re-score: run 46. G2′: run 47. Audit: `coral_mix_audit.json`.

**Controls.**
- μ 0.12 / ×0.25 replays W exactly.
- `ref` replays C5 exactly.
- The runner checks that the μ 0.12 and μ 0.03 mortality cubes reproduce cv2 and run 17 bit for bit.

**Results** (seed pair 1; trough ratio = worse of Lizard and MacGillivray between-wave CPUE divided by the observed bound of 0.054):

| μ | Bleaching multiplier below 4 DHW | Loss v2 | Coral `J` | Trough ratio at `A½ 10` | Trough ratio at `A½ 30` | G2′ at `A½ 10` | Lizard peak adults/ha |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.12 | ×0.25 (**W**) | 6.44 | 0.251 | 1.29 | 0.67 | fail (passes at `A½ 30`) | 29 |
| 0.12 | ×1 | 5.93 | 0.238 | 1.35 | 0.72 | fail | 22 |
| **0.09** | **×0.25** | **6.94** | **0.279** | **0.96** | 0.48 | **pass** | 30 |
| 0.09 | ×1 | 6.32 | 0.262 | 1.11 | 0.57 | fail | 23 |
| 0.06 | ×0.25 | 7.54 | 0.326 | 0.85 | 0.41 | pass | 31 |
| 0.03 | ×1 | 6.51 | 0.371 | 0.74 | 0.36 | pass | 23 |
| C5 (run-17 coral and shares) | – | 6.85 | 0.363 | 0.88 | 0.44 | pass | 29 |

**Findings.**
1. **Chronic background mortality μ is the lever for the trade-off.** Each step down from 0.12 lowers between-wave COTS by about 10–15% and worsens coral `J` by about 0.03–0.05. Peak timing is unaffected.
2. **μ 0.09 with ×0.25 is the compromise.** It passes all G2′ criteria at `A½ 10` (Lizard 0.050, MacGillivray 0.052) with coral `J` 0.279, much better than C5's 0.363.
   - Its loss v2 (6.939) misses the preregistered cut-off (W + 0.5 = 6.938) by 0.001, far inside the seed sd (0.79). It is therefore reported, not selected.
   - It has not yet been replicated across seeds.
3. **The bleaching multiplier below 4 DHW.** Full mortality (×1) gives a slightly better coral fit and loss than ×0.25 at every μ, but a smaller Lizard first wave (22–23 adults/ha).
   - This conflicts with the lead's evidence-based decision (no bleaching mortality evidence below 4 DHW).
   - Per AGENTS, an evidence-based setting is not changed to improve the fit. It is reported only.
4. **Under the stronger low-density under-detection (`A½ 30`), W passes G2′.** This is relevant to the lead's sixth-round decision 3.

## Cull-catch vs manta results (2026-10-14)

**Run.** `analyse_cull_vs_manta.py` → `sandbox/data/derived/cull_vs_manta_lizard/` (git-ignored; licensed data). Summary statistics only:

**Pairing.**
- 1,073 cull visits are paired with 3,231 manta tows over the same cull site within ±30 days.
- 752 visits are paired within ±14 days.
- Domain reefs, 2013–2026.

**Detection at the visit level.**
- Of the visits whose cull caught COTS, **81–82% had zero manta COTS** over the same site.
- By sextile of cull catch (median 0, 0.5, 1.4, 2.7, 4.5 and 9.3 COTS per diver-hour), pooled manta was 0.018, 0.060, 0.076, 0.064, 0.108 and 0.173 COTS per tow (±30 days). The share of visits with zero manta COTS was 95%, 88%, 85%, 83%, 79% and 70%.

**Fit.**
- Poisson fits of manta counts give a linear relation with a small background (`0.0146 · C + 0.035` per tow, ±30 days).
- A power law adds nothing (`k` 1.14 at ±30 days, 0.96 at ±14 days; ΔAIC +1.8 and +2.0 against the linear-with-background fit).
- So, relative to cull catch, manta is a weak, roughly linear signal with no extra low-density curvature.

**Against LTMP manta (by reef and reporting year).**
- Lizard's cull catch in 2024–26 was 1.6–2.3 per diver-hour (2013 outbreak year: 3.7) while LTMP manta read 0.00–0.05 (2013: 0.74).
- Eyrie in 2024: cull 2.2 per diver-hour against manta 0.09 (2013: 10.1 against 0.85).
- North Direction in 2024–25: cull 1.2 against manta 0.00.

**Interpretation.**
- Cull divers kept finding COTS at roughly a third to two-thirds of 2013 catch rates in the 2020s, while manta read near zero.
- Either manta under-detects strongly at low density, or cull catch per hour stays high as density falls (divers target aggregations).
- Both readings imply that reef populations in the 2020s and between waves are larger than manta suggests. This supports the lead's view that more modelled COTS between outbreaks is acceptable.
- The data cannot separate the two mechanisms without an absolute density anchor (the SALAD diver surveys are the candidate anchor).
- The model's 2020s Lizard population (about 2 adults/ha) is likely too low.

**Decision.** Nothing is promoted. Culling stays off pending the lead.


## Lead's decisions (2026-10-14, seventh round)

1. **Background mortality 9%/yr.** Replicate μ 0.09 (bleaching ×0.25 below 4 DHW) across the 5 seed pairs. If confirmed, it replaces the run-17 coral control (C5) in the simulation set, which becomes W, C3a and C6.
2. **Between-wave (trough) test.** The test formally uses the stronger low-density under-detection by manta (`A½ = 30`). `A½ = 10` is reported alongside.
3. **SALAD as the density anchor.** SALAD diver densities anchor the conversion of cull catch to density.
   - **Priority:** the model must reproduce the build-up of the fifth outbreak in 2020–2026 seen in the SALAD and cull data. These are more accurate than manta, which only starts to rise in 2025/26.
4. **Promotion gates.** Held-out and sensitivity gates follow on the final set.

## Fifth-outbreak build-up: diagnosis and diver-anchored targets (2026-10-14)

**Model** (W at seed pair 1; run 37 flows and mechanisms):
- COTS stay below about 1.6 adults/ha on Lizard, North Direction and Eyrie through 2017–2024.
- **Three limits act together:**
  1. **External larval supply collapses.**
     - Lizard: 98–174 per ha per year (2014–16), then 6–14 (2018–24).
     - North Direction: 169–295, then 8–22.
     - The supply is built from the GBR-wide outbreak-probability surface `gbrPredsAdj_20262408.csv`, i.e. the probability that manta CPUE > 0.22 on each upstream reef.
     - Mean probability on northern GBR reefs (12–17.5° S): 0.25 in 2012–14, but **0.018–0.038 in 2018–2025, the lowest in the record**.
     - The boundary forcing therefore inherits manta's low-density blindness.
  2. **Starvation.** Fast coral is near zero from 2015 to about 2020 (second wave, 2016–17 bleaching, cyclones). Adult food survival is 0.05–0.11 until 2021, so maturing adults die.
  3. **The predator pit** removes 2–6 of the 7–26 recruits per ha per year.

**Diver-anchored targets** (`build_2020s_density_targets.py` → git-ignored `sandbox/data/derived/density_targets_2020s/`; summaries only).
- **SALAD–cull pairs:** 21 SALAD reef-months paired with cull dives on the same reef within ±60 days.
- **Conversion:** density ≈ **7.9 × (cull COTS per diver-hour)^0.55**.
  - The power fit is better than linear (AIC 350 against 393).
  - The concave shape is consistent with cull divers targeting aggregations.
  - The Poisson fit ignores overdispersion, so its intervals are too narrow.

Targets (adults/ha, reporting year; SALAD where surveyed, otherwise converted cull catch):

| Reef | 2020–21 | 2022–23 | 2024–25 | Model W, 2020–24 |
| --- | --- | --- | --- | --- |
| Lizard | 3.4 | 6–8 | 9 | 0.4–1.6 |
| North Direction | 3.4–3.6 | 2.5–3.7 | 11–13 | 0.05–1.25 |
| Eyrie | 0–2 | 17–45 | 12–16 | 0.1–1.3 |
| Martin | 7 | 7 | 9–11 | (not output) |

- SALAD itself misses some COTS (published detectability 50–82%), so the targets are conservative.
- The model's years end in 2024. 2025–26 needs forcing extended (DHW, cyclones).

## Preregistered diagnostic: 2017–2024 external-supply counterfactual (2026-10-14)

**Question.** If more larvae arrive from upstream in the 2020s, can the model build toward the diver-anchored densities? Or do starvation and the predator pit still block it?

**Change.** Labelled counterfactual switch, opt-in.
- New ADRIA env: `COTS_EXTERNAL_SUPPLY_WINDOW = "<start>:<end>:<factor>"`. It multiplies the evidence-derived external supply (potential and pelagic) by `factor` in the calendar years `start..end`.
- Unset gives the established path.
- This is a diagnostic only, not a candidate mechanism. The supply itself would later be rebuilt from evidence.

**Arms** (W at seed pair 1, `q` 0.015):
- window 2017–2024 × factor {1, 3, 10, 30};
- window 2019–2024 × factor 10.

Factor 1 must replay W exactly.

**Expected.**
- If supply is the main limit, adults rise in 2021–24 toward 3–13/ha, roughly in proportion to the factor.
- If starvation dominates, the rise is delayed until fast coral recovers (about 2021–22), and is damped.

**Read-out.**
- Per reef: adults/ha in 2020–24 against the targets; the first year above 3/ha;
- the flux decomposition (settlement, juvenile predation, maturation, food survival);
- coral in 2020–24;
- loss v2 and G2′, to check that the 1992–2016 fit is not changed (the window starts after the second wave).

**Failure criteria.** Factor 1 does not replay, or an arm fails.

**Rollback.** The env is unset by default.


## Candidate rerun with μ 0.09 (C6) across 5 seed pairs (2026-10-14)

**Runs.**
- Run `20261014T_48_candidates_1991_v2`: 35 arms, no failures.
- Re-scores: runs 50 and 51. G2′: runs 52 and 53.
- Audit: `candidates_audit.json`, extended with C6 and the `A½ 30` columns.

**Controls.** All 30 original arms replay run 37 exactly, and C6 at seed pair 1 replays run 45 (93 table checks, maximum difference 0.0).

**C6** (W with μ 0.09), over 5 seeds:
- loss v2 7.14 ± 0.85;
- coral `J` 0.279;
- troughs at `A½ 10`: Lizard 0.050 ± 0.001, MacGillivray 0.055 ± 0.005 (bound 0.054);
- G2′ passes in 3/5 seeds at `A½ 10` and 3/5 at `A½ 30`.

**Under the adopted trough standard (`A½ 30`)**, G2′ passes in 5/5 seeds for W, 3/5 for C6 and 1/5 for C5.

**Decision.** C6 is confirmed as the trough-respecting member: its troughs are at the bound on average, with a much better coral fit than C5. Per the lead's conditional decision, it replaces C5. **The simulation set is W, C3a and C6.**

## External-supply counterfactual results (2026-10-14; labelled diagnostic)

**Implementation.**
- ADRIA `COTS_EXTERNAL_SUPPLY_WINDOW` (`parse_cots_supply_window`, `cots_supply_window_factor` in `cots.jl`; hook in `scenario.jl`). Unit test `test_supply_window.jl`: 12/12. The cull adapter test still passes (16/16).
- Run `20261014T_49_supply_window_1991_v1`: 5 arms, no failures. Factor 1 replays W exactly. Re-score: run 54. G2′: run 55. Audit: `supply_window_audit.json`.

**Results** (W at seed pair 1). "2020s error" is the mean log((model + 0.5) / (target + 0.5)) over 18 reef-years in 2020–24 with diver-anchored targets.

| Supply | 2020s error (bias, mean absolute) | Lizard adults 2020–24 | North Direction 2020–24 | Eyrie 2020–24 | Loss v2 | G2′ at `A½ 30` |
| --- | --- | --- | --- | --- | --- | --- |
| ×1 (W) | −1.57, 1.59 (about 5× low) | 0.5, 0.4, 0.8, 1.6, 1.2 | 0.5, 0.1, 0.1, 1.3, 0.7 | 0.1, 0.2, 0.4, 0.9, 1.3 | 6.44 | pass |
| ×3, 2017–24 | −0.63, 0.95 | 1.1, 0.4, 2.2, 8.5, 2.3 | 1.9, 0.3, 3.8, 17.9, 14.9 | 0.2, 0.9, 1.8, 5.3, 5.2 | 6.54 | pass |
| ×10, 2017–24 | 0.00, 1.15 | 3.4, 1.5, 12.8, 35.0, 2.9 | 7.3, 1.7, 19.0, 43.9, 2.8 | 0.8, 2.1, 0.9, 5.8, 1.5 | 8.89 | pass |
| ×30, 2017–24 | +0.42, 1.16 | 11.6, 5.8, 10.2, 14.3, 4.0 | 23.6, 6.1, 16.7, 11.9, 4.2 | 3.0, 7.1, 2.6, 6.0, 4.4 | 12.08 | fail |
| ×10, 2019–24 | −0.27, **0.89** | 0.5, 6.7, 14.9, 16.1, 2.0 | 0.5, 1.0, 2.8, 11.3, 2.5 | 0.1, 6.8, 7.3, 26.3, 5.4 | 7.16 | pass |

Targets: Lizard about 3 → 9, North Direction about 3.5 → 11–13, Eyrie about 2 → 12–45 adults/ha.

**Findings.**
1. **External supply is the main limit on the 2020s build-up.** Restoring it (×3–10) raises 2020–23 populations into the diver-observed range. With factors up to 10 the 1992–2016 fit is unchanged, because the window starts after the second wave.
2. **The model booms and busts where the observations show a sustained build-up.**
   - With extra supply, adults peak in 2023 (Lizard 16–35/ha) and crash in 2024 (2–3/ha) as food runs out. Total coral falls to 0.02–0.04 at ×10–30 under COTS grazing plus the 2024 bleaching.
   - SALAD shows a sustained, moderate build-up (Lizard 6–9 adults/ha from 2022 to 2025).
   - Juvenile predation removes 100–130 per ha over 2017–24, near the saturation cap (`Q` 20 per year).
3. **The manta-based loss v2 rises as the 2020s fit improves** (6.44 → 8.89 at ×10). Manta reads near zero in 2020–24, so the frozen score penalises a build-up that the diver data show. The 2020s need a diver-anchored score term.
4. **The model ends in 2024.** The 2025–26 targets (Eyrie 34/ha and South Direction 23/ha in 2026) need extended forcing.

**Decision.** Nothing is promoted. Rebuilding the 2020s supply from evidence, scoring the 2020s against diver data, and the boom-bust dynamics are for the lead.


## Lead's decisions (2026-10-14, eighth round)

1. **2020s supply.**
   - An interim fitted 2019–2026 external-supply multiplier is acceptable as a labelled stand-in.
   - The outbreak-prediction model cannot be re-run.
   - **Preferred routes:**
     - (a) **Initiation-region build-up.** The domain is in the initiation region, where several years of successful local recruitment, under favourable hydrodynamics and larval nutrition (upwelling), should raise densities. The lead would rather help this along than fiddle with the external supply.
     - (b) **Correct the prediction layer.** Use the cull records and SALAD data to correct the outbreak-prediction layer by spatial interpolation.
2. **2020s score.** A diver-anchored 2020s term is added (new score version). Manta stays only as a weak check after 2019.
3. **Boom-bust dynamics.**
   - The post-2017 coral recovery was unprecedented, and the model needs that fast-growing coral for COTS to thrive as they did.
   - The 2024 bleaching had little impact on Lizard and MacGillivray; plenty of hard coral remained.
4. **Extend to 2026.**
   - **DHW:** extract from RIMREP DMS using the lead's BleachingMortality code.
   - **Cyclones:** 2025–26 cyclones were not influential and can be ignored for now.
   - **Future need:** a way to convert real or simulated cyclone tracks to the 4 m wave-height metric.

## Coral 2012–2024: model against observations (2026-10-14)

W at seed pair 1 (run 37) against LTMP manta hard coral (HC) and photo-transect fast coral (Acropora + Pocilloporidae):

| Reef | Manta HC, 2019 → 2024 | Model total (site), 2019 → 2023 → 2024 | Photo fast coral, 2019 → 2023 → 2024 | Model fast coral |
| --- | --- | --- | --- | --- |
| Lizard | 0.07 → 0.21 | 0.07 → 0.135 → **0.077** | 0.009 → 0.093 → 0.062 | 0.032 → 0.109 → **0.024** |
| MacGillivray | 0.07 → **0.47** | 0.10 → 0.165 → 0.11 | 0.011 → 0.053 → 0.064 | 0.045 → 0.079 → 0.023 |
| North Direction | 0.09 → **0.58** | 0.09 → 0.108 → 0.09 | 0.014 → 0.068 → 0.110 | 0.031 → 0.041 → 0.024 |

- **Domain DHW** (NOAA CRW annual maximum, reef mean) at these reefs: 2020 5.5–6.6; 2021–22 3.9–4.4; 2024 5.9–6.2.
- **The model's coral is wrong in two ways after 2017:**
  1. **Recovery is too slow,** especially at MacGillivray and North Direction.
  2. **Bleaching losses are too large** in 2020 and 2024, where the observed coral kept rising.
- **Consequence:** the model's fast coral, the COTS food supply, collapses in 2024 just as the build-up should be accelerating. Together with the low external supply, this explains the boom-bust in run 49.
- **Candidate causes, to be diagnosed, not tuned yet:**
  - the bleaching scalar β 1.5 applied uniformly in every year;
  - full mortality at or above 4 DHW;
  - no thermal acclimation after the 2016–17 events;
  - the constant chronic μ 0.12/yr (identified on the 1992–2016 period), which drags on recovery;
  - NOAA DHW bias at these reefs (local logger corrections exist in BleachingMortality: `build_local_first_dhw_correction.py`).

## Plan for the 2020s fifth-outbreak build-up (2026-10-14)

1. **Score v3** (preregistered below), so that the 2020s can be scored at all.
2. **Forcing to 2026.**
   - DHW 2025–26 from DMS (`domain_building/extract_dhw_2025_2026.py`, validated against the domain's 2016–24 DHW).
   - Cyclones 2025–26 set to zero (lead).
   - Also needed: the Mellin background table, the COTS connectivity-year sampling, and the external-supply outbreak years (predictions end in 2025).
   - These would be assembled as a new domain version (v0.4); v0.3 is not modified.
3. **Post-2017 coral diagnosis.** A forced-COTS coral screen for 2017–2024(26): the bleaching response (β, the low-DHW scale, a post-2017 acclimation or DHW-correction counterfactual) × chronic μ, against manta and photo data. The aim is the observed recovery and the small 2024 impact at Lizard and MacGillivray.
4. **Supply routes, compared under score v3:**
   - the interim fitted 2019–26 multiplier (labelled);
   - (a) the initiation-region mechanism, i.e. year-varying larval survival or retention from an evidence index (e.g. chlorophyll-a / upwelling from DMS, `extract_sst_chla_features.py`);
   - (b) the corrected prediction layer (outbreak probabilities adjusted by interpolated SALAD/cull residuals) feeding the external supply.
5. **Then:** held-out and sensitivity gates on the simulation set; culling back on once the 2020s densities are plausible.

## Preregistered frozen score v3 (2026-10-14)

**Why.** Manta under-detects the 2020s build-up, so v2 penalises a model that matches diver densities (run 49).

**Definition.** A wrapper; v1 and v2 are unchanged.
1. **Manta part:** score v2 (North Direction raw; Eyrie's unobserved first wave), computed on 1992–2019 only (model and observations truncated at 2019).
2. **Diver part:** per reef with diver-anchored targets (`build_2020s_density_targets.py`, reporting years 2020–2024; extended to 2026 when the forcing is):

   `D_reef = mean over target years of |log((A_model + 0.5) / (T + 0.5))|`

   where `A_model` is adults/ha and `T` is the target. The 0.5 adults/ha offset keeps near-zero years finite.
3. **Combined:** `loss_v3_reef = v2_to2019_reef + 1.5 × D_reef`, and the summed v3 loss is the sum over the four scored reefs.
   - A 2× error then costs about 1.0 per reef; the current 5× error costs about 2.4 (comparable to a missed peak).
4. **Manta after 2019** is reported only, as the 2020–24 mean model CPUE against observed. It is not in v3.

**Implementation.**
- `score_v3.jl`: pure functions, with unit tests in `test_score_v3.jl`.
- `rescore_detectability.jl` adds `loss_v2_to2019`, `diver_error` and `loss_v3`; the evaluator adds `summed_loss_v3`.
- The targets file is licensed and git-ignored. Its hash is recorded, and if it is absent the v3 columns are missing (NaN).

**Expected.**
- W: v3 is about v2(to 2019) + 4 × 1.5 × about 1.6, i.e. about +9.
- The ×10 2019–2024 supply arm improves by about 4 relative to W.

**Control.** All existing v1 and v2 columns are unchanged, and the `A½ = 0` replay still holds.


## Score v3: correction and first results (2026-10-14)

**Correction to the preregistered definition** (made before any interpretation): the manta part is cut at **2021**, not 2019.
- **Why:** Eyrie has no manta survey between 2015 and 2021. A 2019 cut ended its common scoring window in 2015, so its 2014 model peak became undetectable and its loss jumped from 2.68 to 5.38.
- **Why 2021:** it is the first post-2019 survey year on every scored reef. With the 2021 cut the manta part equals the full-period values (7.76 summed before the Eyrie adjustment) and is identical across the supply arms, so manta no longer penalises the 2020s build-up.
- **Diver term:** still scored from 2020.
- **Superseded output:** run 56 used the 2019 cut and is marked SUPERSEDED. Run 57 is the corrected re-score; run 58 is its G2′.
- **Implementation:** `score_v3.jl` (`SCORE_V3_LAST_MANTA_YEAR = 2021`, `SCORE_V3_FIRST_DIVER_YEAR = 2020`), with tests in `test_score_v3.jl` (8/8). Re-score columns: `loss_v2_manta`, `diver_error`, `loss_v3`. Evaluator columns: `summed_loss_v3`, `summed_diver_error`, `summed_loss_v2_manta`. All v1 and v2 columns are unchanged (difference 0.0).

**Supply counterfactual under v3** (`A½ = 0`, seed pair 1):

| Supply | v2 | **v3** |
| --- | --- | --- |
| ×1 (W) | 6.44 | 15.66 |
| ×3, 2017–24 | 6.54 | 12.04 |
| ×10, 2017–24 | 8.89 | 15.36 |
| ×30, 2017–24 | 12.08 | 19.09 |
| ×10, 2019–24 | 7.16 | **11.80** |

- v3 rewards bringing the 2020s up to diver densities.
- v3 penalises the 2023 boom and 2024 crash: ×10 from 2017 is no better than W.

## 2025–2026 DHW extension (2026-10-14)

**Builder:** `domain_building/extract_dhw_2025_2026.py`. NOAA CRW daily 5 km DHW, from the public DMS zarr, at the nearest cell to each domain reef centroid; calendar-year maximum. Output: git-ignored `sandbox/data/derived/dhw_2016_2026_lizard/`.

**Validation against the domain DHW, 2016–2024** (1,017 reef-years):
- correlation 0.9996;
- mean absolute difference 0.018;
- 94.5% within 0.1 DHW;
- maximum difference 1.19.

The method reproduces the domain forcing.

**Extension** (domain reefs, annual maximum DHW):
- 2025: mean 5.47 (range 3.0–8.6);
- 2026: mean 4.37 (2.9–6.2; data to 12 October 2026).

Both are moderate-to-high heat years. With the current bleaching response (which over-predicts the 2020 and 2024 losses at Lizard and MacGillivray), they would collapse the modelled fast coral, so the coral diagnosis (plan step 3) has to precede the 2025–26 runs.


## Lead's decisions (2026-10-15, ninth round)

1. **Coral heat sensitivity.**
   - The bleaching multiplier (β, `ADRIA_DEBUG_BLEACHING_SCALAR` 1.5) should not exist.
   - CoralBlox is suspected to be too sensitive to heat stress: there should be little mortality below 6–8 DHW.
   - Calibrate with existing adaptation and tolerance terms in ADRIA/CoralBlox. No new parameters are to be added to those codebases.
   - A post-2016/17 reduction in chronic background mortality needs evidence from the LTMP disturbance tables (fewer unknown-cause and disease events).
2. **Chlorophyll index.**
   - Reef-based wet-season chlorophyll-a anomaly, plus the wet-season 90th percentile as an indicator of the maximum.
   - The lead's original COTS_Model survival function (`COTSModel_ChlorophyllModel.R`: logit-linear in log₂(chl), fitted to Fabricius et al. 2010) is too strict at low chlorophyll.
   - Years before the satellite record are to be predicted from the Niño 3.4 – chl-a anomaly relationship, with a GAM per reef (the lead's earlier workflow).
3. **Outbreak threshold:** 15 COTS/ha corresponds to about 0.22 COTS per manta tow.
4. **2026 predictions:** `sandbox/data/gbrPreds2026Adj_20262408.csv`. The other 2025–26 inputs use the labelled stand-ins.

## Evidence on post-2016/17 disturbance (2026-10-15)

**Source:** AIMS LTMP disturbance records (`BleachingMortality/data/aims_ltmp/reef_disturbance_complete.csv`; manta and photo-point, 1987–2026).

**Denominators:** LTMP reef-years surveyed (`reef_manta.csv`, `reef_photo_transect.csv`). North = sectors CL, CA, CG, CB.

| Survey | Region | Period | Reef-years | Unknown-cause events | Disease events |
| --- | --- | --- | --- | --- | --- |
| Manta | GBR | 1992–2015 | 1,589 | 8.9% | 0.1% |
| Manta | GBR | 2018–2025 | 711 | **5.1%** | 0.7% |
| Manta | North | 1992–2015 | 682 | 9.2% | 0.1% |
| Manta | North | 2018–2025 | 306 | **4.9%** | 1.6% |
| Photo | GBR | 1992–2015 | 1,037 | 2.0% | 2.7% |
| Photo | GBR | 2018–2025 | 485 | 1.4% | 2.7% |
| Photo | North | 1992–2015 | 441 | 2.3% | 5.0% |
| Photo | North | 2018–2025 | 195 | 1.0% | **6.7%** |

**Reading.**
- Unknown-cause losses roughly halved after 2017.
- Disease did not decline: it was flat GBR-wide and slightly higher in the north.
- This supports a partial reduction in chronic background mortality after 2017, not its removal. A test at half the pre-2017 μ is the upper bound the evidence allows.
- **Caveat:** post-2017 losses may be attributed to bleaching more readily than before.

## Preregistered coral heat-tolerance recalibration (2026-10-15)

**Hypothesis.** CoralBlox bleaching mortality is too sensitive. The calibrated mean tolerances (`dist_mean`, `:calib`) are 3.0–4.1 DHW for Acropora and Pocillopora and 5.0–6.6 for massives, with sd 2.9–5.5. That produces substantial fast-coral loss at 4–6 DHW, against observations of little mortality below 6–8 DHW.
- The extra 1.5× scalar compounds the problem.
- The model's slow post-2017 recovery and its 2020/2024 losses at Lizard and MacGillivray follow from this.

**Change.** No code change: existing scenario parameters only.
- **β removed:** `ADRIA_DEBUG_BLEACHING_SCALAR` unset (1.0).
- **Mean heat tolerance:** every `<coral_id>_dist_mean` scenario value is shifted by Δ ∈ {0, +1, +2, +3} DHW. The spreads (`dist_std`) are unchanged, and the tolerance ceiling moves with the mean (initial + `HEAT_UB`), as in ADRIA.
- **Heritability** of tolerance (existing factor, default 0.3, ADRIA range 0.25–0.5): {0.3, 0.5}.
- **Chronic background μ** (cv2 expression `1 − (1 − m_cyclone)^5 (1 − m_Mellin)(1 − μ)`): {0.06, 0.09, 0.12}.
- **× post-2017 multiplier** on μ: {1, 0.5}. The upper bound is set by the disturbance evidence above. It applies only to the 2018–2024 rows of the mortality cube.
- **Fixed:** bleaching below 4 DHW ×0.25 (the lead's evidence decision); growth untouched; κ 5.
- **COTS:** forced from observed CPUE / 0.015 (the reporting-year table), as in coral re-fit v2, so the screen isolates coral.
- **Arms:** `therm_b100_d<Δ>_h<h>_m<μ>_p<post>`, 4 × 2 × 3 × 2 = 48.
- **Reference:** `therm_b150_d000_h030_m120_p100` must replay coral re-fit v2 `refit2_s025_mu120_rep` (run 20) exactly.

**Read-out.**
- coral `J` 1992–2024 (photo fast and slow, manta total; as before) and `J` 1992–2016;
- **post-2017 recovery error:** RMSE 2018–2024 against manta HC and photo fast coral on Lizard, MacGillivray and North Direction;
- the model's 2024 change at Lizard and MacGillivray (observed: no loss);
- the 2016–17 losses (observed: large), so that heat tolerance does not erase the real mass bleaching.

**Selection rule.**
- Report the front of whole-period `J` against post-2017 recovery error.
- A candidate must keep `J` 1992–2016 within 10% of the reference, and reproduce at least half of the observed 2016–17 fast-coral loss at the three reefs.
- The lead chooses.
- The chosen coral is then re-tested with dynamic COTS (W) under score v3.

**Expected.**
- Δ +2 to +3 removes most of the 2020, 2022 and 2024 losses while keeping the 2016–17 losses (DHW 8–9).
- Heritability 0.5 adds acclimation after 2017.
- A lower post-2017 μ speeds the recovery.

**Failure criteria.** The reference does not replay, or no arm keeps the 2016–17 losses while removing the 2024 loss. In that case the tolerance structure itself is the issue, and it is reported.

**Rollback.** Runner-only.

## Plan update (2026-10-15)

- **Chlorophyll index.** Reef-based wet-season (December–April) chl-a mean anomaly and 90th percentile, 2003–2026, from DMS IMOS OC3 (`extract_sst_chla_features.py` method), for domain and upstream reefs. Back-cast to 1991–2002 with a GAM per reef on wet-season Niño 3.4 (BleachingMortality `extract_enso_indices.py`).
  - The larval-survival response would be preregistered separately, with a relaxed low-chlorophyll floor in place of the Fabricius logit's near-zero survival below about 0.4 µg/L.
- **Prediction-layer correction.**
  - **Threshold:** the outbreak indicator is diver-equivalent density ≥ 15 COTS/ha (lead), i.e. cull catch ≥ (15 / 7.9)^(1/0.55) ≈ 3.2 COTS per diver-hour under the SALAD anchor.
  - **Residuals:** at reef-years with SALAD or cull data, 2018–2026, the residual is (indicator − predicted probability). Residuals are interpolated by inverse distance to all source reefs, and the corrected probability is clamped to [0, 1].
  - **2026:** uses `gbrPreds2026Adj_20262408.csv`.
- **Domain v0.4 (1991–2026).**
  - DHW 2025–26 from DMS (validated);
  - cyclones 2025–26 zero;
  - Mellin background 2025–26 at its 2015–2024 mean (stand-in);
  - outbreak predictions 2026 from the lead's file;
  - v0.3 untouched.


## Coral heat-tolerance recalibration results (2026-10-15)

**Runs.**
- Run `20261015T_01_coral_thermal_1991_v1`: 49 arms, no failures. The reference replays coral re-fit v2 exactly (0.0). Audit: `coral_thermal_audit.json`.
- Diagnostics:
  - run 02: bleaching off;
  - run 03: failed, because `ADRIA_LOG_DHW_TOLS` is reset by `ADRIA.setup()` (FAILED.txt kept);
  - run 04: tolerance log via the new opt-in runner output `OWEN_1991_TOLERANCE_LOG`.

**Results** (forced COTS; μ 0.12 unless stated; Lizard + MacGillivray mean change in manta-equivalent HC from 2023 to 2024; observed +0.063):

| Arm | `J` all | `J` 1992–2016 | Post-2017 recovery error | Change 2023→2024 | Model HC 2024 (observed 0.42) |
| --- | --- | --- | --- | --- | --- |
| β 1.5 (reference) | 0.259 | 0.235 | 0.136 | −0.094 | 0.11 |
| β 1, Δ 0 | 0.258 | 0.236 | 0.143 | −0.071 | 0.18 |
| β 1, Δ +3 | 0.260 | 0.238 | 0.145 | −0.067 | 0.18 |
| β 1, μ 0.06 | 0.311 | 0.309 | **0.128** | −0.080 | 0.20 |
| β 1, μ 0.12, post-2017 μ halved | 0.261 | 0.236 | 0.151 | −0.080 | 0.21 |
| bleaching off (β 0) | – | – | – | **+0.001** (Lizard) | Lizard 0.23, MacGillivray 0.32 |

**Findings.**
1. **The 2024 loss is entirely bleaching.** No cyclone is recorded for 2016–2024 at these reefs (RRN wave-hours; Jasper is absent), and Mellin background is about 0. With bleaching off, Lizard's fast coral holds at 0.31 in 2024; with bleaching it falls from 0.27 to 0.17.
2. **Removing β (1.5 → 1) halves nothing.** Coral `J` is unchanged (0.259 → 0.258); the 2024 loss is reduced by a quarter.
3. **The mean-tolerance shift (Δ) and heritability are washed out by ADRIA's tolerance dynamics.**
   - Run 04: Lizard Acropora mean tolerance rises from 3.9 (1992) to 6.9 DHW (2016) at Δ 0, and from 6.9 to 7.4 at Δ +3. The +3 shift shrinks to about +0.3–0.5 before the 2016–17 events.
   - The cause: after any mortality, `bleaching_mortality!` resets each class mean to `truncated_normal_mean(μ, σ, HEAT_LB, ceiling)`, which lifts every mean towards a common attractor each heat year. The code comment says the truncation is at the current DHW, but the call uses `HEAT_LB`. That is reported, not changed (core ADRIA).
   - **So the effective fast-coral tolerance in the 2020s is already about 6–7 DHW.** The lead's "little mortality until 6–8 DHW" is close to the evolved mean. What kills about 40% of fast coral at DHW 6 is the tolerance spread (`dist_std` about 2.9–3.5 DHW), not the mean.
4. **Halving post-2017 μ helps manta HC** (0.144 → 0.133) **but worsens photo fast coral** (0.141 → 0.169): a mixed result, in line with the disturbance evidence (unknown-cause events halved, disease not).
5. **Recovery is too slow even without bleaching.** With bleaching off, modelled 2024 HC is 0.23 at Lizard and 0.32 at MacGillivray, against manta 0.21 and 0.47 (North Direction 0.58). Growth or reef-scale cover representation also limits the recovery at MacGillivray and North Direction.

**Existing levers that remain** (no new parameters):
- the tolerance spread (`<coral_id>_dist_std` scenario values);
- the existing opt-in switches `ADRIA_BLEACHING_LOW_DHW_THRESHOLD` (currently 4; mortality ×0.25 below it) and `ADRIA_BLEACHING_DHW_FLOOR` (minimum susceptibility).

**Decision.** Nothing is adopted. The options go to the lead.

## Chlorophyll index results (2026-10-15)

**Satellite record.** `extract_chla_wet_season.py`: MODIS-Aqua OC3 wet seasons for 2003–2024 extracted for the 113 domain reefs (`sandbox/data/derived/chla_wet_season_lizard/`). The IMOS product ends on 2024-12-25 (Aqua retired), so 2025 and 2026 are missing.
- **High anomaly seasons** (domain median log anomaly): 2011 +0.18, 2019 +0.13, 2024 +0.27. These are flood years (2010–11 La Niña; February 2019; Jasper December 2023).
- Absolute OC3 values (about 0.8–1.2 mg m⁻³) are inflated over shallow reef water, so only reef-relative anomalies are used.

**Niño 3.4 back-cast** (`build_chla_enso_backcast.R`, per-reef GAMs): **no skill.** Leave-one-season-out CV R² is −0.43 to −0.52 for the wet-season and preceding-spring Niño 3.4, for both mean and p90. The output is flagged NOT_FOR_USE. The flood-driven anomalies are not ENSO-predictable here; river discharge is a more plausible predictor.

**Other DMS sources.** The eReefs biogeochemistry stores (gridded baseline, and the reef-level `cots-ereef-biogeochemistry`) both end in March 2019.

**Outstanding:** a 2025–26 chlorophyll source and a pre-2003 predictor. These go to the lead.


## Lead's decisions (2026-10-15, tenth round)

1. **Bleaching sensitivity:** test a narrower tolerance spread (existing `<coral_id>_dist_std`) together with a raised reduced-mortality threshold (existing `ADRIA_BLEACHING_LOW_DHW_THRESHOLD`, 4 → 6 DHW, ×0.25 below). Judge on the 2016–17 losses, which must remain, and the 2020/2024 non-losses.
2. **Recovery:** check the manta and model cover comparison and fix bleaching first. Revisit growth only if a gap remains.
3. **Chlorophyll:** use NOAA VIIRS, calibrated to MODIS, for 2025–26, and test river discharge as the pre-2003 predictor.
4. **Builds:** the corrected outbreak layer and the 2026 domain are built next (not waiting on the lead).

## Manta against model cover (2026-10-15)

Comparison of the forced-COTS μ 0.12 arms, with bleaching on and off (runs 20261015T_01 and _02), against the manta HC median.

**Lizard.** Manta matches the site-area fraction (`total_coral_site_fraction`; 2007: 0.22 vs manta 0.14). The habitable-area fraction is about 2× too high (0.36), so the current comparator is right.

**MacGillivray and North Direction.** Manta in 2021–24 (MacGillivray 0.19–0.47; North Direction 0.50–0.58) exceeds the model's maximum cover in any year, even with bleaching off (about 0.32–0.37 of habitable area).
- The model's cover ceiling under μ 0.12 is below these unprecedented values.
- North Direction's manta jump from 0.09 (2019) to 0.50 (2021) is faster than coral growth allows, which suggests sampling noise on a small reef (few tows).
- **Conclusion:** comparator scale is not the cause. Bleaching (question 1) and the chronic-background equilibrium are the levers; growth is not revisited yet.

## Preregistered tolerance-spread × low-DHW threshold screen (2026-10-15)

**Hypothesis.** ADRIA ratchets fast-coral mean tolerance to about 6–7 DHW by 2016, so the large losses at DHW 4–6 in 2020 and 2024 come from the tolerance spread (sd about 2.9–3.5) and from full mortality above 4 DHW.
- A narrower spread steepens the dose–response around the evolved mean.
- A threshold of 6 limits mortality below 6 DHW to ×0.25.
- The 2016–17 events (DHW 8–9) remain lethal.

**Change.** No code change; existing scenario values and switches only.
- **Experiment:** `coral_spread`, forced COTS (reporting-year table, `q` 0.015).
- **Fixed:** β unset (1.0); mean tolerance unshifted; heritability 0.3; growth untouched; κ 5.
- **Factors:**
  - every `<coral_id>_dist_std` × s ∈ {1, 0.75, 0.5};
  - low-DHW threshold t ∈ {4, 6} DHW (×0.25 below);
  - μ ∈ {0.12, 0.09}.
- **Arms:** `spread_s<s>_t<t>_m<μ>` (12), plus bleaching-off counterfactuals `spread_off_t40_m<μ>` (2) for attribution.
- **Reference:** `spread_s100_t40_m120` must replay run 20261015T_01 `therm_b100_d000_h030_m120_p100` exactly. `spread_off_t40_m120` must replay run 02 `therm_b000_d000_h030_m120_p100`.

**Read-out**, per arm (Lizard, MacGillivray and North Direction):
- the bleaching-attributable relative loss, `1 − cover / cover_bleaching_off` at the same μ, for 2016–17 (expected large, about 50% of fast coral in the northern GBR 2016 event) and for 2020 and 2024 (observed small);
- coral `J` 1992–2024 and 1992–2016;
- the post-2017 recovery error.

**Selection rule.** Among arms that keep at least half of the reference's 2016–17 bleaching-attributable fast-coral loss, take the lowest 2020+2024 attributable loss; ties go to coral `J`. The lead decides.

**Failure.** A reference does not replay, or an arm fails.


## Tolerance-spread × low-DHW threshold results (2026-10-15)

**Run.** `20261015T_05_coral_spread_1991_v1`: 14 arms, no failures. Both references replay exactly (0.0). Audit: `coral_spread_audit.json`.
- The audit's final console print failed on a display-only column (fixed after the outputs were written). The recorded `script_sha256` is the pre-fix version; the outputs are unaffected.

**Bleaching-attributable fast-coral loss** (Lizard, MacGillivray and North Direction mean, against bleaching off at the same μ; μ 0.12):

| Spread × | Threshold | 2016 | 2017 | 2020 | 2022 | 2024 | `J` all | Post-2017 error |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 (reference) | 4 | 72% | 55% | 59% | 22% | 45% | 0.258 | 0.143 |
| **1** | **6** | **72%** | **55%** | **48%** | **13%** | **16%** | **0.262** | 0.150 |
| 0.75 | 4 | 93% | 71% | 79% | 43% | 69% | 0.259 | 0.137 |
| 0.75 | 6 | 93% | 71% | 65% | 22% | 25% | 0.260 | 0.143 |
| 0.5 | 4 | 100% | 78% | 94% | 80% | 97% | 0.279 | 0.151 |
| 0.5 | 6 | 100% | 78% | 72% | 33% | 36% | 0.269 | 0.154 |

**Findings.**
1. **A narrower tolerance spread increases mortality.** ADRIA's tolerance ratchet (class means reset to the mean truncated at `HEAT_LB` after any mortality) lifts the means less when the spread is small, so the evolved tolerance stays near the DHW of later events.
2. **Raising the existing low-DHW threshold from 4 to 6 DHW** (mortality ×0.25 below it) keeps the 2016–17 mass-bleaching losses unchanged, and cuts the attributable fast-coral loss in 2024 from 45% to 16% and in 2022 from 22% to 13%.
   - 2020 stays at 48%, since its DHW (5.5–6.6) is partly above 6.
   - Coral `J` changes 0.258 → 0.262; the post-2017 error moves 0.143 → 0.150.
3. **μ 0.09 against 0.12** changes the losses little but worsens `J` (as in run 45).

**Selection (preregistered rule):** `spread_s100_t60_m120`, i.e. tolerance spread unchanged, low-DHW threshold 6 DHW (×0.25 below), μ 0.12, no bleaching multiplier. The lead decides; the next step is to re-test it with dynamic COTS under score v3.

## Corrected outbreak-prediction layer (2026-10-15)

**Builder:** `domain_building/build_corrected_outbreak_layer.py`. Output: git-ignored `sandbox/data/derived/outbreak_layer_corrected_2018_2026/`.

**Inputs.** The lead's predictions for 1991–2025 plus the 2026 file. Observed outbreak indicators for reporting years 2018–2026:
- SALAD density ≥ 15 COTS/ha;
- otherwise cull catch ≥ 3.235 COTS per diver-hour (the 15/ha equivalent under the SALAD anchor), with ≥ 2 h of bottom time.

**Method.** Residuals (indicator − predicted probability) are spread by inverse distance within 50 km and shrunk toward zero (a single observation at 50 km contributes half). The result is clamped to [0, 1]. Years 1991–2017 are unchanged.

**Result.**
- 66–151 observed reef-years per year, of which 12–29 per year are outbreaks.
- Northern GBR mean outbreak probability (12–17.5° S), original → corrected:

| Year | Original | Corrected |
| --- | --- | --- |
| 2018 | 0.041 | 0.051 |
| 2019 | 0.018 | 0.012 |
| 2020 | 0.025 | 0.047 |
| 2021 | 0.029 | 0.080 |
| 2022 | 0.024 | 0.033 |
| 2023 | 0.018 | 0.040 |
| 2024 | 0.031 | 0.086 |
| 2025 | 0.038 | 0.138 |
| 2026 | 0.031 | 0.199 |

  For comparison, 2012–14 was 0.25.

**Caveat.** Culling is targeted at outbreaks, so cull data are present mainly where COTS were found. Unculled reefs keep the original probability; they are not treated as evidence of no outbreak.


## Chlorophyll 2025–26 (VIIRS) (2026-10-15)

**Builder.** `domain_building/extract_chla_viirs_wet_season.py`: NOAA CoastWatch `noaacwNPPVIIRSSQchlaDaily` (S-NPP VIIRS Science Quality, 4 km). Same reef pixel and season rules as MODIS. Output: `sandbox/data/derived/chla_index_lizard_2003_2026/`.

**Calibration against MODIS** (overlap 2013–2024, 1,152 reef-seasons):
- reef-season correlation of the log anomalies: r 0.56 (mean) and 0.54 (p90);
- domain-median r: 0.84 and 0.86;
- VIIRS anomalies are scaled to MODIS by no-intercept slopes of 0.43 (mean) and 0.40 (p90).

**Combined index:** MODIS 2003–2024, plus calibrated VIIRS for 2025–2026. Domain medians: 2025 +0.063, 2026 +0.012.

**Pre-2003:** no valid predictor yet (the Niño 3.4 GAMs have no skill). River discharge has not yet been tested.

## Domain v0.4 (1985–2026) and the corrected COTS dataset (2026-10-15)

**`domain_building/build_lizard_domain_v4.py`** (output `sandbox/data/Lizard_Historical_v0.4`, git-ignored): v0.3 copied, plus
- DHW 1985–2026 (2025–26 from the validated DMS extraction);
- cyclones 2025–26 zero;
- `domain_data_version` 0.4.0.

**`extend_mellin_background_2026.py`:** Mellin background extended to 2026 (stand-in: 2015–2024 site mean).

**`build_lizard_cots_connectivity_owen.py`** gains `--predictions` (one or more files) and `--dataset-name`; the defaults are unchanged. v0.4 holds:
- `owen_global_2018_2023`: the lead's predictions 1991–2025 plus the 2026 file. The boundary is identical to v0.3 for 1991–2025 (maximum difference 0.0);
- `owen_global_2018_2023_corrected`: the corrected layer.

  Domain-summed external outbreak boundary (original → corrected):

  | Year | Original | Corrected |
  | --- | --- | --- |
  | 2018 | 5.1 | 5.8 |
  | 2019 | 2.8 | 4.3 |
  | 2020 | 1.8 | 4.2 |
  | 2021 | 1.7 | 5.0 |
  | 2022 | 1.9 | 2.7 |
  | 2023 | 1.7 | 3.2 |
  | 2024 | 2.8 | 16.9 |
  | 2025 | 3.3 | 37.3 |
  | 2026 | 3.2 | 39.0 |

  For comparison, 2016 was 22.5.

**Runner switches.** `OWEN_1991_DOMAIN_VERSION` (`v0.3` default, or `v0.4`; v0.4 is for dynamic-COTS experiments only) and `OWEN_1991_COTS_DATASET`. The re-score tool accepts 1991–2026 trajectories; manta scoring stays within 1992–2024, and diver targets extend to 2026.

**Smoke and replay checks** (W at seed pair 1, old cv2 coral):
- **Run 06** (v0.4, original dataset) replays v0.3 run 37 for 1991–2024: maximum relative difference about 4e-15 in trajectories, flows and mechanisms.
- **Run 07** (v0.4, corrected dataset): the fifth outbreak appears but late and booming.

| Reef | Run 07 adults/ha | Targets |
| --- | --- | --- |
| Lizard | 4.8 in 2024, 16.8 in 2026 | 3.4 → 9 |
| Eyrie | 8.7 in 2024, 14.7 in 2026 | 17–45 from 2022 |
| North Direction | 9.5 in 2026 | 11–13 in 2024–25 |
| MacGillivray | 22 in 2026 | 2.5–4.3 |

2020–2022 stays below 1/ha, against targets of about 3.5/ha.


## Lead's decisions (2026-10-15, eleventh round)

1. **Coral for the 2020s runs:** no bleaching multiplier (β 1); tolerance spread unchanged; bleaching mortality ×0.25 below 6 DHW (existing `ADRIA_BLEACHING_LOW_DHW_THRESHOLD`); chronic μ 0.12/yr.
   - The working configuration's arm prefix `wc_b150` sets β 1.5 through the runner. The new-coral arms use `wc_b100`.
2. **Chlorophyll larval survival.** The proposed form and a new opt-in ADRIA adapter switch are accepted. A plot was requested and made (`figures/larval_survival_curve_proposal.png`, `plot_larval_survival_curve.py`).
3. **Experiment.** Run the new coral × supply dataset × chlorophyll survival experiment to 2026.
   - The lead accepts that the 2020–22 build-up may be partly missed, but relates it to the February 2019 flooding. The 2019 wet season has the second-highest chl-a anomaly after 2011 (domain median ×1.4–1.5 at b = 3).
4. **Pre-2003:** the chlorophyll index is neutral (×1) until a river-discharge predictor is tested.

## Preregistered chlorophyll larval-survival switch and 2020s experiment (2026-10-15)

**Hypothesis (initiation region; lead).** Wet seasons with high chlorophyll (floods, upwelling) raise survival of larvae produced within the domain, so successive good years build local recruitment. The 2019 floods would then lift recruits in 2019 and adults by about 2021–22. The 2024 Jasper floods (the largest anomaly) feed the 2025–26 build-up.

**Mechanism** (opt-in; units dimensionless).
- **Multiplier:** `m(reef, season) = min(5, (chl / reef normal)^b) = min(5, exp(b · log_anom))`, using the reef wet-season chl-a index from `chla_index_lizard_2003_2026` (MODIS 2003–2024, VIIRS-calibrated 2025–2026). It is 1 before 2003 and where missing.
- **Where:** applied at the source to the pre-pelagic larval production of each domain site (`last_local_fecundity`, separated recruitment) in model year = season year (Dec(y−1)–Apr(y) spawning), after the COTS step and before dispersal. So it scales both local retention and export to other domain reefs. External boundary supply is unchanged.
- **ADRIA:** env `COTS_LARVAL_SURVIVAL_CSV` (columns `site_id, year, multiplier`). Unset gives the established path. Applied multipliers are logged as `cots_larval_survival_log`.
- **Tables:** `build_larval_survival_index.py` builds them (git-ignored derived data). The variants are wet-season mean or p90, × b ∈ {1.5, 3}.

**Tests.**
- the loader: coverage, validation, and years outside the table → 1;
- the applier: scales production only;
- a zero-effect table is the identity.

**Experiment** `chl_survival` (domain v0.4, to 2026; working configuration W at seed pair 1; one run per supply dataset {`owen_global_2018_2023` (the lead's predictions plus 2026), `owen_global_2018_2023_corrected`}).
- **Replay arm:** `wc_b150_..._t40_chloff_cv2_pnd`, the old coral (W), must replay smoke runs 06 / 07 exactly.
- **New coral (β 1, threshold 6):** chl off, mean b 1.5, mean b 3, p90 b 1.5, p90 b 3.

**Read-out.**
- score v3, with diver targets 2020–2026;
- adults/ha against targets per reef and year, highlighting 2020–22 (the post-2019-flood signal) and 2024–26;
- recruits in 2019–20;
- coral 2020–26 (boom-bust check);
- v2 and G2′ for 1992–2016, to confirm the earlier waves are not degraded.

**Expected.**
- The new coral removes most of the 2024 coral crash, so the boom-bust is damped.
- Chl survival raises 2019–20 recruitment, giving adults of about 1–3/ha by 2021–22, and amplifies 2025–26.
- The corrected supply adds the 2024–26 build-up.
- b = 3 could overshoot the 2011–13 wave (2011 ×1.7–1.9), so the second-wave height is checked.

**Failure criteria.** A replay fails, an arm fails, or the 1992–2016 v2 loss worsens by more than 1.0 in a selected arm.

**Rollback.** The env is unset by default.


## Chlorophyll larval-survival × new coral × supply results (2026-10-15)

**Implementation.**
- **ADRIA:** `load_cots_larval_survival` and `apply_cots_larval_survival!` (`cots.jl`); env `COTS_LARVAL_SURVIVAL_CSV`; log `cots_larval_survival_log`. Adapter test `test_larval_survival_index.jl`: 10/10.
- **Tables:** `build_larval_survival_index.py` → git-ignored `larval_survival_index_lizard/`. Domain-median multipliers at b 3: 2011 ×1.5–1.8, 2019 ×1.3, 2024 ×1.8–2.2.

**Runs.**
- Runs 20261015T_08 (the lead's predictions + 2026) and _09 (corrected layer): domain v0.4, 6 arms each, no failures.
- Re-scores: 10 and 12. G2′: 11 and 13. Runs 10–13 were renamed for clarity (RUN_RENAMES.csv).
- The replay arms (old coral) reproduce smoke runs 06 and 07 exactly (0.0).
- Audit: `runs/20261015T_09_.../chl_survival_audit.json` and `chl_survival.png`.

**Score v3** (seed pair 1; diver targets 2020–2026; lower is better):

| Supply | Arm | v3 | Manta part (1992–2021) | Diver error (sum of 4 reefs) | v2 1992–2024 | G2′ at `A½ 30` |
| --- | --- | --- | --- | --- | --- | --- |
| original | old coral (W) | 15.61 | 6.51 | 6.07 | 6.44 | pass |
| original | **new coral**, chl off | **11.34** | 6.60 | **3.16** | 6.54 | pass |
| original | new coral, p90 b 3 | 11.21 | 6.42 | 3.19 | 6.37 | pass |
| corrected | old coral (W) | 12.63 | 6.51 | 4.08 | 6.45 | pass |
| corrected | new coral, chl off | 11.34 | 6.60 | 3.16 | 6.59 | pass |
| corrected | **new coral, p90 b 3** | **10.99** | 6.43 | 3.04 | 6.41 | pass |

**Adults/ha 2020–2026** (new coral, chl off) against targets:

| Reef | Original supply | Corrected supply | Target |
| --- | --- | --- | --- |
| Lizard | 5.1 → 3.9 (2024) → 6.1 | 4.9 → 12.7 (2024) → 36 | 3.4 → 8.6 → 5.0 |
| North Direction | 2.3 → 7.0 → 10.5 | 2.3 → 16.6 → 55 | 3.6 → 12.8 → 7.2 |
| MacGillivray | 2.9 → 3.9 → 5.4 | 2.7 → 8.9 → 34 | 2.5 → 4.3 |
| Eyrie | 0.8 → 2.3 → 2.2 | 0.7 → 13.6 → 22 | 1.9 → 45 → 34 |

**Findings.**
1. **The new coral is the main improvement.**
   - It removes the 2020–22 starvation collapse: adults persist at 2–6/ha, within the diver range.
   - Fast coral in 2026 is about 0.17–0.26, against 0.04–0.09 with the old coral.
   - Diver error halves (6.07 → 3.16). The 1992–2016 fit and G2′ are kept (v2 6.44 → 6.54).
   - The lead's 2020–22 concern is largely met through adult persistence rather than new recruitment.
2. **The chlorophyll survival effect is small** (v3 −0.1 to −0.35).
   - 2019–20 recruitment rises only 2–4% on Lizard, MacGillivray and North Direction, because external supply (Ext×5) dominates settlement. Local production in 2019 is small at the 2019 population levels.
   - The 2024 Jasper index raises 2026 adults (Lizard 6.1 → 8.9 at p90 b 3).
   - b 3 slightly raises the 2010s second wave (Lizard 15.5 → 19.5; North Direction 57 → 68).
3. **The corrected layer with Ext×5 overshoots in 2025–26** (Lizard 36–52 vs 5; North Direction 55–78 vs 7; MacGillivray 34–53 vs 2.5–4.3), while fixing Eyrie (13.6–22 vs 12–45, against 2.2 with the original).
   - Ext×5 was set with the original manta-based layer (0.02–0.04 probability in the 2020s). With the corrected layer it multiplies an already larger 2024–26 supply.
4. **Implication.** In the 2020s the external multiplier and the outbreak layer are confounded. A lower multiplier would let local recruitment, and so the chlorophyll mechanism, carry more of the build-up. That is closer to the lead's initiation-region hypothesis.

**Decision.** Nothing is promoted. The lead decides on adopting the new coral and on the next factorial (external multiplier × supply × chlorophyll).


## Lead's decisions (2026-10-15, twelfth round)

1. **Simulation candidates:** adopt the new coral (β 1; ×0.25 below 6 DHW) for all simulation candidates (W, C3a, C6), re-checked over the 5 seed pairs.
2. **External multiplier test.** Test a lower external-supply multiplier.
   - **Priority:** local build-up. If chlorophyll-driven survival can be shown to help start the 4th (2010s) and 5th (2020s) waves, where satellite data exist, then a 1990s-only boost may be considered, but only if chlorophyll cannot be back-cast.
3. **Outbreak layer:** choose between the original and corrected layers after the multiplier test.

## Preregistered external-multiplier × supply × chlorophyll factorial (2026-10-15)

**Hypothesis.** Ext×5 swamps local recruitment, so the chlorophyll larval-survival effect cannot act. With lower external supply:
- local production and retention carry more of each build-up;
- high-chl seasons (2011, 2019, 2024) advance or amplify the 4th and 5th waves;
- the corrected layer no longer overshoots in 2025–26.

The 1990s wave (chl neutral before 2003) is expected to weaken as Ext falls.

**Design.** Experiment `ext_chl`: domain v0.4, 1991–2026, W at seed pair 1, new coral (β 1, ×0.25 below 6 DHW, μ 0.12).
- external multiplier ∈ {1, 2.5, 5} × chlorophyll survival ∈ {off, p90 b 3} = 6 arms per run;
- one run per supply dataset (the lead's predictions + 2026; corrected layer);
- arms `wc_b100_q20_eb_lin_obs015_s1_t60_e<Ext×10>_<chl>_cv2_pnd`;
- **controls:** e50 chl off and e50 p90 b 3 replay the matching run 08 / 09 arms exactly.

**Read-out.**
- score v3 (manta part and diver error) and v2 1992–2024; G2′ at `A½ 30` and 10;
- **wave starts:** per reef, the first year adults reach 3/ha and the peak adults/year in 1993–2001 (wave 3, chl neutral), 2008–2015 (wave 4) and 2019–2026 (wave 5);
- **the chlorophyll effect** = p90 b 3 minus off at the same Ext and supply, in settled recruits (2011–12, 2019–20, 2024–25) and in wave-4 and wave-5 adults and onset years;
- **the local share of settlement** (local + internal immigration over total) in 2010–2014 and 2019–2026.

**Expected.**
- Lower Ext raises the local share, and so the chlorophyll effect, on wave starts.
- Wave 3 weakens at Ext 1.
- The corrected layer at Ext 1–2.5 builds 2024–26 without the ×5 overshoot.

**Criteria for the lead.**
- the chlorophyll effect on wave-4 and wave-5 onset or size;
- v3;
- no loss of waves 3–4 under G2′.

**Failure.** A replay fails, or an arm fails.

**Rollback.** Runner-only (existing switches).


## External-multiplier × supply × chlorophyll results (2026-10-15)

**Runs.**
- Runs 20261015T_14 (original layer) and _15 (corrected layer): 6 arms each, no failures.
- Re-scores: 16 and 18. G2′: 17 and 19. Audit: `runs/20261015T_15_.../ext_chl_audit.json` and `ext_chl.png`.
- The e50 arms replay runs 08 / 09 exactly (0.0).

**Score v3** (seed pair 1):

| Ext | Original, chl off / on | Corrected, chl off / on | Manta part (off) | Matched v2 peaks (off) | G2′ at `A½ 30` |
| --- | --- | --- | --- | --- | --- |
| ×1 | 28.6 / 28.5 | 28.0 / 28.3 | 20.9 | 1–2 of 7 | fail |
| ×2.5 | 16.9 / 16.7 | 16.0 / 16.1 | 10.0 | 7 | fail |
| ×5 | 11.3 / 11.2 | 11.3 / **11.0** | 6.6 | 7 | pass |

**Waves** (Lizard; onset = first year ≥ 3 adults/ha in the window):
- **Wave 3 (1993–2001):** peak 30/ha in 1997 at ×5, but 32 in 2000 at ×2.5 and 12 in 2001 at ×1. Lower Ext delays and shrinks the 1990s wave, which has no chl support.
- **Wave 4 (2008–2015):** peak 15.5 at ×5 against 7.8 at ×1. Adults are already ≥ 3/ha by 2008 in every arm, so onset is not discriminated.

**Local share of settlement**, i.e. 1 − external/settled, increases as Ext falls:
- Lizard 2010–14: 14% at ×5, 22% at ×1;
- Lizard 2019–26 (original layer): 33% at ×5, 78% at ×1.

**Chlorophyll effect** (p90 b 3, on − off):
- **Recruits:**
  - 2011–12: +10.7% at ×5, +15.7% at ×1 (Lizard); +18–27% at MacGillivray;
  - 2019–20: +3–7% at ×5, +19–20% at ×1;
  - 2024–25: +21–35% at ×5, +42–96% at ×1.
- **Wave-4 peak:** +14–27% on Lizard, MacGillivray and North Direction at every Ext.
- **Wave-5 peak:** +21–75%.
- **Onset:** never earlier. In this configuration chl amplifies waves rather than advancing them.

**Findings.**
1. **Ext×5 cannot be reduced without losing the 1990s and 2010s manta fit** (v3 ×5 11.3 → ×2.5 about 16.5 → ×1 about 28). The external multiplier is doing work that local dynamics alone do not.
2. **Chl survival does act on the 4th and 5th waves** (+15–27% recruits in 2011–12; larger effects in 2019–20 and 2024–25 as the local share rises), but on size, not timing.
3. **The 1990s wave has no chl index** (neutral before 2003), so it cannot benefit at lower Ext. Per the lead (twelfth round), a 1990s-only boost is considered only if chl cannot be back-cast.

## In-water chlorophyll 1992–2009 (lead, 2026-10-15)

**Source.** The lead found `sandbox/data/GBR_chlorophyll_monitoring_1992-2009.csv` (AIMS Long Term Chlorophyll Monitoring Program; 19,174 samples; CHL_A µg/L).
- Lizard region (within 0.6° of Lizard Island): 3,154 samples. The Lizard Island transect (Codhole, MacGillivray, North Reef, Eagle Islet, Near Martin Reef, …) was sampled monthly 1993–1999; Cooktown–Osprey 1997–2008.
- **Wet-season (Dec–Apr) regional mean:**

| Season | Mean (µg/L) |
| --- | --- |
| 1993 | 0.25 |
| 1994 | 0.28 |
| 1995 | 0.28 |
| 1996 | 0.34 |
| 1997 | 0.24 |
| 1998 | 0.31 |
| 1999 | 0.33 |
| 2000 | 0.18 |
| 2001 | 0.19 |
| 2002 | 0.25 |
| 2003 | 0.14 |
| 2004 | 0.24 |
| 2005 | 0.24 |
| 2006 | 0.21 |
| 2007 | 0.21 |

  The 1990s were higher than the 2000s.

**Overlap with MODIS** (seasons 2003–2008):
- satellite OC3 is about 3–4× the in-water level (bottom and coastal contamination);
- the season-level interannual correlation is r −0.19 over 6 seasons (2008 has only 4 samples).

A season-level bias correction would therefore be fragile. Same-day matchups (in-water sample against the satellite pixel, 2002–2008) give hundreds of pairs for the calibration.

**Lead's direction.** Make a full 1992–2026 series, with the in-water data as the reference for bias correction, because surface satellite chl-a underestimates what larvae are exposed to.


## Lead's decisions (2026-10-15, thirteenth round)

1. **Satellite-to-in-water correction by same-day matchups** (in-water sample against the satellite pixel, 2002–2008).
2. **Test both survival forms:** the relative anomaly, and an absolute curve on the in-water scale (relaxed Fabricius).
3. **One Lizard-region in-water value per season for 1993–2002,** with each reef's offset carried back from the satellite years.
4. **Re-test the external multiplier** once the full 1993–2026 series exists.

## Preregistered full chlorophyll series (in-water scale) and Ext re-test (2026-10-15)

**Matchups** (`domain_building/build_chla_inwater_matchups.py`).
- **In-water:** AIMS LTCMP samples north of 17.5° S, July 2002 – January 2009. Each station visit gives the mean CHL_A over all sampled depths (the larval-exposure proxy; a near-surface-only sensitivity uses ≤ 10 m).
- **Satellite:** the same-day MODIS-Aqua OC3 at the nearest valid pixel within 5 km.
- **Fit:** `log10(in-water) = a + b · log10(OC3)` by OLS, with n, r, residual sd, and per-transect checks.

**Series** (`build_chla_index_full.py`; wet season Dec(y−1)–Apr(y); per reef):
- **2003–2024:** MODIS season mean and p90 converted with the matchup fit.
- **2025–2026:** VIIRS season values mapped to MODIS by a per-reef log-log fit over 2013–2024, then converted.
- **1993–2002:** the Lizard-region in-water wet-season mean and p90 (stations within 0.6° of Lizard Island; seasons with ≥ 16 samples), times each reef's mean offset from the regional value over 2003–2026 (converted satellite).
- **Anomalies:** relative to each reef's 1993–2026 geometric mean.

**Survival forms.**
1. **Anomaly:** `min(5, exp(3 · log_anom_p90))`, as before but now covering 1993–2026.
2. **Absolute, relaxed Fabricius:** `S(c) = f + (1 − f) · logistic(−0.48 + 2.02 · log2 c)` on the in-water-equivalent p90 c. The multiplier is `min(5, S(c) / S(c_ref))`, with `c_ref` the domain 1993–2026 median.
   - The floor f relaxes the near-zero survival at low chl. A plot of the alternatives goes to the lead; f = 0.1 is used in the run unless the lead changes it.

**Experiment** `ext_chl2` (v0.4, W at seed pair 1, new coral; both supply layers):
- Ext {1, 2.5, 5} × chl {off, anomaly (full), absolute (f 0.1)};
- **controls:** the Ext and chl-off arms replay runs 14 / 15.

**Read-out.** As in the Ext × chl factorial, plus the wave-3 onset and peak with chl now active in the 1990s.

**The question:** does chl-driven local recruitment let a lower Ext keep waves 3–4, i.e. is the ×5 boost replaceable?


## Full chlorophyll series: build notes and correction (2026-10-15)

**Matchups** (`build_chla_inwater_matchups.py`): 834 station visits, of which 323 have a valid same-day MODIS pixel. All-depth OLS: `log10(iw) = −0.604 + 0.330 · log10(OC3)`, r 0.52, residual sd 0.26 (log10). In-water is 0.42× OC3 at the median.
- By transect: Cairns b 0.49 (r 0.80); Cooktown–Osprey b 0.16 (r 0.28, n 119); Port Douglas b 0.27.
- The near-surface sensitivity has only 5 pairs and is not usable.

**Correction to the preregistered conversion** (made before any model run). The p90 forms use the **RMA** conversion (b = 0.33 / 0.52 = 0.63), not OLS.
- Over 2003–2007, the in-water regional means match both conversions (in-water about 0.21; OLS 0.23; RMA 0.25–0.28).
- The in-water p90 (about 0.32) is matched by RMA (0.33–0.36) but compressed by OLS (0.26–0.28). Daily pixel noise attenuates the OLS slope, and the compressed OLS upper tail would make every satellite season look low against the in-water 1990s.
- The mean-anomaly variant keeps OLS.

**Series** (`build_chla_index_full.py`): MODIS-era rows are limited to ≤ 2024. MODIS ends 2024-12-25; its empty 2025–26 rows had shadowed the VIIRS rows and been mapped to the cap by a Python `min(5, NaN)`. That was fixed with a duplicate guard before any run, and the affected derived outputs were deleted and rebuilt.

**Tables** (`build_larval_survival_index_full.py`; `figures/larval_survival_full_series.png`): reef-median multipliers.

| Season | Anomaly, season mean (OLS) | Anomaly, season p90 (RMA) | Absolute relaxed, p90 (RMA; f 0.1, c_ref 0.361) |
| --- | --- | --- | --- |
| 1996 | 2.61 | 3.61 | 1.50 |
| 1998 | 2.05 | 2.27 | 1.26 |
| 1999 | 2.57 | 3.11 | 1.41 |
| 2000 | 0.37 | 0.21 | 0.83 |
| 2011 | 1.11 | 1.40 | 1.09 |
| 2019 | 1.06 | 1.16 | 1.04 |
| 2024 | 1.21 | 1.65 | 1.13 |

1991–1992 have no data (neutral, ×1).

**Experiment `ext_chl2`:** Ext {1, 2.5, 5} × chl {off, anomaly-mean OLS (`amols`), anomaly-p90 RMA (`aprma`), absolute RMA (`abrma`)} × both supply layers. The chl-off arms replay runs 14 / 15 by name.


## Full chlorophyll series × external multiplier results (2026-10-15)

**Runs.**
- Runs 20261015T_22 (original layer) and _23 (corrected layer): 12 arms each, no failures. Runs 20–21 failed at load (`ext_chl2` was missing from an experiment list; FAILED.txt kept).
- Re-scores: 24 and 26. G2′: 25 and 27. Audit: `runs/20261015T_23_.../ext_chl_audit.json` (generalised `audit_1991_ext_chl.py`).
- The chl-off arms replay runs 14 / 15 exactly (0.0).

**Score v3**, range over the chl variants:

| Ext | Original layer | Corrected layer | Best |
| --- | --- | --- | --- |
| ×1 | 27.2–28.6 | 26.7–28.0 | p90 anomaly (RMA) |
| ×2.5 | 16.5–17.4 | 15.8–16.7 | absolute (RMA) |
| ×5 | 11.17–11.34 | **10.98**–11.45 | p90 anomaly (RMA) |

- G2′ at `A½ 30` passes only at ×5.
- The 1990s manta part at ×1 improves with the p90 anomaly (20.9 → 19.5) but remains about 3× the ×5 value.

**Wave 3 (1990s; chl now from in-water data)**, Lizard / North Direction:
- ×1: onset 1998 / 1996 in every chl arm. The peak rises from 12 to 20/ha (Lizard) and from 23 to 42/ha (North Direction) with the p90 anomaly (+60% / +85%), but still peaks in 2001.
- ×2.5: onset 1996 / 1995; peak in 2000.
- ×5: onset 1994; peak in 1997 / 1998.

**Chl never advances an onset.** The high-chl seasons (1996 ×3.6, 1998 ×2.3, 1999 ×3.1) follow the observed onset (about 1994–95), and 1993–95 are near neutral (×0.8–1.85).

**Chl effects at ×5** (p90 anomaly, RMA):
- wave-4 peak +48% at Lizard, +23% at North Direction;
- 1996–97 recruits +16–22%;
- v3 −0.2 to −0.4.

**Findings.**
1. **The back-cast chlorophyll series amplifies but does not initiate** the 1990s wave. At lower Ext it raises the wave-3 peak strongly but cannot correct its timing, so it cannot replace Ext×5.
2. **At ×5, chl is a modest, consistent improvement,** amplifying waves 3–5 after high-chl seasons. The p90 anomaly with RMA is the best variant.
3. **At ×2.5, chl worsens v3** (+0.5): amplified waves overshoot the manta heights.
4. Per the lead's twelfth-round condition, the 1990s question is now: chl can be back-cast, but does not start that wave. A 1990s-only external boost, with a lower Ext afterwards, is the remaining option if the lead wants local dynamics to carry waves 4–5.


## Lead's decisions (2026-10-15, fourteenth round) and next steps

1. **Keep chl larval survival at Ext×5** (p90 anomaly, RMA conversion).
   - The lead doubts that the 1990s multipliers (1996 ×3.6, 1998 ×2.3, 1999 ×3.1) are really that much higher than the satellite era. It is most likely an artefact: there is no high-chl season in the in-water/satellite overlap (2003–2008), so the splice cannot calibrate the peaks.
   - **Bias-correct the 1990s segment:** bring the 1990s p90-anomaly peaks into line with the 2011 and 2024 satellite-era peaks, or vice versa.
2. **No 1990s-only boost test.** Ext×5 is kept consistently across all years.
3. **Carry both outbreak layers** (the lead's predictions + 2026, and the corrected layer) into the candidate set as a stated uncertainty.
4. **Run the candidate re-check over 5 seed pairs** with the new coral (and chl), followed by the remaining promotion gates.

**Next steps (in order).**
1. **1990s chl peak correction** (`build_chla_index_full.py` or a follow-on builder; new output directory, never overwrite).
   - Scale the 1993–2002 log anomalies (p90, RMA) by `k = mean(domain-median anomaly 2011, 2024) / mean(domain-median anomaly 1996, 1999)` (about (0.112 + 0.167) / (0.428 + 0.378) ≈ 0.35). This brings the 1990s peaks in line with the satellite-era peaks; apply the same k to the whole 1990s segment, both signs.
   - Alternative (lead's "vice versa"): scale the satellite era up by 1/k.
   - Recommended: scale the 1990s down. The in-water regional p90 spans stations, depths and days (spatial spread), which inflates its variance, and the satellite era has more seasons.
   - Rebuild the p90-anomaly survival table (cap 5, b 3), and plot the old and corrected multipliers for the lead.
2. **Candidate re-check** (new runner experiment, domain v0.4 to 2026, Ext×5):
   - candidates W (shares min/q50/min), C3a (Lizard q50) and C6 (μ 0.09), all with the new coral (β 1, `wc_b100`; ×0.25 below 6 DHW) and chl on (the corrected p90 table);
   - × both outbreak layers × 5 seed pairs;
   - scored with v3 (G2′ at `A½ 30`, with 10 reported); `audit_1991_candidates.py` adapted;
   - replay controls: seed-pair-1 arms against the matching runs 22 / 23 or 08 / 09 where configurations coincide.
3. **Remaining promotion gates** on the chosen set:
   - held-out years (fit 1992–2012, check 2013–2026);
   - leave-one-reef-out;
   - ±25% sensitivity (Q, μ, low-DHW scale, Ext, shares, chl b);
   - metadata and documentation;
   - update `candidates/candidates.toml`.
4. **Later:** culling back on (once the 2020s densities are plausible); the cyclone track → 4 m wave-height tool; a river-discharge predictor for pre-2003 chl, if still needed.


## 1990s chlorophyll peak correction: result (2026-10-15)

**Builder:** `build_larval_survival_index_corrected.py` writes to `sandbox/data/derived/larval_survival_index_lizard_full_1990s_corrected/` (new directory; the uncorrected tables are untouched). Figure: `figures/larval_survival_1990s_correction.png`.

- **Correction:** `k = mean(domain-median p90 anomaly 2011, 2024) / mean(1996, 1999) = 0.3459` (RMA conversion).
  - Applied to every reef's 1993–2002 p90 log anomaly, both signs. The reef long-term geometric mean is kept (no re-centring; the full-series mean anomaly shifts by about −0.01).
- **Spread of the domain-median anomaly:** 1990s SD 0.295 → 0.102. The satellite era (2003–26) is 0.066.
- **Control:** the satellite-era and missing-year multipliers reproduce the uncorrected `ap_rma` table exactly (max |Δ| 0, checked by the builder).

| Season | Before | After |
| --- | --- | --- |
| 1993 | 0.82 | 0.93 |
| 1994 | 1.31 | 1.10 |
| 1995 | 1.85 | 1.24 |
| 1996 | 3.61 | **1.56** |
| 1998 | 2.27 | 1.33 |
| 1999 | 3.11 | **1.48** |
| 2000 | 0.21 | 0.58 |
| 2001 | 0.43 | 0.75 |
| (2011, 2024, unchanged) | 1.40, 1.65 | 1.40, 1.65 |

Values are reef-median multipliers, `min(5, exp(3 · anomaly))`.


## Preregistered candidate re-check (2026-10-15, fourteenth round, step 2)

**Experiment** `recheck`: domain v0.4 (to 2026), Ext×5, Q 20, gradual Allee, `eb` seeds, q 0.015. New coral: β 1 (`wc_b100`), ×0.25 below 6 DHW. Chl on, using the 1990s-corrected p90 table (`aprmak`). Each layer (`OWEN_1991_COTS_DATASET` original / corrected) gets 16 arms:
- **W** `wc_b100_q20_eb_lin_obs015_s<k>_t60_e50_aprmak_cv2_pnd`: shares Lizard min, MacGillivray q50, North Direction min; μ 0.12;
- **C3a** `..._e50_aprmak_cv2_pndl50`: Lizard q50 share;
- **C6** `..._e50_m090_aprmak_cv2_pnd`: W with chronic μ 0.09;
- each at seed pairs k = 1–5 (`CANDIDATE_SEEDS`);
- **replay control:** `wc_b100_q20_eb_lin_obs015_s1_t60_e50_aprma_cv2_pnd` (uncorrected table) must replay the same arm of run 20261015T_22 / _23 exactly. Against W at s1, it also gives the effect of the correction.

**Scoring:**
- `rescore_detectability.jl` (score v3), and `evaluate_gate_g2prime.py` (G2′ at `A½ 30`, with 10 reported);
- audit `audit_1991_recheck.py`: per candidate × layer, the seed mean and range of v3, the manta part (v2 to 2021) and the diver error; G2′ passes out of 5; wave-3/4/5 onset and peak per reef; 2025–26 model adults against the diver/cull targets.

**Decision rules (fixed before running):**
1. **Gate.** A candidate is retained for the promotion gates if it passes G2′ at `A½ 30` in at least 4 of the 5 seed pairs on both layers.
2. **Ranking.** Retained candidates are ranked by seed-mean v3. A candidate more than 1.5 above the best seed mean is reported but not carried forward unless the lead decides otherwise.
3. **Correction effect.** If the correction worsens v3 at seed 1 by more than 0.5, it is reported to the lead before the gates. It is not reverted automatically, since it is a data correction rather than a tuning choice.
4. **Seeds.** Seed spread (max − min v3 over seeds) is reported; a spread larger than the gap between candidates means they are not distinguishable.


## Candidate re-check results (2026-10-15)

**Runs and evaluations:**
- runs 20261015T_28 (original layer) and _29 (corrected layer): 16 arms each, no failures;
- re-scores 30 / 31 (control replay delta 0.0);
- G2′ in [34](runs/20261015T_34_gate_recheck_original_lizmacg/gate_g2prime.json) / [35](runs/20261015T_35_gate_recheck_corrected_lizmacg/gate_g2prime.json), with the trough gated on Lizard and MacGillivray (lead, 2026-10-14);
- audit: [`runs/20261015T_29_.../recheck_audit.json`](runs/20261015T_29_recheck_corrected_1991_v1/recheck_audit.json) and `recheck.png` (`audit_1991_recheck.py`).

The replay arm reproduces run 22 / 23 exactly (0.0) and gives the same G2′ verdict as runs 25 / 27.

**Process error, caught and corrected.** Gates 32 / 33 were first evaluated without `--trough-reefs`, so North Direction's trough was gated as well. Every arm failed, including the replay arm, which had passed in run 25. Both gates are kept and marked SUPERSEDED. The audit outputs built from them are kept in run 29 under the prefix `SUPERSEDED_all_reef_trough_`.

**Results by candidate** (seed mean over 5 seed pairs, with the range):

| Candidate | v3, original layer | v3, corrected layer | Manta part | Diver error (orig / corr) | G2′ at `A½ 30` (orig / corr) |
| --- | --- | --- | --- | --- | --- |
| W | 12.29 (11.15–13.48) | 11.54 (10.66–12.67) | 6.83 | 3.64 / 3.14 | **5/5 / 5/5** |
| C3a (Lizard q50) | 12.15 (10.99–13.31) | 11.45 (10.38–12.49) | 6.73 | 3.61 / 3.15 | 0/5 / 0/5 |
| C6 (μ 0.09) | 13.04 (11.97–13.89) | 11.98 (10.95–13.05) | 7.23 | 3.87 / 3.16 | 2/5 / 2/5 |

**Findings.**
1. **Only W passes the preregistered gate (≥ 4/5 seeds on both layers); W is retained.**
   - **C3a:** with the new coral, Lizard's 1990s wave peaks late (1999, detected 0.5–0.65 COTS/tow) and stays high through the 2000s. The trough is 0.11–0.19 against an upper bound of 0.054. In 3 of 5 seeds the 2010s wave does not separate.
   - **C6:** Lizard's 2010s wave is faint (detected 0.04–0.06), and in 3 seeds a peak is missed.
   - C3a's v3 is marginally lower than W's (−0.1), which is well inside the seed spread (about 2.0).
2. **Seed spread (about 2.0 in v3) exceeds the gaps between candidates.** v3 alone cannot separate them; the qualitative gate does.
3. **The 1990s chl correction is neutral to slightly favourable:** W at seed 1, v3 −0.01 (original) and −0.05 (corrected); manta part +0.03; diver error −0.03 / −0.06.
4. **The 2020s by layer** (W, seed mean; model ÷ diver/cull target):

   | Layer | Reef | 2025 | 2026 |
   | --- | --- | --- | --- |
   | Original | Lizard | 0.35 | 0.97 |
   | Original | North Direction | 0.49 | 1.13 |
   | Original | Eyrie | 0.16 | 0.08 |
   | Corrected | Lizard | 1.8 | 7.9 |
   | Corrected | North Direction | 1.3 | 7.2 |
   | Corrected | Eyrie | 1.3 | 0.8 |

   - **Original layer:** flat at about 3/ha through 2020–26 at Lizard (targets 6–9). There is no build-up at Eyrie.
   - **Corrected layer:** tracks the 2020–25 build-up at Lizard, North Direction and Eyrie, then overshoots in 2026: 40–55/ha at Lizard, MacGillivray and North Direction, against targets of 5–7. The observed 2025 → 2026 densities fall at Lizard (9.1 → 5.0) and North Direction (10.6 → 7.2).
   - Culling is off in these runs. Direct removals on Lizard's main reef were about 1 COTS/ha/yr in 2024–26, too small to explain the gap alone. Several nearby reefs were culled at 10–70 COTS/ha/yr (e.g. 15-028, 15-084, 15-071c, 14-163), which would reduce larval supply reaching the calibration reefs.


## Lead's decisions (2026-10-15, fifteenth round)

1. **Test observed culling on W before the promotion gates** (W × 5 seeds × both layers, cull on vs off).
2. **Promotion-gate definitions for screened candidates agreed:**
   - **held-out years:** redo the screening choices on the 1992–2012 fit only, check W is still chosen, and report 2013–26;
   - **leave one reef out:** the same, dropping one calibration reef at a time;
   - **±25% sensitivity:** on W, 5 seeds.
3. **New test: an elevated 1991 starting density.** The lead notes that the model starts well below its own inter-wave level ("it's not like COTS came from nowhere"). Start nearer trough-like densities and check whether the 1990s peaks come forward.
4. Work up to the re-check was committed first (`9f6a252e`).

**Prior evidence.** The initial-state re-baseline (2026-10-13, old coral, β 1.5, before shares, new coral and chl) found the start barely shaped first-wave timing or height, and a large North Direction seed (7.6/ha) grazed its coral early and delayed its wave. The model has changed substantially since, so the test is repeated on W.


## Preregistered starting density × observed culling (2026-10-15, fifteenth round)

**Hypotheses.**
- **Start.** In 1989–91 the calibration reefs held an inter-outbreak remnant population that manta tows under-detect at low density. Reading the 1989–91 tows with the same low-density detection curve used by the G2′ observation model raises the start to trough-like densities (about 2–6 adults/ha, around the 3/ha Allee scale). Local production can then begin before external supply arrives, so the 1990s onset and peak come earlier.
- **Culling.** The observed cull removals (2002–2026; 21,000–23,000 COTS domain-wide in each of 2025 and 2026) reduce adults where they are applied and larval output from heavily culled source reefs. This should lower the corrected layer's 2026 overshoot.

**Change.** No model change. New runner experiment `init_cull` on W (new coral, β 1, ×0.25 below 6 DHW; 1990s-corrected chl `aprmak`; Ext×5; Q 20; q 0.015; v0.4 to 2026):
- **start inversion** of the `eb` 1989–91 surface:
  - `lin`: `A = C / q` (current);
  - `h10`, `h30`: `A = (C/q + sqrt((C/q)² + 4·A½·C/q)) / 2` with `A½` 10 or 30;
  - juvenile and subadult preloads keep the 6:3:1 ratio;
- **culling:** off, or on with the observed removals at σ 210 m (`removals_sigma210.csv`; pre-spawning culls before the COTS step, post-spawning culls after it);
- × 5 seed pairs × both outbreak layers, giving 30 arms per layer;
- arm names `wc_b100_q20_eb_<lin|h10|h30>_obs015_s<k>_t60_e50[_cull210]_aprmak_cv2_pnd`.

**Approximate reef-mean starting adults/ha** (from the 1991 seeds; the runner reports the exact values):

| Reef | `lin` | `h10` | `h30` |
| --- | --- | --- | --- |
| Lizard | 0.47 | about 2.3 | about 4 |
| MacGillivray | 0.35 | about 2 | about 3.5 |
| North Direction | 0.82 | about 3 | about 5 |

**Controls.** The `lin`, no-cull arms carry the exact names of the run 28 / 29 W arms and must replay all five seeds exactly (trajectories, flows, mechanisms).

**Expected.**
- `h10` / `h30`: adults ≥ 3/ha 1–2 years earlier; first-wave peak year earlier by 0–1 years, most at North Direction (observed 1995) and MacGillivray (1996); larger 1990s peaks.
  - Risk: early coral depletion at North Direction (as on 2026-10-13), giving a smaller or later first wave.
  - Little effect after 2000.
- Culling: no effect before 2002; small effects 2012–16; lower 2025–26 adults on the corrected layer; lower diver error if the 2026 overshoot falls.

**Read-out** (audit `audit_1991_init_cull.py`):
- seed mean and range of v3, the manta part and the diver error;
- G2′ passes at `A½ 30` (10 reported), with the trough gated on Lizard and MacGillivray;
- per reef: start adults, first-wave onset (≥ 3/ha), peak year and height, and the mean absolute first-peak-year error against the observed raw peaks (Lizard 1998, MacGillivray 1996, North Direction 1995);
- 2025–26 model ÷ target;
- removals requested vs applied (unmet fraction) for the calibration reefs and the domain.

**Decision rules (fixed before running):**
1. **Start.** An elevated start is recommended to the lead only if all of these hold:
   - it passes G2′ at `A½ 30` in ≥ 4/5 seeds on both layers;
   - its seed-mean v3 is no more than 0.5 worse than `lin`;
   - it reduces the mean absolute first-peak-year error.
2. **Culling.** Report its effect on the 2020s diver error and the 2026 overshoot. As in run 20261014T_42, an unmet fraction above 25% on a calibration reef limits interpretation, and the result goes to the lead.
3. **Interaction.** Only reported. The start acts in the 1990s and culling from 2002, so little interaction is expected.


## Starting density × observed culling results (2026-10-15)

**Runs and evaluations:**
- runs 20261015T_36 (original layer) and _37 (corrected layer): 30 arms each, no failures;
- re-scores 38 / 40 (control replay delta 0.0);
- G2′ in 39 / 41, trough gated on Lizard and MacGillivray;
- audit: [`runs/20261015T_37_.../init_cull_audit.json`](runs/20261015T_37_init_cull_corrected_1991_v1/init_cull_audit.json) and `init_cull.png` (`audit_1991_init_cull.py`).

All ten `lin`, no-cull arms replay runs 28 / 29 exactly (0.0).

**Start** (no cull; seed means; both layers are identical before 2018):

| Start | Lizard / MacGillivray / North Direction start (adults/ha) | First-peak year (observed 1998 / 1996 / 1995) | First-peak height (adults/ha) | Mean first-peak-year error | v3 (orig / corr) | G2′ `A½ 30` |
| --- | --- | --- | --- | --- | --- | --- |
| `lin` | 0.47 / 0.35 / 0.82 | 1997.0 / 1997.6 / 1997.2 | 30 / 46 / 71 | 1.60 y | 12.29 / 11.54 | 5/5, 5/5 |
| `h10` | 2.41 / 2.05 / 3.30 | 1996.8 / 1997.2 / 1996.8 | 25 / 39 / 66 | 1.40 y | 12.11 / 11.37 | 5/5, 5/5 |
| `h30` | 3.98 / 3.42 / 5.39 | 1996.0 / 1996.6 / 1995.8 | 19 / 33 / 46 | 1.13 y | 13.14 / 12.40 | 4/5, 4/5 |

- **`h10` meets all three preregistered rules:** G2′ 5/5 on both layers; v3 −0.17; peak-year error −0.2 y.
  - Paired by seed: v3 is lower in 4 of 5 seeds, and the peak-year error is lower or equal in every seed.
  - Onset (≥ 3/ha) comes a year earlier at North Direction (1993) and MacGillivray.
- **`h30` fails the v3 rule (+0.86)** although its timing is best.
  - In seed 2, North Direction loses its first wave entirely (its large start grazes the coral early, as on 2026-10-13), and the v3 for Lizard and North Direction doubles. G2′ fails in that seed.
  - It also moves Lizard's peak earlier than observed (1996 vs 1998).
- **Trade-off:** a higher start brings first waves earlier but makes them smaller. Effects vanish after about 2000 (wave-4 and wave-5 peaks within a few percent).

**Culling** (σ 210 m, against no cull):

| | Original layer | Corrected layer |
| --- | --- | --- |
| v3 change (seed mean, `lin`) | +0.76 | +0.32 |
| Diver-error change | +0.43 | +0.14 |
| Lizard 2026 model ÷ target | 0.97 → 0.82 | 7.9 → 7.6 |
| North Direction 2026 | 1.13 → 0.95 | 7.2 → 6.9 |
| Unmet removals 2022–26: Lizard / Eyrie / domain | 34% / 82% / 73% | 6% / 57% / 51% |
| Unmet removals 2012–16: Lizard / North Direction / domain | 3% / 21% / 57% | 3% / 21% / 57% |

- **Culling worsens the fit on both layers and barely touches the 2026 overshoot.**
- The preregistered limit is triggered: unmet removals exceed 25% at Eyrie (both layers) and at Lizard (original layer). On the original layer the model holds too few COTS at the cull sites for the removals to be taken.

**Source of the 2026 overshoot (corrected layer).** It is external supply, not local build-up. At Lizard, external settlement (W, seed 1) runs 14–20 recruits/ha per year in 2024–26 on the original layer, but 60, 133 and 148 on the corrected layer. North Direction shows the same: 86, 163, 179 against 22–30.
- These are the years of the heaviest culling nearby, so the cull-informed correction raises outbreak probability most there, and Ext×5 then multiplies it.
- Recruits settled in 2024–25 mature into the 2026 adults.

**Status.** Nothing is promoted. The decisions go to the lead:
1. adopt the `h10` start into W;
2. keep culling off in calibration;
3. how to treat the corrected layer's 2024–26 external supply.


## Lead's decisions (2026-10-15, sixteenth round) and the external-supply clarification

**Decisions:**
1. **Adopt the `h10` start into W.** The 1989–91 `eb` surface is inverted with the low-density detection curve at `A½ = 10`.
2. **Culling stays off in calibration for now,** but it is the most important management lever, so the model must be made to work with it. Open item: the model must hold the COTS that were actually culled; unmet removals were 57–82% at Eyrie in 2022–26.
3. **Corrected layer:** carry both layers as a bracket for now, and test option (b), where the ×5 boost is not applied to the part of the supply added by the diver/cull correction.

**Clarification: where external larvae come from** (diagnostic `analyse_external_sources.py`; outputs in the git-ignored `sandbox/data/derived/external_supply_sources/`).
- The recomputed boundary matches the dataset files ADRIA reads (max relative difference 8e-16).

**What the code does:**
- **External supply comes only from reefs that are not among the 113 modelled domain reefs.** For each domain reef `i` and year `y`, the builder sums over every non-domain GBR reef `j` with RME area and survival data: `B[i, y] = Σ_j C[j→i] · A_j / A_i · s_j · P_j(y)`.
  - `C`: Owen connectivity for the hydrodynamic year drawn that year (seeded);
  - `A`: habitable area;
  - `s`: RME larval survival;
  - `P`: predicted probability that manta CPUE > 0.22.
- ADRIA multiplies `B` by the outbreak production constant × Ext (×5) and adds the result as settling recruits.
- **The predictions are never used for the 113 domain reefs.** Their COTS are simulated. Internal supply (between domain reefs) comes from simulated densities through the same connectivity matrices.
- **The corrected layer** adds IDW-spread residuals (within 50 km) from SALAD and from GBR-wide cull dives to every GBR reef's predicted probability. It changes supply only through the non-domain reefs.

**"External" does not mean far away.**
- The 113 modelled reefs are a subset of the reefs in the area. 24 unmodelled reefs lie inside the footprint of the 113 (their convex hull), including **Lizard Island Reef (North West), 14-116a, 2.4 km from Lizard**, and Nymph Reef (23 km). These reefs are boundary sources driven by the predictions, not by simulated dynamics.
- The modelled reefs reach 118 km from Lizard.
- Of Lizard's outbreak-free connection weight, 9% comes from within 50 km, 56% from within 100 km and 92% from within 200 km. The strongest connections run from reefs to the north-west (14-0xx, 50–100 km) and to the south (15-0xx–16-0xx, 100–200 km).
- Unmodelled reefs inside the footprint supply 0–30% of a calibration reef's external larvae (e.g. Lizard: 25% in 1995, 8% in 2005, 20% in 2012, 18–20% in 2025). The rest comes from beyond the footprint, almost all within 200 km. Sources beyond 500 km are negligible.

**External supply dominates settlement at Ext×5.** Share of settling recruits from external supply (W, seed mean; the remainder is internal, own-reef retention or background):

| Reef | 1992–2000 | 2001–11 | 2012–16 | 2017–21 | 2022–26, original / corrected |
| --- | --- | --- | --- | --- | --- |
| Lizard | 73% | 71% | 87% | 76–79% | 63% / 74% |
| MacGillivray | 67% | 60% | 84% | 62–67% | 55% / 70% |
| North Direction | 72% | 70% | 88% | 80% | 64% / 74% |
| Eyrie | 85% | 84% | 92% | 87–91% | 76% / 86% |

At Ext×1 the local share was much higher (Lizard 2019–26: 78% local against 33% at ×5; run 20261015T_15).

**What drives the 2025–26 overshoot on the corrected layer.** The top Lizard sources in 2025 are Endeavour (15-089, 123 km), U/N 15-077 (100 km), Batt (16-029, 196 km), Pickersgill (15-093, 134 km) and Evening (15-095, 138 km). The correction raises their predicted outbreak probabilities from 0.01–0.04 to 0.41–0.98 (heavily culled reefs between Cooktown and Port Douglas). The top 5 sources then carry 57% of Lizard's external supply, against 29% on the original layer. MacGillivray, North Direction and Eyrie share the same southern sources.

**Implication for option (b).** Because `B` is linear in `P`, (b) can be built without any model change, by a new derived prediction layer `P* = P_orig + (P_corr − P_orig) / 5`. Under Ext×5 this gives `5·P_orig + (P_corr − P_orig)`: the manta-based predictions keep the ×5 boost, and the diver/cull correction enters at face value. `P*` stays within [0, 1].


## Preregistered face-value correction, option (b) (2026-10-15, sixteenth round)

**Question (lead):** how large is the double count when the ×5 external boost multiplies the diver/cull-corrected part of the outbreak layer?

**Input** (no model change):
- `build_outbreak_layer_face_value.py` builds `P* = P_orig + (P_corr − P_orig) / 5` (`sandbox/data/derived/outbreak_layer_face_value_x5_2018_2026/`).
- `build_lizard_cots_connectivity_owen.py` builds v0.4 dataset `owen_global_2018_2023_facevalue` from it, with the same policies as the other two layers.
- **Checks, all passed:**
  - the external boundary equals `B_orig + (B_corr − B_orig)/5` for every reef, year and hydrodynamic year (max relative difference 2e-14);
  - every other forcing file (internal connectivity, survival, site map, external supply coefficients) is byte-identical to the original dataset;
  - years before 2018 are unchanged.

**Arms:** W with the adopted `h10` start, no cull, × 5 seed pairs on the face-value layer. This uses runner experiment `init_cull` with `OWEN_1991_ARM_FILTER='eb_h10_obs015_s\d_t60_e50_aprmak_cv2_pnd'`. The comparators are the same arms on the original and corrected layers (runs 36 / 37).

**Expected:**
- 2024–26 external settlement at Lizard falls from about 60 / 133 / 148 (corrected) toward 20–40 recruits/ha;
- the 2026 overshoot falls from about 7.8× to roughly 1.5–3× target;
- the 2020–25 build-up is weaker than on the corrected layer but stronger than on the original;
- the diver error lies between the two layers;
- nothing changes before 2018.

**Read-out:** the seed mean and range of v3, the manta part and the diver error; G2′ (trough on Lizard and MacGillivray); 2022–26 model ÷ target per calibration reef; external settlement 2022–26. Diagnostic only; nothing is adopted without the lead.


## Face-value correction, option (b): results (2026-10-15)

**Run and evaluations:**
- run 20261015T_42 (face-value layer, W + `h10`, no cull, 5 seeds, no failures);
- re-score 43, G2′ 44 (trough on Lizard and MacGillivray);
- audit [`facevalue_audit.json`](runs/20261015T_42_facevalue_h10_1991_v1/facevalue_audit.json) and `facevalue.png` (`audit_1991_facevalue.py`);
- comparators: the same arms in runs 36 / 37.

| Layer | v3 (seed mean, range) | Manta part | Diver error | G2′ `A½ 30` |
| --- | --- | --- | --- | --- |
| Original | 12.11 (11.08–13.12) | 6.73 | 3.59 | 5/5 |
| **Face value (b)** | **11.64 (10.71–12.71)** | 6.73 | 3.28 | 5/5 |
| Corrected | 11.37 (10.16–12.27) | 6.73 | 3.09 | 5/5 |

Paired by seed, face value lies between the two layers in 4 of 5 seeds (best of all in seed 1).

**Adults/ha, seed mean (target in brackets):**

| Reef | 2022 | 2023 | 2024 | 2025 | 2026 |
| --- | --- | --- | --- | --- | --- |
| Lizard, original / **face value** / corrected | 3.3 / **3.3** / 3.9 (8.3) | 3.2 / **3.6** / 6.3 (6.0) | 3.0 / **3.9** / 9.7 (8.6) | 2.9 / **4.9** / 17.4 (9.1) | 4.9 / **10.7** / 44.0 (5.0) |
| North Direction | 3.2 / **3.3** / 3.4 (3.7) | 4.3 / **4.5** / 6.0 (2.5) | 4.8 / **5.5** / 9.2 (12.8) | 5.0 / **6.5** / 14.4 (10.6) | 8.3 / **15.5** / 57.4 (7.2) |
| Eyrie | 1.5 / **1.8** / 3.7 (16.9) | 1.6 / **2.3** / 7.0 (45.2) | 1.7 / **2.9** / 11.1 (16.3) | 1.6 / **3.6** / 17.1 (12.2) | 2.2 / **5.6** / 28.9 (34.0) |

External settlement at Lizard (recruits/ha) in 2024 / 2025 / 2026: original 11 / 15 / 14; face value 21 / 33 / 33; corrected 62 / 104 / 112.

**Findings.**
1. **The double count is large.** Most of the corrected layer's 2020s build-up and of its 2026 overshoot comes from applying ×5 to the diver/cull correction.
   - At face value, Lizard 2026 falls from 44 to 10.7 adults/ha, about 75% of the corrected layer's excess over the original.
   - External settlement in 2024–26 falls to about a third.
2. **At face value the 2020s build-up is mostly missing.** Lizard, North Direction and Eyrie stay well below the 2022–25 targets, then rise late and overshoot in 2026 (about 2×).
   - Eyrie reaches 2–6/ha against 12–45.
   - The shape is consistent with supply arriving too late: settlers need about two years to mature, so the observed 2022–24 build-up needs recruitment in about 2019–22, where the correction is small.
3. Face value is intermediate on v3 (11.64), best on the Lizard diver error (0.53 vs 0.58 / 0.67), and passes G2′ in 5/5 seeds.


## Self-recruitment: why it is low (2026-10-15)

**Lead's question:** is low self-recruitment coming from the connectivity layer, or from a parameter we set? Diagnostic `analyse_self_recruitment.py`; outputs in the git-ignored `sandbox/data/derived/self_recruitment_diagnostic/`.

**Two different things:**

**(1) Own-reef self-recruitment is set by the hydrodynamic layer alone.** In COTSMod (separated recruitment), own-reef settlement = own larval output × the Owen diagonal `C[i, i]` × settlement. No parameter scales it: the immigration scalar and settlement probability are 1, and the settlement gate is the food/condition gate. Values (mean of 2018–22 hydrodynamic years):

| Reef | Own larvae returning to the reef | Settling on any GBR reef | Settling on the 113 domain reefs |
| --- | --- | --- | --- |
| Lizard (3 units combined; single units 0.2%) | 0.33% | 16% | 7.1% |
| MacGillivray | 0.07% | 15% | 7.0% |
| North Direction | 0.06% | 15% | 7.9% |
| Eyrie | 0.33% | 17% | 5.0% |

- The earlier retention screen (2026-10-13) raised Lizard's retention up to 90% of its surviving larvae (source-conserving) and changed Lizard's first wave by only +3%.

**(2) The low local share of settlers (external 60–90%) comes mainly from parameters, not connectivity.** If every reef produced the same larvae per hectare, the connectivity alone gives each calibration reef's arrivals as:

| Reef | Self | Other domain reefs | External reefs |
| --- | --- | --- | --- |
| Lizard | 2.1% | 74% | 24% |
| MacGillivray | 0.3% | 83% | 17% |
| North Direction | 0.2% | 79% | 21% |
| Eyrie | 3.6% | 59% | 37% |

- So the hydrodynamics say **the domain itself should supply about two-thirds to four-fifths of the larvae.** External supply dominates in the model because of the production asymmetry:
  - **Local:** 37.9 pre-pelagic recruits per adult per year, further reduced by body condition², the gradual Allee term (×0.37 at 1/ha, ×0.64 at 15/ha) and the Ricker term.
  - **External:** an outbreak reef at P = 1 produces 2,798 recruits/ha/yr, the output of about 74 adults/ha at full condition and no Allee reduction. Ext×5 makes it about 369 adults/ha. At P = 0.5 that is still about 185 adults/ha.
  - Domain reefs between waves hold 2–8 adults/ha.
- The external production constant was matched in an earlier screen (`matched_boundary_pre_pelagic_ha_year`) and migrated with the density unit (k ≈ 19.6); Ext×5 is the lead's working assumption (2026-10-12).

**Implication.** The model's local dynamics are under-weighted by the size of the external boundary term, not by the hydrodynamic connectivity. Option (b) shows that removing the double count leaves the 2020s build-up short. Together, these point to local (domain) production rather than more external supply as the missing source. This is consistent with the lead's initiation-region hypothesis. A labelled diagnostic could test it: lower Ext with higher local fecundity, or the external production constant re-derived as an outbreak density (e.g. 15–50 adults/ha) instead of the matched value.


## Lead's decisions (2026-10-15, seventeenth round)

1. **The face-value layer (`owen_global_2018_2023_facevalue`) is the main outbreak layer.** The original layer is kept as the low case, and the corrected layer as a labelled high case only.
2. **Run the local-vs-external production diagnostic** (below).


## Preregistered local vs external production diagnostic (2026-10-15, seventeenth round)

**Hypothesis (lead's initiation-region view).** The hydrodynamics give the domain 59–83% of each calibration reef's arrivals under equal production. The model nonetheless relies on external supply because an outbreaking non-domain reef counts as about 74 adults/ha (×5 ≈ 369/ha), while domain reefs produce 37.9 recruits per adult. Shifting weight from external to local production should let local build-up, recruitment years ahead of the peak, start the 2020s outbreak earlier with less external supply.

**Change.** No model change. New runner experiment `local_ext` on W + `h10` start, no cull, face-value layer, v0.4 to 2026:
- external boost Ext {×5, ×2.5};
- × local larval output multiplier on `COTS_LARVAL_FECUNDITY` {×1, ×2, ×3};
- × 5 seed pairs, giving 30 arms;
- arm names `wc_b100_q20_eb_h10_obs015_s<k>_t60_e<50|25>[_f<200|300>]_aprmak_cv2_pnd` (no `_f` token at ×1).

The multiplier scales potential larvae only; the Allee, Ricker, condition and settlement terms are unchanged. **Labelled diagnostic.** Local output per adult is a calibrated, not measured, value, so nothing is adopted without evidence and the lead's decision.

**Control.** The Ext×5, ×1 arms carry the run 42 names and must replay them exactly.

**Expected:**
- **Local output ×2–3 at Ext×5:** larger and earlier local build-up in all waves; overshoot risk in waves 3–4 (manta part worsens); the 2022–25 build-up arrives earlier and the 2026 overshoot grows.
- **Ext×2.5 at ×1:** as in run 20261015T_22, a later 1990s wave and G2′ failures.
- **Ext×2.5 with ×2–3:** local output partly replaces external supply. The 1990s onset may partly recover; the local share of settlers rises.
- **The key comparison:** whether any Ext×2.5 arm matches or beats Ext×5 ×1 on v3 and passes G2′.

**Read-out** (audit `audit_1991_local_ext.py`):
- seed mean and range of v3, the manta part and the diver error;
- G2′ at `A½ 30` (trough on Lizard and MacGillivray);
- per reef: wave-3/4/5 onset (≥ 3/ha), peak year and height;
- 2022–26 adults against targets;
- local (domain-produced) share of settlers by period.


## Local vs external production diagnostic: results (2026-10-15)

**Run and evaluations:**
- run 20261015T_45 (30 arms, face-value layer, no failures);
- re-score 46, G2′ 47 (trough on Lizard and MacGillivray);
- audit [`local_ext_audit.json`](runs/20261015T_45_local_ext_facevalue_1991_v1/local_ext_audit.json) and `local_ext.png` (`audit_1991_local_ext.py`).

The Ext×5 ×1 arms replay run 42 exactly (0.0).

| Ext | Local output | v3 (seed mean, range) | Manta part | Diver error | G2′ `A½ 30` | First-peak-year error | Lizard domain share of settlers 2019–26 | Lizard 2024 / 2025 / 2026 ÷ target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ×5 | ×1 | **11.64** (10.71–12.71) | 6.73 | 3.28 | **5/5** | 1.40 y | 25% | 0.47 / 0.55 / 2.0 |
| ×5 | ×2 | 13.13 | 7.65 | 3.65 | 4/5 | 1.27 y | 42% | 0.62 / 0.79 / 3.2 |
| ×5 | ×3 | 14.71 | 8.80 | 3.93 | 1/5 | 1.20 y | 52% | 0.71 / 0.99 / 4.5 |
| ×2.5 | ×1 | 15.26 | 8.88 | 4.26 | 0/5 | 2.93 y | 40% | 0.42 / 0.37 / 0.9 |
| ×2.5 | ×2 | 15.72 | 9.47 | 4.16 | 1/5 | 2.40 y | 55% | 0.58 / 0.55 / 1.7 |
| ×2.5 | ×3 | 16.00 | 9.95 | 4.03 | 0/5 | 1.67 y | 64% | 0.68 / 0.72 / 2.7 |

v3 is worse than Ext×5 ×1 in every seed for every other arm.

**Findings.**
1. **More local production does not bring the 2020s build-up forward; it amplifies it after it starts.**
   - The 2022–25 deficit narrows only slightly (Lizard 2024: 0.47 → 0.71 of target).
   - The 2026 overshoot grows (2.0 → 4.5×).
   - Wave-4 peaks fall and wave-3 peaks rise; the manta part worsens.
   - Local output does not substitute for external supply in the 1990s either: at Ext×2.5 the first peaks stay late (Lizard about 1998–99, North Direction 1998).
   - Eyrie stays at 0.05–0.4 of target in every arm.
2. **The binding constraint in 2019–23 is the juvenile predator pit, not larval supply.**
   - At Lizard (W, seed 1, face value), 10–40 recruits/ha settle each year in 2020–23. About 59% survive their first year into subadults.
   - Each year the pit (`Q` = 20 COTS/ha/yr, removal `Q·J/(J + Q)`) eats 54–57% of subadults, background mortality takes about 31%, and only 12–15% mature: about 1 adult/ha/yr.
   - Adults hold at about 4/ha until juveniles swamp the pit (2026: pit share 36%).
   - In 2015–18 the limit is starvation instead (coral 5–8%, food survival 0.06–0.09).
   - This is why the double-counted corrected layer appeared to work: its supply swamped the pit.
3. `Q` = 20 was chosen in the predation-capacity grid (2026-10-12) and re-checked in the timing experiment (`Q` {10, 15, 20}, 2026-10-14) under older configurations (old coral, no chl, v0.3). It keeps the inter-wave troughs low (G2′ component 4), so lowering it trades trough fit for the 2020s build-up.

**Status.** Labelled diagnostic, nothing adopted. Local output stays ×1 and Ext stays ×5. Next for the lead: re-check the predator-pit capacity on the current configuration (W + `h10`, face-value layer).


## Lead's decisions (2026-10-16, eighteenth round)

1. **Re-check the predator-pit capacity `Q`** on the current configuration.
2. **Consider a weaker corrected layer.** The CPUE-to-density conversion behind the correction may be too strong: cull CPUE reflects densities before culling, not what was left to reproduce. A partly reduced correction might let COTS build up without overshooting, and culling, once implemented, should mediate the rest.


## Preregistered predator-pit capacity × correction weight (2026-10-16, eighteenth round)

**Hypotheses:**
- The juvenile predator pit (`Q` 20) is the binding constraint on the 2019–23 build-up (it eats 54–57% of Lizard subadults a year). A weaker pit lets moderate supply build up sooner.
- A correction weight between face value and fully boosted gives enough supply without the overshoot.
- Both act on the same build-up, so they are crossed.

**Correction weight `w`.**
- `build_outbreak_layer_face_value.py` gains `--weight w`: `P_w = P_orig + w·(P_corr − P_orig)/5`, so supply under Ext×5 = `5·P_orig + w·(P_corr − P_orig)`.
  - The default (`w` = 1) reproduces the face-value layer exactly (checked).
  - `w` = 5 would reproduce the fully boosted corrected layer.
- New v0.4 datasets `owen_global_2018_2023_corr_w2` and `_corr_w3`. The boundary identity holds (max relative difference 2e-14), and all other forcing files are byte-identical to the original dataset.

**Arms.** Runner experiment `pit_q`:
- W + `h10` start, no cull, local output ×1, Ext×5;
- × `Q` {10, 15, 20, 25} COTS/ha/yr × 5 seed pairs;
- arms `wc_b100_q<Q>_eb_h10_obs015_s<k>_t60_e50_aprmak_cv2_pnd`;
- run once per dataset: face value (`w` 1), `corr_w2`, `corr_w3`. That is 60 arms.

**Control.** The `Q` 20 arms on the face-value dataset must replay run 42 exactly.

**Expected:**
- **Lower `Q`:**
  - more subadults mature, so adults rise earlier in 2020–23 and the 2022–25 deficit shrinks;
  - higher inter-wave densities (G2′ trough risk at Lizard and MacGillivray);
  - larger wave peaks (manta part may worsen);
  - a larger 2026 overshoot.
- **Higher `w`:** more 2022–26 supply, an earlier build-up and a larger 2026 overshoot; no change before 2018.
- A combination such as `Q` 15 with `w` 2 may reach the 2022–25 targets with a modest 2026 overshoot.

**Read-out** (audit `audit_1991_pit_q.py`):
- per (`w`, `Q`): seed mean and range of v3, the manta part and the diver error;
- G2′ passes at `A½ 30` (trough on Lizard and MacGillivray), with the model inter-wave CPUE at Lizard and MacGillivray against the bound;
- 2022–26 adults against targets at Lizard, North Direction, MacGillivray and Eyrie;
- the share of subadults eaten by the pit at Lizard, 2020–23.

**Selection rule (fixed before running).** A combination is recommended only if all of these hold:
- it passes G2′ at `A½ 30` in ≥ 4/5 seeds;
- its seed-mean v3 is lower than the current configuration (`w` 1, `Q` 20: 11.64);
- its diver error is lower.

Among the qualifiers, prefer the smallest change from `Q` 20 and `w` 1. Culling remains off; the 2026 overshoot is read in that light.


## Predator-pit capacity × correction weight: results (2026-10-16)

**Runs and evaluations:**
- runs 20261016T_01 / _02 / _03 (`w` 1 / 2 / 3; 20 arms each, no failures);
- re-scores 04 / 06 / 08, G2′ 05 / 07 / 09 (trough on Lizard and MacGillivray);
- audit [`pit_q_audit.json`](runs/20261016T_01_pit_q_w1_1991_v1/pit_q_audit.json) and `pit_q.png` (`audit_1991_pit_q.py`).

The `Q` 20 arms on the face-value dataset replay run 42 exactly (0.0).

| `w` | `Q` | v3 (seed mean) | Manta part | Diver error | G2′ `A½ 30` | Lizard 2023 / 2024 / 2025 / 2026 ÷ target | Eyrie 2025 ÷ target | Qualifies |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 10 | 13.24 | 8.13 | 3.41 | 4/5 | 0.80 / 0.67 / 0.86 / 3.2 | 0.50 | no |
| 1 | 15 | 11.70 | 6.75 | 3.30 | 5/5 | 0.74 / 0.56 / 0.68 / 2.5 | 0.39 | no |
| **1** | **20** | **11.64** (current) | 6.73 | 3.28 | 5/5 | 0.63 / 0.47 / 0.55 / 2.0 | 0.32 | — |
| 1 | 25 | 11.53 | 6.48 | 3.37 | 5/5 | 0.59 / 0.41 / 0.47 / 1.6 | 0.27 | no |
| 2 | 15 | 11.48 | 6.75 | 3.16 | 5/5 | 0.82 / 0.69 / 0.98 / 4.0 | 0.63 | yes |
| **2** | **20** | **11.36** | 6.73 | 3.09 | 5/5 | 0.70 / 0.59 / 0.82 / 3.4 | 0.53 | **yes** |
| 2 | 25 | 11.19 | 6.48 | 3.14 | 5/5 | 0.65 / 0.51 / 0.69 / 2.8 | 0.45 | yes |
| 3 | 15 | 11.42 | 6.75 | 3.12 | 5/5 | 0.89 / 0.83 / 1.31 / 5.7 | 0.91 | yes |
| 3 | 20 | 11.26 | 6.73 | 3.02 | 5/5 | 0.76 / 0.72 / 1.12 / 4.9 | 0.77 | yes |
| 3 | 25 | **11.01** | 6.48 | 3.02 | 5/5 | 0.70 / 0.63 / 0.96 / 4.2 | 0.67 | yes |

The `w` 2 and 3 rows at `Q` 10 match the `w` 1 pattern (v3 about 13.1–13.2, G2′ 4/5).

**Findings.**
1. **Lowering the pit capacity does not fix the 2020s.** `Q` 10 halves the capacity but only lowers the share of Lizard subadults eaten in 2020–23 from 56% to 50%. It lifts Lizard 2022–25 by 15–40%, but:
   - the manta part worsens (6.73 → 8.13): the 2010s wave ends early and MacGillivray sits above its 2020–25 targets;
   - seed 2 collapses (v3 17.5) and fails G2′.
   - `Q` 25 is slightly better than 20 on v3 (−0.1 to −0.25, through the manta part) but slows the 2020s build-up.
   - **Correction to the previous section:** the pit is a large sink, but not the single binding constraint. The 2022 deficit at Lizard (model about 3.5–4 against 8.3) is unchanged in every arm; it reflects the few recruits of 2019–20, after starvation in 2015–18.
2. **A stronger correction weight improves the fit consistently.** Paired by seed at `Q` 20, `w` 3 beats `w` 1 in all 5 seeds, and `w` 2 in 4.
   - `w` 3 brings Lizard 2023–25 close to target (0.76 / 0.72 / 1.12) and Eyrie from 0.32 to 0.77 of target in 2025.
   - The 2026 overshoot grows (Lizard 2.0 → 3.4 → 4.9×), with culling still off.
3. **Under the preregistered rule** (≥ 4/5 G2′, and v3 and diver error both below `w` 1 `Q` 20), the qualifiers are `w` {2, 3} × `Q` {15, 20, 25}. The rule's smallest change from the current configuration is **`w` 2, `Q` 20** (v3 11.36, diver error 3.09). The lowest v3 is `w` 3, `Q` 25 (11.01).

**Status.** Nothing is adopted. Correction weight `w` 2 (or 3) with `Q` 20 goes to the lead.


## Lead's decisions (2026-10-16, nineteenth round)

1. **Adopt correction weight `w` 2** (`owen_global_2018_2023_corr_w2`) as the main outbreak layer, with `w` 3 as a sensitivity. Predator-pit `Q` stays 20.
2. **Re-check observed culling first**, then the promotion gates.
3. **Refresh the LTMP manta tables.** There are new manta-tow data for the Lizard sector; check `COTS-Coral-OutcomesV2/R/aims_api.R` first, then refresh.
4. **Extend the calibration set:**
   - **Martin and Linnet reefs:** the 2020s increase shows up in manta tows there.
   - **Yonge, Carter and No Name reefs:** these almost never have COTS; check the model can reproduce that absence.


## Preregistered cull re-check on the `w` 2 layer (2026-10-16, nineteenth round)

**Question.** With more COTS in the 2020s (`w` 2), can the model now take the observed cull removals, and do they bring the 2026 overshoot down? In run 20261015T_36 / _37, 57–82% of Eyrie's 2022–26 removals were unmet.

**Arms.** Runner experiment `init_cull` with `OWEN_1991_ARM_FILTER='eb_h10_obs015_s\d_t60_e50_(cull210_)?aprmak_cv2_pnd'` on dataset `owen_global_2018_2023_corr_w2`:
- W + `h10`, `Q` 20, Ext×5;
- culling {off, observed σ 210 m} × 5 seed pairs = 10 arms.

**Control.** The no-cull arms must replay the `Q` 20 arms of run 20261016T_02 exactly.

**Expected:**
- unmet removals on the calibration reefs fall below the 2026-10-15 levels;
- 2025–26 adults fall where culls are heavy (Eyrie, and the reefs near Lizard);
- the 2026 overshoot shrinks but does not vanish (direct removals on Lizard's main reef are about 1 COTS/ha/yr).

**Read-out:** as `audit_1991_init_cull.py`, plus unmet removals by reef. Scoring waits for the decision on the refreshed LTMP tables (the score uses manta observations).


## LTMP refresh (2026-10-16 / 19)

**Source:** COTS-Coral-OutcomesV2 `scripts/01_download_aims_ltmp.R` (client `R/aims_api.R`), snapshot `data/raw/aims-ltmp/2026-10-16`. It is identical, table for table, to the 2026-10-15 snapshot.

**Builder:** `domain_building/refresh_ltmp_tables.py` writes `sandbox/data/reef_cots.csv` and `reef_manta.csv` in the existing layout. The previous tables are archived (hash-checked) under `sandbox/data/archive_ltmp/0d8871dd90a5/` and `839d857a5c02/`. Report: `sandbox/data/derived/ltmp_refresh_2026-10-16/`.

**COTS surveys:**
- 2,863 existing surveys are unchanged (year, tows, counts);
- 297 added;
- 5 removed, all Hope Reef (Burdekin): 2 renamed to "Hope Reef (19046)", 3 dropped by AIMS (1989–91).

The additions are new Lizard-sector surveys from April and Sept–Oct 2026, plus back-filled histories for reefs newly in the AIMS list. The Sept–Oct 2026 surveys carry AIMS reporting year 2027:

| Reef | COTS/tow (late 2026) |
| --- | --- |
| Linnet | 0.51 |
| Martin | 0.35 |
| North Direction | 0.12 |
| MacGillivray | 0.06 |
| Helsdon | 0.06 |
| Lizard | 0.02 |
| Yonge, Carter, No Name | 0 |

**Coral (manta, AIMS modelled):** all estimates were re-fitted by AIMS. At the calibration and candidate reefs the changes are ≤ 0.006 cover; elsewhere up to 0.06. 21 rows were removed elsewhere, none at calibration or candidate reefs.

**Check:** re-scoring run 42 with the refreshed table reproduces run 43 exactly (run 20261019T_01). The score uses manta reporting years 1992–2024 (v3: to 2021), so no existing score changes. Reporting year 2027 lies beyond the model's end year (2026).


## Cull re-check on the `w` 2 layer: results (2026-10-19)

**Runs:** run 20261016T_10 (10 arms, no failures); re-score 20261019T_02; G2′ 03; audit `cull_w2_audit.json` (`audit_1991_cull_w2.py`). The no-cull arms replay run 20261016T_02 exactly (0.0).

| | No cull | Culling (σ 210 m) |
| --- | --- | --- |
| v3 (seed mean) | 11.36 | 11.90 (worse in all 5 seeds) |
| Diver error | 3.09 | 3.38 |
| G2′ `A½ 30` | 5/5 | 5/5 |
| Lizard 2024 / 2025 / 2026 adults/ha (targets 8.6 / 9.1 / 5.0) | 5.1 / 7.5 / 18.7 | 4.5 / 6.6 / 17.2 |
| Eyrie 2023 / 2025 (targets 45 / 12) | 3.2 / 6.4 | 2.4 / 5.2 |
| Unmet removals 2022–26: Lizard / Eyrie / domain | — | 15% / 69% / 61% |

**Finding.** The model's 2020s build-up is about two years late, so culling still fails:
- observed culls removed COTS from dense 2022–24 populations that the model does not yet hold;
- in the model, culling deepens the 2022–25 shortfall;
- the model's surge comes in 2025–26, after the observed population had been culled down, and culling there removes only about 1 COTS/ha/yr on Lizard.

Making culling work requires the build-up to arrive on time. Culling stays off in calibration.


## Preregistered extra reefs (2026-10-19, nineteenth round)

**Change.** Runner switch `OWEN_1991_EXTRA_REEFS=true` writes `trajectories_extra.csv` and `flows_extra.csv` for Martin (14-123), Linnet (14-126), Yonge (14-138), Carter (14-137) and No Name (14-139), all already domain reefs. Default off; existing outputs are unchanged. Report-only, not part of any score.

**Run.** Experiment `pit_q` with `OWEN_1991_ARM_FILTER='wc_b100_q20_'` (W + `h10`, `Q` 20, 5 seeds) on `corr_w2` (main) and `corr_w3` (sensitivity), with extra reefs on.

**Control.** `trajectories.csv`, `flows.csv` and `mechanisms.csv` must replay the same arms of runs 20261016T_02 / _03 exactly.

**Read-out** (`audit_1991_extra_reefs.py`), against LTMP manta by reporting year. Model detected CPUE uses q 0.015, linear and with the `A½ 30` low-density detection curve.
- **Martin and Linnet:** wave peak years and heights (observed 1996–99 and 2013, then 0.35 / 0.51 in late 2026, reporting year 2027, beyond the model's 2026 end); the 2022–26 level (observed 0.03–0.10).
- **Yonge, Carter and No Name:** model years with detected CPUE ≥ 0.22, max detected CPUE, and mean adults (observed max ≤ 0.01 COTS/tow over 35 years).

**Expected.** Martin and Linnet behave like the calibration reefs. The absence reefs are a real test: they receive the same external supply rules, so false outbreaks there would point to supply or habitat assumptions.


## Extra reefs: results (2026-10-19)

**Runs:**
- runs 20261019T_04 (`corr_w2`) and _05 (`corr_w3`): 5 arms each, no failures;
- audit [`extra_reefs_audit.json`](runs/20261019T_04_extra_reefs_w2_1991_v1/extra_reefs_audit.json) and `extra_reefs.png` (`audit_1991_extra_reefs.py`).

The standard outputs replay runs 20261016T_02 / _03 exactly (0.0), so the switch does not alter the model.

**Absence reefs (Yonge, Carter, No Name; observed ≤ 0.01 COTS/tow in 30–35 surveys since 1986):**

| | Yonge | Carter | No Name |
| --- | --- | --- | --- |
| Model 1990s peak, adults/ha (year) | 23 (2001) | 32 (2002) | 22 (2000) |
| Model 2014 peak, adults/ha | 20 | 19 | 16 |
| Max detected COTS/tow (`A½ 30` / linear) | 0.15 / 0.35 | 0.25 / 0.49 | 0.14 / 0.33 |
| 2022–26 adults/ha | 0.2–2.2 | 0.6–3.2 | 0.2–1.8 |

- **The model produces false outbreaks on all three in the 1990s and 2014 waves.** It matches their 2020s absence.
- Nothing in the model distinguishes them from outbreak reefs:
  - **settlement** (W, `w` 2) of 32–51 recruits/ha/yr, similar to Martin's 39–43 and 59–76% external;
  - **coral** 17–26%, similar.
- **Chlorophyll does not distinguish them either:** median wet-season p90, in-water scale 2003–26, is 0.35–0.47 µg/L, against North Direction 0.23, MacGillivray 0.30, Martin 0.55 and Linnet 0.67.
  - The anomaly form used in the model removes between-reef differences in any case.

**Outbreak reefs (Martin, Linnet):**
- Wave timing is about right: model 1990s peak 2000 / 1998 (observed 1998–99 / 1998); 2014 (observed 2013).
- Detected heights are low (`A½ 30` max 0.25 / 0.31 against observed 0.46–0.97).
- 2022–26: model 2–10 adults/ha (detected 0.003–0.04) against observed 0.03–0.10; then 0.35 / 0.51 in the late-2026 surveys (reporting year 2027, beyond the model's end).
- `w` 3 raises Martin and Linnet in 2023–26 by about 30–50%.

**Status.** Report-only; nothing is scored or adopted. The false outbreaks on the never-COTS reefs are a structural gap: the model has no reef-level property that keeps them free. The ecological explanation goes to the lead.


## Lead's decisions (2026-10-19, twentieth round)

1. **Never-COTS reefs.** Most likely a connectivity issue keeping COTS at bay, plus strong wave exposure on the eastern flank. Check the wave data (GBR ReefGuidance processed rasters: wave height and period) and the in-water logger chl out there, in case the satellite series obscures lower levels.
2. **Run the promotion checks,** noting that a way must still be found to bring the 2020s outbreak earlier. The lead suspects slow coral dominance (recovery of fast-growing coral) in this CoralBlox calibration and asks for that behaviour to be examined.


## Never-COTS reefs: data checks (2026-10-19)

Diagnostic `analyse_never_cots_reefs.py`; outputs in the git-ignored `sandbox/data/derived/never_cots_diagnostic/`.

| | Yonge / Carter / No Name | Outbreak reefs (Martin, Linnet, Lizard units, MacGillivray, North Direction, Eyrie) |
| --- | --- | --- |
| Wave power index (Hs² Tp) at the most exposed site | **15.0 / 14.1 / 16.1** (80–85th domain percentile) | 2.4–5.7 (5th–43rd) |
| Mean peak wave period Tp | **7.0–7.2 s** | 4.3–5.4 s |
| Mean Hs (reef, area-weighted) | 0.65–0.71 m | 0.49–0.80 m |
| In-water chl, nearest LTCMP station, median 1992–2009 | **0.13 µg/L** (Cod Hole, 2–9 km) | 0.21–0.23 µg/L |
| Satellite in-water-equivalent p90 (2003–26) | 0.35–0.47 | 0.19–0.67 |
| Larval arrivals per ha under equal production (domain percentile) | 0.11–0.16 (41st–55th) | 0.04–0.29 (4th–78th) |
| Arrivals from domain reefs / external | 80–85% / 12–15% | 53–85% / 15–43% |

(Linnet lies in a gap between the two wave rasters.)

**Findings.**
1. **Wave exposure separates them clearly.** Their exposed flanks receive ocean swell (Tp about 7 s), with wave power 3–6× that of the outbreak reefs' most exposed sites.
2. **In-water chl separates them clearly; the satellite series inverts it.**
   - The outer-shelf Cod Hole stations average 0.13 µg/L, about 40% below the mid-shelf Lizard-area stations, and below the Fabricius larval-development range.
   - The satellite-derived values rank these reefs higher, which is likely optically shallow water or bottom reflectance on the ribbon reefs.
3. **Reef-scale connectivity does not separate them.** Under equal production they receive mid-range larval arrivals, mostly from other domain reefs; Martin (an outbreak reef) receives fewer.
   - In the model their settlers come from simulated domain outbreaks, so a supply-side explanation would need lower survival of larvae travelling in low-chl outer-shelf water, not lower hydrodynamic connection.

**Implications (for the lead):**
- **Waves:** a post-settlement or habitat effect of wave exposure, such as reduced settlement or juvenile survival on exposed reef, would be a new reef-level mechanism, evidence-based on these rasters.
- **Chlorophyll:** the absolute larval-survival form fed by in-water-calibrated chl would suppress production from, and survival in, outer-shelf water. That needs an outer-shelf correction to the satellite series, since the current one inverts the cross-shelf gradient.


## Preregistered promotion checks (2026-10-19, twentieth round)

**Configuration under test (W-current):** W shares, `h10` start, `Q` 20, Ext×5, new coral (β 1, ×0.25 below 6 DHW, μ 0.12), chl `aprmak` (b 3), `owen_global_2018_2023_corr_w2`, culling off, q 0.015, score v3.

**Gate 5: ±25% sensitivity, one parameter at a time, 5 seeds.**
- **Parameters:**
  - Ext 3.75 / 6.25;
  - chronic μ 0.09 / 0.15;
  - low-DHW bleaching scale 0.1875 / 0.3125;
  - chl slope b 2.25 / 3.75 (new tables from `build_larval_survival_index_corrected.py --b`, whose default reproduces the b 3 table exactly);
  - initial fast shares ×0.75 / ×1.25 (capped at 0.9);
  - `Q` 15 / 25, already run (20261016T_02).
- **Runner experiment `sens`:** 55 arms, including 5 base arms that must replay 20261016T_02 `Q` 20 exactly.
- **Pass rule:** each perturbation passes G2′ (`A½ 30`, trough on Lizard and MacGillivray) in ≥ 4/5 seeds, and its seed-mean v3 is no more than 1.0 above the base (11.36). The result is "not confined to a narrow region" if every perturbation passes. Failures are reported with the qualitative change.

**Gates 4a / 4b: held-out years and leave-one-reef-out (re-scoring existing runs; no new model runs).**
- **Held-out years:** for each screening decision in the current chain, recompute the loss on the fit window only (manta 1992–2012) and the check window (manta 2013–2021 plus diver 2020–26). The decisions are: starting state (runs 36 / 37), correction weight × `Q` (runs 20261016T_01–03), and Ext × local output (run 45).
  - **Pass:** the chosen option ranks first, or within the seed spread of first, on the fit window, and does not rank last on the check window.
- **Leave-one-reef-out:** the same decisions re-ranked with each calibration reef dropped from the summed v3.
  - **Pass:** the chosen option stays first or within the seed spread for every dropped reef.

**Coral dominance (lead's hypothesis).** Diagnostic, in parallel: model fast-coral share and total cover at the calibration reefs 2015–26 against LTMP photo-transect composition (Acroporidae and Pocilloporidae share) from the 2026-10-16 snapshot.


## Promotion checks: results (2026-10-19)

### Gate 5: ±25% sensitivity

Run 20261019T_06 (55 arms, no failures); re-score 10, G2′ 11; audit [`sens_audit.json`](runs/20261019T_06_sens_w2_1991_v1/sens_audit.json) and `sens.png` (`audit_1991_sens.py`). The base arms replay run 20261016T_02 exactly (0.0). Base v3 is 11.36.

| Parameter | ×0.75: Δv3, G2′ | ×1.25: Δv3, G2′ | Notes |
| --- | --- | --- | --- |
| Pit capacity `Q` (15 / 25) | +0.12, 5/5 ✓ | −0.17, 5/5 ✓ | Lower `Q` gives earlier, larger 2020s build-up |
| External multiplier (3.75 / 6.25) | +0.81, 4/5 ✓ | **+1.28, 5/5 ✗** | ×1.25: one seed collapses (v3 17.5); 2026 overshoot grows |
| Chronic coral mortality μ (0.09 / 0.15) | +0.75, **3/5 ✗** | **−0.50**, 4/5 ✓ | Higher μ fits better (manta part 6.73 → 6.22) |
| Low-DHW bleaching scale | 0.00, 5/5 ✓ | 0.00, 5/5 ✓ | Not identified by these data |
| Chl slope b (2.25 / 3.75) | +0.03, 5/5 ✓ | 0.00, 5/5 ✓ | Changes the 2026 level only |
| Initial fast shares | −0.10, 5/5 ✓ | +0.73, 5/5 ✓ | |

**Gate 5 fails on 2 of 12 perturbations** (Ext ×1.25 on v3; μ ×0.75 on G2′). The result is not confined to a knife edge, but it is sensitive to more external supply and to lower background coral mortality. Peak years in the 1990s move by at most about 0.8 y across all perturbations.

### Gates 4a / 4b: held-out years and leave-one-reef-out

`rescore_holdout.jl` re-scores with the frozen functions on truncated windows. Its full-window self-check reproduces the stored v3 manta part exactly (0.0; run 20261019T_07). Held-out scores: run 08. Gate audit: run 09 (`audit_promotion_gates.py`).

"Within the seed spread of first" is made operational, after the preregistration, as: the paired seed difference from the first-ranked option is positive in fewer than 4 of 5 seeds.

| Decision (chosen) | 4a: fit 1992–2012 rank | 4a: check 2013–26 rank | 4b: rank with each reef dropped | Verdict |
| --- | --- | --- | --- | --- |
| Ext ×5, local ×1 (`e50_f100`) | 2 of 6 (within spread) | 2 of 6 | 1 for every reef | **pass** |
| Start `h10` (both layers) | 2 of 3 (within spread) | 2 of 3 | 1, except without North Direction (2 of 3; `lin` better in 4/5 seeds) | **marginal**: its advantage is mostly North Direction's earlier 1990s wave |
| Correction weight `w` 2, `Q` 20 | 4 of 12 (`Q` 25 better in 4/5 seeds) | 8 of 12 | 4–6 of 12; `w` 3 `Q` 25 (or `Q` 15) first | **fail** |

The fit window ranks only `Q` (`w` acts from 2018), and favours `Q` 25. On the full v3 and every leave-one-out, `w` 3 with `Q` 25 (v3 11.01) ranks first. The earlier selection rule picked `w` 2 `Q` 20 as the "smallest change", not as the best, and these gates do not support that preference.

### Coral dominance (lead's hypothesis)

Diagnostic `analyse_coral_dominance.py` (outputs in run 20261016T_02: `coral_dominance.*`). It compares the model fast share (fast / total coral) with LTMP 9 m photo-transect (*Acropora* + Pocilloporidae) / hard coral.

| Reef | Model fast share 2019 / 2021 / 2023 | Observed 2019 / 2021 / 2023 | Model total cover 2021 | Manta HC 2021 |
| --- | --- | --- | --- | --- |
| Lizard | 0.25 / 0.55 / 0.65 | 0.13 / 0.23 / 0.45 | 0.22 | 0.13 |
| MacGillivray | 0.26 / 0.59 / 0.69 | 0.11 / 0.20 / 0.30 | 0.23 | 0.19 |
| North Direction | 0.27 / 0.30 / 0.47 | 0.27 / 0.58 / 0.59 | 0.15 | **0.50** |

(Model cover is a fraction of habitable area; manta cover is a fraction of reef surface.)

**Finding.** Fast-coral dominance is **not** slow at Lizard or MacGillivray: the model's fast share recovers faster and higher than observed, and food survival is 1.0 from 2021. North Direction's recovery is much too slow in the model (observed manta HC 0.50 and fast share 0.58 in 2021). Before 2016 the model strips fast coral harder than observed (2013 share 0.02–0.10 against about 0.32). So slow coral recovery does not explain the late 2020s build-up at Lizard; it may contribute at North Direction.

### Status

W-current is **not promoted**. It fails gate 4 on the correction-weight × `Q` choice and gate 5 on two perturbations, with the late 2020s build-up and the never-COTS false outbreaks as known structural gaps. Decisions go to the lead.


## Coral dominance: correction (2026-10-19)

**Error.** The coral-dominance section above estimated observed fast cover as the photo-transect fast share × *manta* hard coral, and stated the API photo data were group-level only. That is wrong: `HARD CORAL` rows with `purpose == "COMPOSITION"` carry taxa in `reefpage_category` (the lead had pointed this out before). Manta (reef-wide, shallow) and 9 m photo transects differ greatly at North Direction: hard coral 0.50 against 0.11 in 2021.

**Fix:**
- `refresh_ltmp_tables.py --photo` refreshed `sandbox/data/reef_photo_transect.csv` from the 2026-10-16 snapshot. It adds mainly 2026; the previous copy is archived; revisions at our reefs are ≤ 0.024 cover.
- `analyse_coral_dominance.py` (version 2) now compares like with like: photo hard coral (`GROUP_LEVEL`) and photo *Acropora* + Pocilloporidae (`COMPOSITION`).
- Outputs `coral_dominance_photo.*` in run 20261016T_02; version 1 outputs `coral_dominance.*` are kept and superseded.

| Reef | Model fast cover 2019 / 2021 / 2023 / 2025 | Photo fast cover 2019 / 2021 / 2023 / 2025 | Model / photo fast share 2021 |
| --- | --- | --- | --- |
| Lizard | 0.026 / 0.12 / 0.21 / 0.17 | 0.008 / 0.023 / 0.093 / 0.057 | 0.55 / 0.22 |
| MacGillivray | 0.027 / 0.14 / 0.23 / 0.19 | 0.011 / 0.027 / 0.053 / 0.064 | 0.59 / 0.20 |
| North Direction | 0.034 / 0.044 / 0.093 / 0.14 | 0.017 / 0.061 / 0.069 / 0.073 | 0.30 / 0.56 |

(Model cover is a fraction of habitable area; photo cover is a fraction of the 9 m transect.)

**Corrected finding.** There is no evidence of slow fast-coral recovery in the model:
- at Lizard and MacGillivray model fast cover recovers earlier and 2–5× higher than observed;
- at North Direction it is close (0.04–0.12 against 0.06–0.11 in 2021–24); the previous "far too slow" statement is withdrawn.

Before 2016 the model strips fast coral harder than observed (2013 share 0.02–0.10 against 0.12–0.32). Food therefore does not explain the late 2020s build-up.


## Lead's decisions (2026-10-19, twenty-first round) and the ribbon connectivity check

**Decisions (all agreed):**
1. Adopt correction weight `w` 3 with `Q` 25 as the working configuration (W + `h10`), and rerun the ±25% sensitivity on it.
2. Model the never-COTS reefs first with wave exposure: an opt-in reef-level effect reducing settlement or juvenile survival where exposure is high, driven by the GBR ReefGuidance wave rasters.
3. Diagnose the 2016–22 juvenile-to-adult pipeline at Lizard and Eyrie against observed size classes.

Work through round 20 was committed (`df024da1`) before this check.

**Ribbon connectivity check** (lead: the never-COTS reefs draw most larvae from within the domain; are they connected to the reefs that actually had COTS?). Diagnostic `analyse_ribbon_connectivity.py`; outputs in the git-ignored `sandbox/data/derived/ribbon_connectivity/`.
- **The internal COTS connectivity is hydrodynamic, not the distance-based stand-in.** The dataset's `reef_connectivity_<year>.csv` equal the Owen global sub-matrix (max difference 1e-16), sampled by hydrodynamic year. The distance-based `connectivity/Lizard_Connectivity.csv` is ADRIA's coral connectivity.
- **Outbreak-weighted larval pressure,** `Σ_j C[j,i] A_j/A_i s_j × (P_j, LTMP outbreak, or culled ≥ 1/ha/yr)`, relative to Martin (an outbreak reef):

  | Reef | Predicted P (W3 / W4 / W5) | LTMP outbreak reefs (W3 / W4) | Culled reefs (W4 / W5) |
  | --- | --- | --- | --- |
  | Yonge / Carter / No Name | 2.3–3.0 / 2.4–3.1 / 3.3–4.5 | 1.6–1.8 / 1.4–1.5 | 1.6–1.8 / 2.7–3.1 |
  | Linnet / Eyrie | 1.6–1.7 / 1.5–1.7 / 1.5–1.7 | 1.4–1.6 / 1.4–1.8 | 1.5–1.9 / 1.3–2.0 |
  | Lizard / MacGillivray / North Direction | 4.1–7.6 / 4.2–7.9 / 3.7–6.7 | 4.4–8.2 / 4.1–7.0 | 4.9–8.4 / 5.9–11.1 |

- **Findings:**
  1. The ribbon reefs are **not** poorly connected to reefs known to have COTS. They receive more outbreak-weighted larvae than Martin, Linnet or Eyrie (all outbreak reefs), including from Mackay Reefs (15-024, LTMP outbreak in W3) and U/N 15-022 / 14-149 (P 0.5–0.7).
  2. 57–59% of their larval connection weight comes from **Ribbon No 10 Reef (14-146)**, a domain reef with P about 0.015 and no outbreak record. About 96% comes from domain reefs, mostly the outer ribbon cluster.
  3. In the model the ribbon cluster can therefore seed itself: false outbreaks on Ribbon No 10 would feed Yonge, Carter and No Name. A reef-level effect on the whole outer ribbon cluster (wave exposure, low in-water chl) remains the best-supported explanation; connectivity alone is not.
  - Next: add Ribbon No 10 (14-146) to the extra reported reefs.

**Next steps after compaction (in order):**
1. **Working configuration W-next:** W + `h10`, `Q` 25, `owen_global_2018_2023_corr_w3`, Ext×5, new coral, chl `aprmak`, culling off.
   - Rerun gate 5 (±25%) on it. The runner's `sens` experiment hard-codes `q20` in its arm names and was run on `corr_w2`, so parameterise its base `Q`. `Q` ±25% of 25 is 18.75 / 31.25; the arm-name parse is integer, so either add a decimal token or use 19 / 31 and record it.
   - Extend the extra reefs with Ribbon No 10.
2. **Wave-exposure mechanism** (preregister before coding).
   - **Hypothesis:** high ocean-swell exposure lowers COTS settlement or juvenile survival.
   - **Implementation:** opt-in ADRIA adapter hook, like `COTS_LARVAL_SURVIVAL_CSV`: a per-site multiplier on settling recruits from all sources, applied after dispersal and external supply.
   - **Input:** built by a checked-in script from `never_cots_diagnostic/site_waves.csv`; fill Linnet's raster gap by nearest valid site and record it.
   - **Choices for the lead:** site-level versus reef-level exposure (reef mean Hs does not separate the ribbons; most-exposed-site power and Tp do), and the functional form and floor.
   - **Controls:** the switch off replays W-next exactly; test the extra reefs plus the calibration score. COTSMod is unchanged unless the hook needs it.
3. **2016–22 juvenile-to-adult pipeline diagnostic** at Lizard and Eyrie: model recruits, subadults and adults by year against observed size classes (cull cohorts 1–4 in `cull_reef_year.csv`; SALAD sizes if available).
4. **Later:** outer-shelf in-water correction of the satellite chl (absolute form); culling once the build-up is on time; a cyclone-track tool.


## Preregistered: gate 5 on W-next (2026-10-19, twenty-first round, item 1)

**Configuration under test (W-next):** W shares, `h10` start, `Q` 25, Ext×5, new coral (β 1, ×0.25 below 6 DHW, μ 0.12), chl `aprmak` (b 3), `owen_global_2018_2023_corr_w3`, culling off, q 0.015, score v3.

**Runner change.** `OWEN_1991_SENS_Q` sets the base `Q` of experiment `sens` (default 20 reproduces run 20261019T_06's arm set). Any other base adds `Q` × f through the token `_sxq<f>_`, which multiplies the parsed `Q`, so 25 × 0.75 / 1.25 = 18.75 / 31.25 exactly (no rounding). Ribbon No 10 Reef (14-146; no LTMP surveys) is added to the opt-in extra reported reefs (`OWEN_1991_EXTRA_REEFS=true`), which writes `trajectories_extra.csv` / `flows_extra.csv` only.

**Design.** 65 arms (5 seeds): 5 base arms, which must replay the `q25` arms of 20261016T_03 (`pit_q`, corr_w3) exactly, and ×0.75 / ×1.25 on Ext, μ, low-DHW scale, chl slope b, `Q` and the initial fast shares.

**Pass rule (as for W-current).** Each perturbation passes G2′ (`A½ 30`, trough on Lizard and MacGillivray) in ≥ 4/5 seeds, and its seed-mean v3 is no more than 1.0 above the base (11.01 on run 20261016T_03). The result is "not confined to a narrow region" if every perturbation passes. Failures are reported with their qualitative change.

**Also reported (not scored):** Ribbon No 10 against Yonge / Carter / No Name. If Ribbon No 10 has false outbreaks that lead the other ribbon reefs, the ribbon cluster self-seeds in the model.


## Results (2026-10-19): gate 5 on W-next, Ribbon No 10, wave evidence, juvenile pipeline

### Gate 5 (±25%) on W-next

Run 20261019T_12 (65 arms, no failures); re-score 13 (control replay 0.0); G2′ 14 (trough on Lizard and MacGillivray); audit `sens_audit.json` / `sens.png` in run 12 (`audit_1991_sens.py --base-q 25`). Base arms replay the `q25` arms of 20261016T_03 exactly (0.0); base v3 is 11.01.

| Parameter | ×0.75: Δv3, G2′ | ×1.25: Δv3, G2′ |
| --- | --- | --- |
| Pit capacity `Q` (18.75 / 31.25) | +0.30, 5/5 ✓ | +0.01, 5/5 ✓ |
| External multiplier (3.75 / 6.25) | **+1.55, 3/5 ✗** | +0.57, 5/5 ✓ |
| Chronic coral mortality μ (0.09 / 0.15) | +0.58, **3/5 ✗** | −0.17, 4/5 ✓ |
| Low-DHW bleaching scale | −0.08, 5/5 ✓ | +0.10, 5/5 ✓ |
| Chl slope b (2.25 / 3.75) | +0.02, 5/5 ✓ | 0.00, 5/5 ✓ |
| Initial fast shares | −0.11, 5/5 ✓ | +0.60, 5/5 ✓ |

**Gate 5 fails on 2 of 12** (Ext ×0.75; μ ×0.75 on G2′). Compared with W-current, the Ext failure has moved from more supply to less. W-next depends on external supply at about the ×5 level, and on background coral mortality of at least about 0.12. 1990s peak years move by at most 0.8 y. W-next is **not promoted**.

### Ribbon No 10 (extra reefs, base arms)

Ribbon No 10 does **not** lead the other ribbon reefs. It peaks at 16 adults/ha in 2002 and 7.8 in 2015, against Yonge / Carter / No Name at 21–30 in 2001–03 and 14–18 in 2014. Settlement on the ribbon reefs in 1993–2016 is 62–65% external supply (from non-domain reefs) and 33–37% from domain reefs. Ribbon No 10 is 64% external.

So the false ribbon outbreaks are fed mainly by external supply, not by the ribbon cluster seeding itself. A sink-side (settlement) effect reaches all sources; a source-side effect on domain reefs would reach only about a third.

### Wave exposure: independent check against all LTMP reefs in the domain

Mean wave period Tp, area-weighted over each reef (GBR ReefGuidance rasters), against LTMP manta history 1986–2027 for the 29 domain reefs with wave data:
- **Tp ≥ 6.98 s (9 reefs):** Carter, No Name, Yonge, Harrier, Ribbon Reefs 1 / 3 / 6 / 9 and Lena never reached 0.22 COTS/tow (max ≤ 0.06). The exception is Irene (Tp 7.62; max 1.13, 2 years ≥ 0.22).
- **Tp ≤ 6.59 s (20 reefs):** 18 had at least one year ≥ 0.22 (all except Fly and Hilder). They include Marx (6.59), Forrester (6.03) and Mackay (5.72).

Most of these reefs are outside the score, so this is evidence beyond the calibration fit. Tp marks ocean-swell exposure and so also outer-shelf position. It does not separate waves from other outer-shelf properties such as low in-water chl.

In the domain, 52 of 113 reefs have Tp ≥ 6.0, 30 have Tp ≥ 7.0, and 27% of reef area has Tp ≥ 6.5. Martin, Linnet and Maxwell have 29 sites missing raster values; these take the nearest valid site's values (≤ 4.2 km). Martin is 5.33 s and Linnet 5.25 s.

### Wave hook (implemented, off by default)

- **ADRIA `COTS_SETTLEMENT_MODIFIER_CSV`** (`site_id, multiplier` in [0, 1]; every site exactly once). After the COTS transition and before dispersal, it multiplies each site's settlement gate, and so scales local retention, internal immigration, external supply and pulses exactly once. It also scales the background immigration already settled in the transition (`N[1]` and its flux). It requires separated recruitment.
- **Unset is the established path.** Test `test_settlement_modifier.jl` passes 15/15. COTSMod is unchanged.
- **Builder `domain_building/build_wave_settlement_modifier.py`** (`--metric tp_reef | power_reef_max | power_site`, linear ramp from 1 at `e0` to `floor` at `e1`). Not yet run; it awaits the lead's choice.

### Juvenile-to-adult pipeline, 2010–26

Diagnostic `analyse_juvenile_pipeline.py` (W-next arms of 20261016T_03, 5 seeds); outputs in the git-ignored `sandbox/data/derived/juvenile_pipeline/` (summary statistics only). Maturation yield = adults matured in year t / subadults at the end of year t−1.

| Reef | Subadults 2015–17 (/ha) | Maturation yield 2012–14 → 2015–18 → 2020–24 | Pit removals 2015–18 (total /ha) | Min food survival 2015–19 | First year model adults ≥ 5/ha |
| --- | --- | --- | --- | --- | --- |
| Lizard | 88 | 0.20 → **0.025** → 0.12 | 68 | 0.06 | 2024 |
| Eyrie | 50 | 0.31 → **0.076** → 0.11 | 56 | 0.06 | 2024 |
| MacGillivray | 88 | 0.31 → **0.024** → 0.09 | 68 | 0.05 | 2026 |
| North Direction | 132 | 0.27 → **0.028** → 0.12 | 77 | 0.05 | 2019 |

**Observed:**
- Eyrie: 531 culled in 2021, 71% over 25 cm.
- Lizard: culls in 2021–23 are 75–91% over 25 cm, and the SALAD median diameter is 350 mm in 2021 and 2022 (5% under 150 mm).
- North Direction: the SALAD median is 300–500 mm in 2022–25.

With approximate growth (about 15–25 cm at age 2, over 25 cm at age 3+), the adults present in 2021 settled in about 2018 or earlier, during or before the coral trough. Small COTS are under-detected by divers (detection below 0.5 under 15 cm), so cohort 1 shares are not used.

**Finding.** The model settles a large 2012–16 cohort, then loses it at the juvenile-to-adult step in 2015–18. The predator pit removes about 55–77/ha. The remainder starves at maturation, because food survival is about 0.05 while fast coral is near zero after the 2016–17 bleaching. The model's 2020s adults therefore have to come from new settlement after 2019, which is low (Lizard 7–16/ha/yr in 2019–21) and arrives about 3 years late.

Deaker et al. 2020 show juveniles can stay herbivorous on coralline algae for years and resume growth once coral returns, the "juvenile bank". In this model, juveniles that cannot find coral die instead of waiting.

The model already has this mechanism behind a switch (COTSMod `preferred_maturation`: maturation fraction √F̄, retained juveniles exempt from starvation). It failed on 2026-10-12 because fast coral on the old coral calibration never recovered enough to release the bank. On the new coral, fast coral recovers by 2019–21, so the test is no longer the same.

### Chlorophyll (scoping, not built)

The current chl index (`aprmak`) averages 1.0 on every reef: it carries year-to-year anomalies only, with no cross-shelf level. Owen source survival is off (`COTS_APPLY_LARVAL_SURVIVAL=false`).

In-water LTCMP medians for 1992–2009 are:
- outer shelf: Cod Hole 0.13 / 0.13, Osprey 0.12–0.13, Cairns shelf-break 0.15, Saxon 0.15 µg/L;
- mid shelf, Lizard area: 0.21–0.27 µg/L;
- inner-shelf river influence: Port Douglas and Daintree 0.35–0.62 µg/L.

A static reef-level factor `(chl_iw / chl_ref)^b`, capped at 1, applied at the source through the existing `COTS_LARVAL_SURVIVAL_CSV` path (multiplied into `aprmak`), would lower outer-shelf larval production. Ecologically, larvae spawned in oligotrophic outer-shelf water develop poorly. It needs a per-reef in-water climatology (interpolated from LTCMP stations by cross-shelf position) and a decision on whether it also scales the external boundary sources. Because about two thirds of ribbon settlement is external, chl at the domain's sources alone would do little for the ribbons. Combined with waves, it would mainly cut outer-shelf production that leaks to mid-shelf reefs.

### Status

W-next is not promoted (gate 5 fails 2/12). Decisions requested from the lead: (1) the wave-mechanism form; (2) whether to retest the juvenile bank on the new coral; (3) whether to build the cross-shelf chl term to combine with waves.


## Lead's decisions (2026-10-19, twenty-second round)

1. **Wave mechanism agreed as proposed:** reef-average wave period, a ramp from 6.6 to 7.0 s, and floors 0.1 / 0.25 / 0.5.
2. **Juvenile bank agreed:** retest `preferred_maturation` on W-next, with the pit on.
3. **Cross-shelf in-water chl term agreed, after the wave test.** Also test whether to reduce supply only from domain sources (inside) or also from the external boundary sources (outside).

Order: waves → juvenile bank → chl (inside versus inside + outside).


## Preregistered: wave-exposure settlement mechanism (2026-10-19, twenty-second round)

1. **Hypothesis.** On reefs exposed to ocean swell, COTS settlement and early post-settlement survival are low (dislodgement, unstable habitat and fewer sheltered juvenile microhabitats). These reefs therefore stay below outbreak densities even with high larval supply. The evidence is the LTMP cross-section in the domain: mean wave period ≥ 7 s → 8 of 9 reefs never reached 0.22 COTS/tow, against 18 of 20 outbreak reefs at ≤ 6.6 s.

2. **Change.** Each site's settlement gate is multiplied by `m = 1 − (1 − floor) · clamp((Tp_reef − 6.6) / (7.0 − 6.6), 0, 1)`, where `Tp_reef` is the reef's area-weighted mean peak wave period (s).
   - The multiplier is time-invariant and shared by all sites of a reef.
   - It acts after the COTS transition and before dispersal, so it scales settlement from all sources (local, internal, external, background) exactly once.
   - Switch: `COTS_SETTLEMENT_MODIFIER_CSV`; unset is the established path.
   - Inputs: `build_wave_settlement_modifier.py --metric tp_reef --e0 6.6 --e1 7.0 --floor {0.1, 0.25, 0.5}`, giving `sandbox/data/derived/wave_settlement_modifier/tp66_70_f{010,025,050}/`. Twenty-nine sites without raster values (Martin, Linnet, Maxwell) take the nearest valid site's values (≤ 4.2 km).
   - Affected: 37 of 113 reefs, including Yonge, Carter, No Name, Ribbon No 10, the other ribbon reefs and Irene. Unaffected (m = 1): all score reefs, Martin, Linnet, Marx, Forrester, Mackay and 15-022.

3. **Expected consequences.**
   - Ribbon reefs: false 1990s, 2000s and 2014 outbreaks are much reduced. With floor 0.25 and about 63% external settlement, settlers fall to about a quarter; through the pit and the Allee effect, adults should fall by more than that.
   - Score reefs: affected only through less larval production on exposed reefs reaching them (Lizard receives about 8% of its connection weight from Ribbon No 10). Expect a small reduction in 1990s and 2014 peak heights, and little change in timing.
   - The late 2020s build-up is not expected to be fixed by this.

4. **Control and factorial.** W-next (`wc_b100_q25_eb_h10_obs015_s<k>_t60_e50_aprmak_cv2_pnd`, `corr_w3`, Ext×5, culling off) × {off, floor 0.5, 0.25, 0.1} × 5 seed pairs = 20 arms. Runner experiment `wave`; arms carry a suffix `_wv<floor×100>` (none = off).
   - The off arms must replay the `q25` arms of 20261016T_03 exactly.
   - Extra reefs are on (Martin, Linnet, Yonge, Carter, No Name, Ribbon No 10).

5. **Parameters and evidence.** The ramp ends (6.6 and 7.0 s) come from the LTMP cross-section above (largest outbreak reef Tp 6.59; smallest never-COTS reef Tp 6.98). The floor has no field estimate; 0.1–0.5 are labelled sensitivities.

6. **Readouts and failure criteria.**
   - **Ribbon reefs** (Yonge, Carter, No Name): seed-mean model adults and detected CPUE (q 0.015 with `A½ 30`) against LTMP (never ≥ 0.22).
     - **Success:** detected CPUE stays below 0.22 in every year 1992–2026, in at least 4 of 5 seeds.
   - **Calibration:** v3 and G2′ (`A½ 30`, trough on Lizard and MacGillivray).
     - **Failure:** seed-mean v3 rises by more than 0.5, or G2′ passes in fewer than 4/5 seeds.
   - Martin and Linnet: must keep their 1990s and 2014 outbreaks.
   - Fluxes: settled recruits and their external / internal split on the extra reefs.
   - **Rollback:** the switch is unset by default.


## Wave-exposure screen: results (2026-10-19)

Run 20261019T_15 (20 arms, no failures); re-score 16 (control replay 0.0); G2′ 17 (trough on Lizard and MacGillivray); audit `wave_audit.json`, `wave_summary.csv` and `wave.png` in run 15 (`audit_1991_wave.py`). The off arms replay the `q25` arms of 20261016T_03 exactly (0.0).

| Floor | Seed-mean v3 (Δ) | G2′ | Max detected COTS/tow (seed mean), Yonge / Carter / No Name | Seeds with no year ≥ 0.22 (Y / C / N) | Ribbon settlers 1993–2016 (/ha/yr) |
| --- | --- | --- | --- | --- | --- |
| off | 11.01 | 5/5 | 0.18 / 0.27 / 0.16 | 4 / 1 / 5 | 38.6 |
| 0.5 | 10.90 (−0.11) | 5/5 | 0.05 / 0.10 / 0.04 | 5 / 5 / 5 | 17.7 |
| 0.25 | 11.01 (+0.00) | 5/5 | 0.006 / 0.018 / 0.003 | 5 / 5 / 5 | 8.6 |
| 0.1 | 10.89 (−0.12) | 5/5 | 0.002 / 0.003 / 0.001 | 5 / 5 / 5 | 3.6 |

(LTMP maximum on these reefs is 0.00–0.01.)

**Findings.**
- **Ribbon reefs:** every floor meets the preregistered success criterion; floor 0.25 is closest to the observed levels.
- **Calibration:** unchanged (Δv3 within ±0.12; G2′ 5/5). Score-reef peak years are identical in every wave. Peak heights change by ≤ 1 adult/ha, except North Direction's 2013 peak, which rises from 60 to 65–69.
- **Martin and Linnet:** unchanged. The preregistered "keep the 1990s and 2014 outbreaks" check fails in every arm, including off: model Martin reaches 0.15 detected COTS/tow in 2014 against an observed 0.51. That is an existing under-prediction, not a wave effect.
- **The 2020s are not changed**, as expected.

**Status.** The mechanism passes its screen; the floor is the lead's choice (0.25 recommended). Not promoted until the combined configuration passes the gates.


## Preregistered: juvenile bank retest on the new coral (2026-10-19, twenty-second round)

1. **Hypothesis.** Juveniles that cannot find preferred coral stay herbivorous on coralline algae and wait, instead of starving at maturation (Deaker et al. 2020). The 2014–16 cohort would then survive the 2016–17 coral trough and mature as fast coral returns in 2019–21, bringing the 2020s adults forward by about 2–3 years.

2. **Change.** This uses the existing COTSMod switch `preferred_maturation` (`COTS_PREFERRED_MATURATION=true`; runner treatment code `PDA` instead of `PA`), with no code change.
   - The maturing fraction is `g = √F̄`, where F̄ is remembered preferred (fast) cover (habitable fraction).
   - Matured juveniles face the food factor `f`; retained juveniles `(1 − g)` stay in N2, exempt from starvation and still subject to `m2` and the predator pit each year.
   - The 2026-10-12 test of this switch failed on the old coral, where F̄ never recovered. This is a retest on the new coral, not a new parameter.

3. **Expected consequences.**
   - **Benefit:** 2015–18 maturation yield stops collapsing to about 0.03, and adults at Lizard and Eyrie reach 5/ha by about 2021–22 instead of 2024.
   - **Risk 1:** g is about 0.25–0.45 at typical F̄ (0.06–0.2), so maturation is slower in normal years and the 1990s and 2014 waves may be delayed or lowered, as in the earlier test.
   - **Risk 2:** the pit erodes the bank each year, which may cancel the benefit.

4. **Control and factorial.** W-next × bank {off, on} × wave {off, floor 0.25} × 5 seed pairs = 20 arms.
   - Runner experiment `bank`; arms `wc_b100_q25_eb_h10_obs015_s<k>_t60_e50[_jb]_aprmak_cv2_pnd[_wv025]`.
   - The bank-off arms must replay run 20261019T_15 (wave off: also 20261016T_03) exactly.

5. **Evidence and bounds.** The form and its √ follow CoCoNet; the retention-without-starvation evidence is Deaker et al. 2020 (6.5 years on coralline algae). No new parameters.

6. **Readouts and failure criteria** (paired by seed and wave setting).
   - **Required improvement:** a lower diver error (2020–26).
   - **Failure:**
     - seed-mean v3 rises by more than 0.5, or G2′ passes in fewer than 4/5 seeds;
     - a 1990s peak (Lizard, MacGillivray, North Direction) moves later by more than 1 year;
     - the 2020s are not earlier: the first year with model adults ≥ 5/ha at Lizard and Eyrie does not move earlier.
   - **Also reported:** maturation yield 2015–18, bank size (subadults) and pit removals, and the ribbon reefs.
   - **Rollback:** the switch is off by default.
   - If the 1990s are delayed, the next form to propose is a bank that operates only when coral is scarce. That would be a new preregistration.


## Juvenile bank retest: results (2026-10-19)

Run 20261019T_18 (20 arms, no failures); re-score 19 (control replay 0.0); G2′ 20; audit `bank_audit.json`, `bank_summary.csv` and `bank.png` in run 18 (`audit_1991_bank.py`). The bank-off arms replay run 15 exactly (0.0).

| Wave | Bank | v3 | Manta part | Diver error | G2′ | 1990s peak year (Lizard / MacGillivray / North Direction) | First year adults ≥ 5/ha (Lizard / Eyrie) | Lizard subadults 2018 | Lizard adults 2022 / 2026 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| off | off | 11.01 | 6.48 | 3.02 | 5/5 | 1997.0 / 1997.4 / 1997.0 | 2023.2 / 2024.2 | 21 | 3.5 / 23.7 |
| off | on | **16.50** | 9.60 | **4.60** | **0/5** | **1999.8 / 2000.2 / 1999.8** | 2019.0 / 2020.4 | 94 | 11.1 / 10.4 |
| 0.25 | off | 11.01 | 6.32 | 3.13 | 5/5 | 1997.0 / 1997.4 / 1997.0 | 2024.2 / 2024.2 | 22 | 3.2 / 23.6 |
| 0.25 | on | **16.51** | 9.59 | **4.62** | **0/5** | 1999.8 / 2000.2 / 1999.8 | 2019.0 / 2020.4 | 94 | 11.1 / 10.3 |

**Verdict: fail** on every preregistered criterion except "earlier 2020s":
- Δv3 is +5.5 and G2′ passes 0/5;
- the 1990s peaks are 2.8 y later, and the 2014 peaks are lower;
- the diver error is worse (+1.5 to +1.6).

**What happened.**
- **The bank works as intended through the trough.** About 94 subadults/ha survive to 2018 at Lizard, against 21, and are released as fast coral returns.
- **Timing and size are wrong.** Adults reach 15–20/ha in 2020–21 at Lizard, MacGillivray and North Direction, against observed 3–8/ha, then decline. That is too early and too high; Eyrie stays too low (about 5/ha against 16–40).
- **The 1990s delay comes from the form.** With `g = √F̄` and F̄ ≈ 0.05–0.2, only 22–45% of juveniles mature in an ordinary year.

Waves do not interact with the bank. The test shows that juvenile survival through the trough can move the 2020s, but this form is too strong in normal years and releases too much too early.


## Preregistered: scarcity-only juvenile bank (2026-10-19, twenty-second round)

Agreed with the lead as the follow-up form if the retest delayed the 1990s.

1. **Hypothesis.** Juveniles wait on coralline algae only when they cannot find coral; when coral is available they mature as normal. The bank therefore forms in coral troughs and costs nothing in normal years.

2. **Change** (new opt-in COTSMod switch `starvation_bank`; ADRIA env `COTS_JUVENILE_STARVATION_BANK`, default false). Of the background and pit survivors J:
   - a fraction `f` (the existing food-survival factor) matures;
   - the fraction `1 − f`, which currently starves, stays in N2 as herbivorous juveniles (subject again to `m2` and the pit next year);
   - nothing starves.
   - Off reproduces the current path exactly. It is exclusive with `preferred_maturation` and `juvenile_storage`. No new parameter; the logged maturation fraction is `f`.

3. **Expected consequences.**
   - **1990s:** f is close to 1 in normal years (0.85–1.0 at Lizard in 2010–13), so 1990s timing should be nearly unchanged.
   - **Trough:** in 2015–18 (f about 0.05) about 95% of juvenile survivors wait. The bank is eroded by `m2` (about 0.32/yr) and the pit, and is released as f recovers in 2019–20.
   - **Risk:** release is still in 2019–21, possibly earlier than the observed rise (Lizard 2022; Eyrie 2021–23), and the size of the release is not tuned.

4. **Control and factorial.** W-next × {scarcity bank off, on} × wave {off, 0.25} × 5 seeds = 20 arms. Runner experiment `bank2`; arms carry `_sb_`. Off arms replay run 20261019T_15.

5. **Evidence.** Deaker et al. 2020 (juveniles persist on coralline algae up to 6.5 y and resume growth when coral is supplied); Wilmes et al. 2020 (the coral diet switch is needed to mature).

6. **Failure criteria:** as for the bank retest (Δv3 > 0.5, G2′ < 4/5, 1990s peak > 1 y later, no earlier 2020s, or diver error not lower).
   - **Tests:** COTSMod unit tests (off replays legacy; on gives matured `J·f`, retained `J·(1−f)`, no starvation; exclusivity) and an ADRIA env parse check.


## Scarcity-only juvenile bank: results (2026-10-19)

**Implementation.**
- COTSMod `starvation_bank` (opt-in; exclusive with `preferred_maturation` and `juvenile_storage`); COTSMod tests 395/395 pass.
- ADRIA env `COTS_JUVENILE_STARVATION_BANK`; runner experiment `bank2`.

Run 20261019T_21 (20 arms, no failures); re-score 22 (control replay 0.0); G2′ 23; audit `bank_audit.json` and `bank.png` in run 21 (`audit_1991_bank.py`). The off arms replay run 15 exactly (0.0).

| Wave | Bank | v3 | Manta part | Diver error | G2′ | 1990s peak years | First year adults ≥ 5/ha (Lizard / Eyrie) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| off | off | 11.01 | 6.48 | 3.02 | 5/5 | 1997.0 / 1997.4 / 1997.0 | 2023.2 / 2024.2 |
| off | on | **12.41** | 6.71 | **3.80** | 5/5 | unchanged | 2019.4 / 2021.8 |
| 0.25 | on | **12.37** | 6.65 | **3.81** | 5/5 | unchanged | 2019.4 / 2021.8 |

**Per reef** (v3, bank off → on, wave off):

| Reef | v3 | Diver error | Detail |
| --- | --- | --- | --- |
| Lizard | 2.47 → 2.42 | 0.53 → 0.51 | Adults 2021–25 are 8.5 / 8.1 / 8.6 / 9.4 / 11.3 against targets 3.4 / 8.3 / 5.9 / 8.6 / 9.1 (off: 3.3–9.0) |
| MacGillivray | **2.55 → 3.57** | 0.74 → 1.47 | Adults 22 / 13.5 / 11.7 in 2019–21 against targets 2.5 (2020) and 0 (2021) |
| North Direction | 2.40 → 2.74 | | Mostly the manta part |
| Eyrie | 3.60 → 3.68 | | Still too low: 4.7–6.9 in 2021–24 against 17–45 |

**Verdict: fail** (Δv3 +1.4 > 0.5; diver error +0.7 to +0.8). The 1990s are untouched, which removes the main defect of the √F̄ form, and the bank brings Lizard's 2020s on time.

The failure is MacGillivray's early release. The bank releases when modelled fast coral (food survival) recovers, and MacGillivray's model fast coral recovers 2–5× faster than its photo transects (fast share 0.59 against 0.20 in 2021; see the coral-dominance correction above). So the bank's timing is only as good as the coral recovery. The 2026 overshoot (21–32 adults/ha against 5–7) is unchanged by either bank.
