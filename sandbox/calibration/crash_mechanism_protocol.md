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

The [completed 1991 cohort factorial](runs/20261008T_crash_factorial_1991_v1/metadata.toml) passed an exact archived state/flux control replay and invariant external-supply audit. None of the three treatments met the qualitative gate. Stage storage delayed Lizard's first peak ten years and worsened the worst trough ratio (`0.594→0.691`); `m3=0.30` delayed it eight years and cut its height by 56% while the worst trough remained `0.515`. The combination lost Lizard's first detected wave altogether. See [per-reef diagnostics](runs/20261008T_crash_factorial_1991_v1/crash_factorial_audit.csv) and [paired plot](runs/20261008T_crash_factorial_1991_v1/cots_coral_crash_factorial.png). **No arm is promotable.** The observed mechanism is compensation: storage retains immature COTS and raises later local/internal recruitment; constant mortality acts before as well as after the desired first peak. The model coral discrepancy persists, with the stated observation-support caveat.

For this **diagnostic capacity test**, hold `rho=0.5` (one-year memory weight), `B_crit=1.5 adult COTS/ha`, `w_B=0.5 adult COTS/ha`, `C_crit=0.15` coral fraction and `w_C=0.03` coral fraction. Cross density-only and density×food with `h_max` of `0.75` and `1.5` per annual step, plus a zero-hazard control: five arms, one recorded hydrodynamic seed, no optimizer. The higher hazard corresponds to at most `1-exp(-1.5)=0.777` extra loss of *background-and-food survivors* in one year, not 77.7% of all adults. These are **stress-test bounds, not field-estimated pathogen parameters**. They were selected before implementation: the archived 1991 site-level adult seed has median `0.363`, 95th percentile `0.993`, and 99th percentile `3.039 COTS/ha`; modeled control reef adult 90th percentiles over 1992–2024 are `2.410` at Lizard and `2.891` at MacGillivray. Thus `B_crit` is above most seed sites but within simulated wave densities, and is deliberately unrelated to the `3 COTS/ha` fertilization half-saturation. The new run must retain site-level burden/hazard/death summaries because archived site-level *dynamic* abundance was not saved; these settings cannot be called biological bounds.

Before implementing the lagged hazard, write pure core tests for (a) zero hazard reproducing legacy transitions bit-for-bit; (b) burden recurrence and units, including `B_0=0`; (c) increasing hazard only when the configured density/food trigger increases; (d) nonnegative, mass-accounted adult survivors and added-death diagnostics; and (e) deterministic replay. Add ADRIA adapter tests for switch, validation, dimensions and per-site logging. Hold the frozen score and boundary; record first-wave peak preservation, duration of the gap, second-wave reseeding, coral and fluxes. If the hazard fails because external settlement imposes a floor, test a separate labelled boundary sensitivity rather than silently suppressing outside-domain production. No generic outbreak timer is eligible as a biological mechanism.

## Lagged-hazard decision and next structural question

The [five-arm hazard run](runs/20261008T_lagged_hazard_1991_v1/metadata.toml) and [independent audit](runs/20261008T_lagged_hazard_1991_v1/lagged_hazard_audit.json) replay the archived control COTS/coral and all 12 recruitment-flow channels exactly. All arms retain identical external immigration. The core unit tests passed `128/128`; the ADRIA focused tests passed `23/23`. The site-level audit checks `B_0=0`, the annual recurrence to `7.1e-15 COTS/ha`, nonnegative states, and that extra deaths do not exceed pre-transition adults. See [reef outcomes](runs/20261008T_lagged_hazard_1991_v1/lagged_hazard_audit.csv), [site-year activation](runs/20261008T_lagged_hazard_1991_v1/site_hazard_summary.csv), and [COTS/coral plot](runs/20261008T_lagged_hazard_1991_v1/cots_coral_lagged_hazard.png).

No arm passes. Control worst two-wave trough/smaller-peak ratio is `0.593935` with 6/7 observed peaks matched. Density-only `h_max=1.5` gives `0.600000` and 6/7, delaying Lizard's first detected peak from 1997 to 2001. Density×food at `1.5` gives `0.593493` and 6/7, also moving Lizard's first peak to 2001; the tiny trough change is not an improvement. At `0.75`, both mode variants lose Lizard's first detected wave and match only 5/7. Coral is spared slightly, but remains low relative to 9 m observations, with the acknowledged support difference.

This is a *phase failure*, not a silent implementation failure. In the strong density-only treatment, hazard-active sites comprise roughly `33%` in 1995, `51%` in 2005, but `0.5%` in 2008 and `0.3%` in 2012. The mechanism acts near the first wave and switches off in the gap; adult density and continuing settlement then maintain the floor. The food multiplier does not rescue timing because the density trigger has already collapsed. An offline recurrence check on the **fixed control trajectory** (not a new ecological simulation) shows that changing only `rho` from `0.5` to `0.75` or `0.9` at `B_crit=1.5` leaves just `0.7%` or `1.0%` of sites above threshold by 2012; `rho=0.9` also delays substantial activation past the first wave. Do not launch a broad optimizer on this formulation. A future bounded test needs a separately justified persistence/recovery state or boundary mechanism that operates *after* the first peak while allowing later reseeding.

## Cohort audit and ReefMod duration precedent (2026-10-09)

The [readout-only 1991 replay](runs/20261009T_cohort_audit_1991_v1/metadata.toml) reproduces archived adult/coral trajectories and all 12 flux channels exactly. Its [annual area-weighted cohorts](runs/20261009T_cohort_audit_1991_v1/cohorts.csv) and [Lizard plot](runs/20261009T_cohort_audit_1991_v1/cohort_floor_diagnostics.png) show that the 1999–2007 sawtooth coincides with the realized food-survival factor alternating between about `0.27` and `1.0`. Food survival already multiplies N2-to-adult maturation and adult carryover, but does not multiply N1 recruits. In 2008–12, Lizard adult food survival averages `0.986`; adult abundance averages `0.805 COTS/ha`, of which `0.667` is carryover and `0.138` is newly matured. External settlement averages `0.542 COTS/ha/year` in those years. Thus maturation contributes to the floor but does not dominate it, and ongoing settlement maintains the future cohort pipeline. The fixed Lizard plotting scale is `0.2942 COTS/tow per adult COTS/ha`; `0.1 COTS/tow` corresponds to `0.340 adult COTS/ha` under that mapping. The provisional `q=0.015` used for 1991 initialization is distinct.

The published [ReefMod 7.0 settings](https://github.com/ymbozec/REEFMOD.7.0_GBR/blob/main/settings/settings_COTS.m) specify a randomly checked `4–6` **year** outbreak duration and comment a fixed two-year Lizard calibration option. The [transition code](https://github.com/ymbozec/REEFMOD.7.0_GBR/blob/main/functions/f_runmodel.m) confirms that adult density must remain above the threshold over the checked duration; it redraws the candidate duration at a qualifying check, then resets older age classes to background while leaving younger recruits. A [regional report](https://www.barrierreef.org/uploads/CCIP-R-04-Final-Report-Regional-Modelling.pdf) describes a different-version `2–5`-year disease window, an eight-year maximum age, and a separate `5%` preferential-prey starvation trigger. These are version-specific modelling conventions, not field estimates of natural pathogen timing. ReefMod's densities and reset values are in COTS per `400 m²` grid, so numerical transfer to this model's `COTS/ha` state and provisional CPUE mapping needs explicit justification.

The next bounded comparator should be opt-in and state-triggered, with a frozen density threshold in adult `COTS/ha`: control versus adult-only reset versus adult-plus-N2 reset at two preregistered durations. Log clock age, threshold crossings, pre/post state, removed adults and N2, surviving N1, later maturation, all settlement sources, coral, and peak/trough metrics. The adult-plus-N2 arm tests whether immature refill defeats adult-only collapse; keep external N1 supply unchanged so reseeding remains a constraint. Before implementation, specify whether the trigger requires **continuous** above-threshold density as in ReefMod or simply starts a countdown on outbreak onset; those are different hypotheses. Do not fit the clock to observed peak years or call a successful timer proof of epizootics. Failure includes losing the 1997 Lizard peak, a one-year-only trough followed by immediate N2 refill, an absent second peak, or unacceptable coral/flux degradation. Any one-seed improvement remains experimental pending replicated and held-out gates.

Rollback is `COTS_LAGGED_ADULT_HAZARD=false` (the default) or `COTS_HAZARD_MAX=0`; core COTS states, observation units, the Allee threshold, and legacy defaults are unchanged. The mechanism remains experimental and **must not** be described as confirmed pathogen mortality.

## Preregistered structural comparison: preferred-prey starvation × adult senescence (2026-10-09)

Status: **preregistered before implementation; no default or objective change.** This supersedes the ReefMod timer as the *next* screen. The timer remains eligible as a later phenomenological comparator. The reason is that ReefMod adds its timer on top of an eight-year maximum age and preferred-prey starvation. All three published GBR comparators share those two foundations, and this model lacks both: ReefMod 7.0, [CoCoNet v3.4](https://research.csiro.au/coconet/wp-content/uploads/sites/486/2025/11/CoCoNet-user-guide-and-technical-summary-v3.4.pdf) Table A.1, and MICE ([Rogers & Plagányi 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9085818/)). Porting only the timer would borrow the least mechanistic part.

### Structural diagnosis motivating the test

Values below are from the archived anchor (`a_F=1.18`, `a_S=0.165`, `C_max=0.498`, `m3=0.195`) and the [cohort audit](runs/20261009T_cohort_audit_1991_v1/cohorts.csv).

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

[Run](runs/20261009T_structural_crash_1991_v1/metadata.toml) and [audit](runs/20261009T_structural_crash_1991_v1/structural_crash_audit.json): the control replay is exact, and external immigration is identical in every arm. **No arm passes.**

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

[v2 run](runs/20261009T_unit_migration_q015_1991_v2/metadata.toml) and [audit](runs/20261009T_unit_migration_q015_1991_v2/unit_migration_audit.json): the compatibility replay is exact (max relative delta `3e-15`); external immigration is identical in every arm; simulated CPUE equals `q·A`. **No arm passes.**

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
