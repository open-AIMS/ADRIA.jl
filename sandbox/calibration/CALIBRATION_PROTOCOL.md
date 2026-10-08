# COTS Calibration Protocol

This protocol defines the gates for promoting a COTS parameter set. Peak timing,
peak count, and outbreak amplitude are primary. Pointwise error and correlation
are supporting diagnostics.

## Observed Lizard Island Peak Contract

The current annual reef means and calendar-aware detector identify:

| Reef | Outbreak peak years |
|---|---|
| Lizard Island Reef | 1998, 2013 |
| MacGillivray Reef | 1996, 2015 |
| North Direction Reef | 1995, 2013 |
| Eyrie Reef | 2013 |

These values are regression-tested. Changing the detector or peak contract
requires inspecting the empirical series and recording the scientific reason;
an optimizer result alone is not justification.

## Promotion Criteria

A production candidate should satisfy all of the following on the calibration
reefs and then again on held-out data:

- detected peak count equals the observed count;
- median absolute matched-peak timing error is at most 3 years;
- each matched peak height is within 30% of the observed reef maximum;
- each matched peak prominence is within 35% of observed prominence;
- inter-peak period error is at most 3 years where two peaks are observed;
- no reef is flatlined and no unmatched high-amplitude outbreak is introduced;
- the stochastic ensemble retains the peaks in at least 70% of replicates; and
- observations fall inside the ensemble p10-p90 interval at least 60% of the
  time, reported as a diagnostic rather than optimized directly.

Failure of any primary criterion prevents language such as calibrated,
validated, or best-fit from being used without a qualifier.

## Search Sequence

1. **Focused pilot:** tune `a_F`, `a_S`, `IMM`, and initial density. Its purpose
   is pipeline validation and localization, not final calibration.
2. **Expanded no-pulse search:** tune demographic, mortality, recruitment, and
   food-limitation parameters with external supply disabled. Use multiple fixed
   seeds and compare cumulative-best curves between seeds.
3. **Mechanism test:** only if expanded internal dynamics cannot meet the peak
   contract, compare pulse and connectivity mechanisms using the same seeds,
   budgets, and objective.
4. **Production search:** run the selected mechanism with replicated candidate
   evaluations and a predeclared budget. Do not change weights after this point.

At least three optimizer seeds are required for a convergence claim. Report the
best, median, and spread of each objective component, not only total loss.

## Sensitivity and Identifiability

Run Sobol or PAWN analysis only after the expanded pilot identifies a plausible
region. Analyze peak-1 year, peak-2 year, peak heights, prominence, inter-peak
period, and coral-cover response separately. Parameters with weak total effects
should be fixed before the production search. Strongly interacting or
non-identifiable parameters should be reported rather than hidden by a single
best vector.

## Connectivity and External-Supply Experiment

The existing coral connectivity proxy is the control. A COTS-specific matrix
cannot enter calibration until it has:

- the same ordered location identifiers as the domain;
- explicit source/destination orientation and units;
- non-negative finite weights and documented normalization;
- tests for self-connectivity treatment and mass/supply scaling; and
- provenance and a distributable licence or download procedure.

Use a factorial comparison: coral proxy versus COTS matrix, crossed with pulse
off versus the selected pulse model. Hold parameters, seeds, and run budget
constant. Prefer the simpler mechanism unless held-out peak and amplitude
criteria improve materially.

The Lizard COTS matrix now meets the technical entry criteria. It is the
arithmetic mean of six 3,806-reef ReefMod spawning-season matrices, subset to
113 Lizard reefs and expanded to 2,914 sites by distributing each sink-reef
probability uniformly among its sites. The builder validates source-row /
sink-column orientation, finite non-negative values, ordered site identifiers,
and reef-to-reef mass preservation, and records source/output SHA-256 hashes.

The V2 migration is a separate domain treatment, not a replacement calibration.
`Lizard_Historical_v0.2` has 2,895 canonical `indexV2` sites, recomputed polygon
areas, and corrected nonzero historical DHW. It also carries a legacy COTS
control expanded to the V2 sites. The Owen 2018–2022 connectivity dataset is an
explicit opt-in treatment and remains below the promotion gate while non-finite
source rows, 56 unmapped external sources, and the absence of same-year
2018-2022 survival modifiers require provisional policies. A v0.2 result is not
pooled with v0.1 results; paired controls must use identical biological
parameters, observation mapping, objective, and seeds.
The area-weighted reef observation mapping is now frozen for bounded V1/V2
comparisons in `LIZARD_OBSERVATION_CONTRACT.md`. It uses the parent reef ID,
polygon areas, and a single V1-control-fitted CPUE scale. The V2 BlackBoxOptim
driver remains gated because it still uses the legacy unweighted objective and
no V2 treatment passes the qualitative peak contract. Three-year calendar
smoothing of observed CPUE is an explicitly labelled sensitivity; raw scoring
remains available and the underlying survey records are unchanged.

## Validation Design

Validation proceeds in increasing cost:

1. Leave one Lizard reef out, fit the remaining reefs, and score the held-out
   reef without refitting its biological parameters.
2. Fit observation scale on the early survey period and assess later peak
   timing and relative amplitude without re-estimating that scale.
3. Freeze parameters and validate Moore/Cairns regional domains against AIMS
   tow observations.

The calibration study repository is extracted only after inputs have stable
licensing/download instructions, all run configuration is file-backed, and a
fresh checkout can regenerate the selected tables and figures.

## Current Gate Result

No tested configuration passes the peak contract:

- Paired V1/V2/Owen sampled-year gate
  `20261001T113401_lizard_domain_gate`: three seeds, fixed biological candidate,
  fixed Allee threshold and CPUE scales, raw and three-year-smoothed survey
  scores. With smoothed observations, mean losses are 5.779 (V1 legacy closed),
  5.770 (V2 legacy closed), 5.083 (V1 legacy boundary), 4.648 (V2 legacy
  boundary), 7.122 (V2 Owen closed), and 4.489 (V2 Owen boundary). Owen plus
  boundary matches 12 of 21 observed reef-seed peaks, but its main 2017–18
  peaks remain only 0.026–0.190 COTS/tow. The independently defined outbreak
  threshold is 0.22 COTS/tow. Smoothing does not change the failed gate.
  `cycle`-mode seed repetitions were identical as expected; the `sample` run
  is a hydrodynamic-year stress test, not a historical hindcast.
- Owen source audit `20261001T115318_owen_coverage_audit`: the 56 unmatched
  sources contribute a mean 0.14% but up to 10.6% of *raw* external connection
  weight at a sink-year, mostly in the 2022 matrix. The 2010–17 mean mapped
  survival gives a median 1.84× the 2017-weighted external coefficient. No
  same-year survival or missing-source habitable area is available, so neither
  gap is evidence-resolved and Owen remains experimental.

- V2 same-seed legacy/Owen by closed/evidence-boundary pilot
  `20261001T103104_lizard_v2_connectivity_seed20260930`: mean losses 5.799,
  5.515, 7.148, and 4.511 respectively. Owen/closed detected no peaks;
  Owen/evidence-boundary detected five across four reefs, matching one observed
  peak per reef, but the main peaks were late (2017–18) and below 0.22
  COTS/tow after the frozen conversion. The upstream production ceiling and
  Owen missing-data/survival policies are provisional, so this is a mechanism
  screen, not validation or a promoted parameter estimate.

- focused `peak_amplitude_pilot_v3`: mean 0.75 matched peaks per reef, timing
  penalty 0.8125, amplitude penalty 0.894;
- expanded no-pulse `expanded_peak_pilot_seed20260930`: loss 4.819, mean 0.75
  matched peaks, and no second simulated peaks;
- expanded pulse `expanded_pulse_peak_pilot_seed20261001`: loss 5.034, mean
  0.75 matched peaks, and no second simulated peaks;
- same-seed COTS-connectivity `expanded_cotsconn_pilot_seed20260930`: best
  loss 4.942, mean 0.75 matched peaks, and no second simulated peaks; and
- post-2005 holdout: zero matched peaks across all four reefs, with mean
  amplitude penalty 1.0.

Across all six paired expanded candidates, COTS-specific connectivity lowered
mean loss by 0.262 and slightly lowered timing/amplitude penalties on average,
but neither connectivity mode produced a second peak in any candidate/reef
row. This is not a biologically material improvement.

Formal Sobol/PAWN analysis and multi-seed production optimization are deferred:
there is no plausible two-peak parameter region to analyze. A six-candidate
exploratory rank screen points first to `a_ricker`, then initial density and
immigration, but is suitable only for prioritizing mechanistic recurrence tests.

The age-stage diagnostic further shows a low-density/Allee trap. In 2024, the
four calibration reefs have high body condition (0.91-0.96) and recovered coral
cover (0.41-0.46), but only 0.094-0.114 adults. Against the best candidate's
Allee threshold of 3.20, the fecundity multiplier is approximately
0.00085-0.00127. The next mechanism pilot must therefore verify that all model
states use COTS/ha, fix `COTS_ALLEE_THRESHOLD=3.0`, and separately log local
fecundity, background immigration, internal dispersal, and external boundary
supply. Lower thresholds are counterfactual sensitivity cases, not calibration
parameters.

The six-treatment recurrence factorial has now been run with the same
demographic candidate and seed and with the Allee threshold fixed at 3.0
COTS/ha. Static and annual connectivity each produced one early peak per reef.
Source survival reduced the total to three peaks; adding the tested external
boundary did not change that result. Juvenile storage retained only one peak,
and storage plus the open-habitat proxy retained none. Mean losses were 4.941,
4.945, 5.175, 5.175, 6.425, and 7.125 respectively. No treatment produced a
second simulated peak, so all fail the mechanism gate.

The q3baseline survival factors reduce mean internal immigration from 0.216 to
0.000234 COTS/ha/year under the frozen candidate. The tested external boundary
adds only `1.55e-5` COTS/ha/year on average. `COTS_EXTERNAL_SOURCE_DENSITY` is a
legacy environment-variable name: the supplied quantity is pre-dispersal
recruit production in recruits/ha/year, not adult COTS density. It must be
derived from an upstream production model or observations before another
boundary experiment. The open-habitat term is an unvalidated proxy and must not
be described as rubble-dependent recruitment.

The recruitment audit has selected an opt-in separated model: pre-pelagic
fecundity, source-specific pelagic survival, and sink settlement are now
distinct stages. A time-varying boundary uses 1991-2025 GBR predictions of the
probability that CPUE exceeds 0.22 COTS/tow. Probability is treated as expected
source activity and combined with source-to-sink connectivity and habitable
source/sink area ratios. The prediction-to-ReefMod mapping, aliases, one missing
source, ignored prediction-only IDs, and all hashes are recorded.

The bounded separated-production pilot does not pass promotion. With
hydrodynamic forcing averaged across the six available seasons, a
scale-preserving closed treatment produces no peaks. An evidence boundary based
on production at the Allee density also produces no peaks. Production at the
Ricker maximum produces two Eyrie peaks (2001 and 2017) but none at the other
three reefs; the Eyrie peaks are too early/late and an order of magnitude too
small after observation scaling. Cycling annual matrices amplifies the result
and yields Eyrie peaks in 2001 and 2019, demonstrating sensitivity to repeated
use of the strong 2010 hydrodynamic season.

The mechanism can therefore generate recurrence without lowering the Allee
threshold, but it fails timing, amplitude, and cross-reef criteria. The next
gate is validation of area normalization and of the lag from predicted upstream
adult outbreak through spawning and settlement to destination adult abundance.
Do not tune settlement or juvenile maturation to compensate for those unresolved
source and transport contracts. No production optimization or formal
sensitivity analysis starts yet.

Leave-one-reef-out runs are supported through `BBO_EXCLUDE_REEFS`, but should
not be promoted as cross-validation until a training fit passes the peak
contract. Regional validation requires a historical Moore/Cairns domain; the
available Moore package begins in 2025 and cannot test the historical tow
series.

The 2026-10-01 V1 old-wave replay (`20261001T132602_lizard_peak_mechanism_replay`)
used frozen area-weighted reef mapping and CPUE conversion. The old candidate
with coral-proxy transport matched all seven smoothed observed peaks only in
the labelled `1 COTS/ha` Allee counterfactual. With the threshold fixed at
`3 COTS/ha`, that treatment matched two, and switching to the COTS matrix
matched three. The latter treatments did not restore the late waves. The
archived environment is incomplete, so the replay is not certified identical
to the historical run. The counterfactual must not be promoted: the threshold
remains `3 COTS/ha` in further mechanism experiments.

Coral-cover diagnostics add a second observable, not a post-hoc term in the
frozen COTS score. The initial exact-model-name lookup was **incorrect**: it
found only Eyrie because the LTMP files use the frozen COTS observation reef
aliases. The corrected `plot_lizard_ltmp_coral.py` mapping uses reef-level
`HC`/`MANTA` rows in `reef_manta.csv` for all four reefs and reef-level
`HARD CORAL`/`GROUP_LEVEL` rows in `reef_photo_transect.csv` for Lizard,
MacGillivray, and North Direction. The photo composition rows are excluded;
Eyrie has no matching photo-transect series. Both methods are 9 m observations
and remain separate on the revised plots, with source intervals and hashes.

In 2015, manta medians were 0.079, 0.084, 0.068, and 0.105 for Lizard,
MacGillivray, North Direction, and Eyrie; the available photo medians were
0.085, 0.118, and 0.088. The sampled V2 Owen/boundary seed 20260930 modeled
0.280, 0.283, 0.235, and 0.308 as area-weighted whole-reef cover. The
multi-reef gap is substantial, but the 9 m survey and whole-reef model are not
like-for-like measurement supports. Earlier years show reef-specific
differences in both directions. In closed fixed-Allee runs modeled coral rises
while adult COTS remain low, so modeled food scarcity is not the immediate
cause of their missing wave. The corrected survey coverage elevates coral
response and spatial observation mapping to an explicit validation gate;
it does not authorize fitting coral into the frozen COTS score post hoc.

The fixed-3 stage-survival screen (`20261001T133241_lizard_stage_survival_screen`)
restored the older default `m2=0.2` and/or `m3=0.1` as labelled sensitivities
on the provisional V2 Owen/boundary path. Lower adult mortality produced a
larger late Lizard peak (up to 0.254 COTS/tow with both older values), but
reduced matched peaks from four to three and worsened mean frozen loss from
4.643 to 4.864. With both older mortality defaults, 2015 modeled coral
remained 0.198, 0.206, 0.152, and 0.247 across the four reefs, above their
9 m manta medians (0.079, 0.084, 0.068, and 0.105). The screen used
Owen-mean normalization of pre-pelagic
fecundity (4532 recruits/adult/year), rather than the paired gate's V2-legacy
normalization (6389); interpret only within-screen contrasts. No stage
settings pass promotion or authorize a broader optimizer.

The next technical audit, `20261001T170000_owen_flux_audit`, replays the
paired gate's frozen V2 Owen/boundary candidate and verifies the implemented
production -> source survival -> reef transport -> site settlement -> recruit
state accounting to maximum absolute residual about `3.4e-13 COTS/ha`.
All 156 target reef-years exactly reproduce the archived gate's selected
hydrodynamic years and adult/coral trajectories to floating-point precision.
It independently checks that every local 2018-22 Owen forcing column is the
RME q3baseline **2017** source factor, and that Lizard's three Owen nodes are
aggregated only after site flux reconstruction. It finds that 2008-18
external pelagic supply (`0.153-0.309 COTS/ha/year` across target reefs)
dominates local-network supply (`0.0020-0.0085`). This is an arithmetic
contract, not validation that Owen matrices exclude mortality or that reused
2017 factors describe 2018-22. An equal-site source-production mean is the
current equation; at North Direction an area-weighted mean would be about
`1.15` times it in the median 2008-18 year. Test that proposed change only
as an opt-in, mass-checked treatment against the frozen control.

The separate `20261001T172000_ltmp_coral_mapping_audit` confirms paired
manta/photo LTMP reef-years (`23`, `24`, `23` for Lizard, MacGillivray, North
Direction) and no Eyrie photo series. Method offsets differ in sign and
magnitude; do not pool or rescale them with a universal conversion. The
frozen observation policy is now: manta is a four-reef survey diagnostic;
photo is a three-reef independent check; model cover is area-weighted
whole-reef; neither 9 m LTMP method is yet a direct observation of the V2
model, whose 2,895 sites all have a confirmed placeholder median depth 5 m.
These cover comparisons
must remain out of the frozen COTS objective.

An explicitly assumption-first screen now brackets the unresolved Owen
mortality contract without changing defaults. Treating Owen connectivity as
already including larval mortality, and globally matching mean incoming flux
to the paired control, produced two formal model maxima at each reef when
upstream boundary production was doubled. Jointly quadrupling local and
upstream production matched all seven observed peak years under each of three
recorded hydrodynamic seeds. These are exploratory scaling treatments, not
estimates of fecundity or external production. The frozen score was unchanged.

Peak count is not a sufficient qualitative gate. The three-year-smoothed
modeled inter-peak minimum divided by the smaller adjacent peak remained
0.41-0.53 in the joint-fourfold treatment, versus 0-0.005 in the observed
Lizard, MacGillivray, and North Direction series. An oracle 2000-10 shutdown
of upstream boundary flux lowered modeled ratios only to 0.17-0.21; local
recruitment and maturation persisted. Raising constant adult mortality and/or
Ricker density dependence to their configured upper bounds produced ratios
about 0.30-0.47, with timing or match penalties. These tests show that the
current equations have recurrent-maximum capacity under a simple transport
assumption, but the tested production, boundary, and constant-loss settings
do not reproduce the observed outbreak/crash contrast. They do not establish
that every existing-parameter combination fails or that a new mechanism is
required. Run IDs, fixed controls, flux normalization, and plots are in the
calibration README.

Before a broader fit, prospectively freeze an inter-peak depletion diagnostic
with a survey-noise tolerance; retain the existing COTS objective unchanged
for comparisons already made. A bounded next factorial should discriminate
sustained internal recruitment, adult persistence, and intermittent external
supply using the same Allee threshold, observation operator, and coral
guardrails. This exploratory shape test does not replace the *promotion*
evidence gate: obtain Owen methods to determine whether matrix entries include
larval mortality, resolve 2018-22 source survival and unmatched external
source areas, and build a depth/habitat-aware coral observation operator.
Do not tune maturation to compensate for provisional boundary or
support-mismatched coral residuals. No default change or production optimizer
is authorized by this screen.

### Frozen inter-wave crash screen (2026-10-02)

Before running the next diagnostic, freeze this *screening* definition: for
each of Lizard, MacGillivray, and North Direction, detect two model peaks
with the existing three-calendar-year smoothing and peak detector, then take
the minimum three-year-smoothed modeled COTS/tow strictly between them and
divide by the smaller detected peak. A convincing crash requires this ratio
to be at most 0.10 on **all three** reefs, with two detected peaks and no
loss of observed-peak matches relative to the embedded joint-fourfold
reference (seven of seven at the seed used here). The observed
ratios are 0-0.005. The 0.10 cutoff is a deliberately generous qualitative
allowance for survey/tow variation, not a statistical confidence bound or a
replacement for peak height, timing, coral, and mass-balance gates. Eyrie has
only one observed peak and is excluded from this *trough* test, but remains in
the other comparisons. Zero or tiny simulated peaks cannot pass by a ratio
alone: retain the frozen score and report peak heights in COTS/tow.

Hypothesis: continuous internal settlement and/or external supply sustain the
mid-cycle adults; if both stop, residual adults and maturation may still
prevent a near-zero trough. The diagnostic intervention sets only the
internal larval-settlement scalar to zero during calendar years 2000-10,
without changing adult mortality, maturation, production, or matrix units.
Separately, the existing boundary coefficients may be set to zero over the
same years. Compare an exact baseline replay, external-only blackout,
internal-only blackout, and both blackouts at one seed; add both blackouts
with configured-upper adult mortality only as an adult-persistence bound.
The window and upper mortality are oracle counterfactuals, not inferred
historical forcings or evidence-based parameter values. Log external,
internal, local retention, maturation, adult density, coral cover, seed,
hydrodynamic year, input hashes, and failures. The expected signature is a
selective reduction in inter-wave recruitment followed by a lower adult
trough; unchanged trough despite both blackouts implicates adult carry-over
or pre-existing juveniles. Reject the mechanism as an outbreak explanation
if it only passes under the oracle schedule, loses credible peak heights, or
degrades coral. Roll back by leaving the opt-in blackout years unset; the
legacy/default path must replay exactly.

The completed seed-20260930 factorial (`20261002T120000_owen_recruitment_blackout`)
passed adult/coral replay to `1e-12` and exact hydrodynamic-year replay in
its reference treatment, and verified source-flux zeros in the blackout
treatments. Relative to the
joint-fourfold reference's trough ratios 0.405-0.518, external-only blackout
gave 0.169-0.198, internal-only gave 0.270-0.366, and both off gave
0.051-0.071. The double blackout lost one of seven matched observed peaks;
adding configured-upper `m3=0.3` yielded 0.037-0.048 and restored seven
matches. Only this last *oracle* treatment passes the new trough screen.
Its mean frozen loss worsened from 2.923 to 3.494, first peaks were delayed
and smaller, Eyrie's observed high peak remained far underpredicted, and
modeled 2015 whole-reef coral cover increased. The screen therefore fails
promotion despite showing that the equations can crash when recruitment is
interrupted long enough. It does not establish that the actual upstream and
internal supplies stopped for eleven years, or that upper-bound mortality is
biologically justified. Do not run a production optimizer or promote the
blackout switch. The next discriminating work is to derive independently
defensible time variation in local settlement and external supply, preserve
peak amplitude/decline timing and coral guardrails, then replicate credible
treatments across seeds.

### Upstream probability and food-survival discrimination (2026-10-02)

Before changing another equation, reconstruct the external boundary from the
source outbreak probabilities, source areas, and annual Owen source-to-sink
matrices. Compare linear `p`, `p^2`, and `p^4` as **offline counterfactuals**,
using one global scalar for each nonlinear response to preserve the four-reef
2012-15 mean incoming flux. This normalization must not be fitted separately
to reefs or trough years. Keep the archived joint-fourfold treatment,
hydrodynamic-year selection, `3 COTS/ha` Allee threshold, and observation
mapping unchanged. Record source/data hashes and check the reconstruction
against the archived boundary and logged external flux. The hypothesis is
that predicted outbreak probability is too linearly mapped to upstream
production, allowing quiet-year sources to sustain the trough. The expected
signature is reduced 2004-09 boundary flux without unacceptable loss of the
1994-97 first-wave supply. Reject the proposed map if it mainly suppresses
the first peak or cannot address continuing internal recruitment. Roll back
by retaining the linear map; no default is changed by this offline audit.

In parallel, treat low-food survival as a *separate* hypothesis: at low coral
cover the current starvation factor multiplies both juvenile maturation and
adult survival. A future single-seed, fixed-forcing screen may compare the
reference with only `C_max`, only `p_tilde`, and both at their configured
upper bounds, then report trough ratio, peak heights/years, frozen loss,
coral, and stage fluxes. These bounds are not independent biological evidence.
Do not vary the Allee threshold, production, boundary strength, or scoring
rule in that screen. `eta_starve` and `fecundity_gate` must **not** be tuned as
currently exposed: they are passed into COTSMod but are not used by the
transition. Wiring either setting requires a separately specified, opt-in
core change and its own compatibility tests.

The completed fixed-seed `20261002T150000_owen_food_survival_screen` replayed
the archived control and tested those three configured-upper treatments.
None met the 0.10 trough criterion. `p_tilde=1.0` left trough ratios nearly
unchanged (0.409-0.519 on the three two-wave reefs); `C_max=1.0` gave
0.407-0.502 where two peaks remained, removed MacGillivray's first formal
peak, and reduced matches from seven to six. At Lizard the absolute trough
fell, but the first peak also shrank and shifted three years later; 2015
modeled whole-reef coral rose from 0.097 to 0.166 against the 9 m manta
median 0.079, subject to the observation-support caveat. Internal settled
immigration fell strongly but the unchanged external supply continued.
This rejects **constant, globally stronger food-survival loss alone** as the
missing crash explanation within this bounded screen; it does not exclude a
different evidence-supported, time-varying mechanism. Do not tune these
bounds further or promote them on a trough ratio alone. Next specify the
smallest independently motivated factorial that can make internal and
external recruitment jointly intermittent while preserving both waves;
retain the frozen objective, observed-peak mapping, coral guardrail, and
fixed Allee threshold.

### External-source condition proxy test (pre-registered 2026-10-02)

The available GBR prediction supplies annual source-reef probability of
exceeding `0.22 COTS/tow`, but not adult density or maternal condition for all
3,692 external sources. LTMP coral observations cover only a subset. Therefore
the exact shared density/condition production equation cannot yet be tested
for the full external boundary. As a *structural counterfactual only*, define
an upstream condition proxy for each external source `u`:

`q[u,1991]=0.8`;
`q[u,t]=(1-1/tau)*q[u,t-1]+(1/tau)*(1-p[u,t-1])`, with
`tau=8.333153725198475` from the frozen local candidate. Replace the linear boundary
source weight `p[u,t]` with `p[u,t]*q[u,t]^2`. Apply one global multiplier to
preserve mean incoming four-reef 2012-15 boundary flux, not reef-specific or
trough-specific multipliers. Keep Owen matrices, area conversion, seed,
local/internal equations and parameters, settlement, observation operator,
and `3 COTS/ha` Allee threshold fixed. The `1-p` food signal is an explicit
proxy assumption, **not** observed coral or measured larval health. The
upstream density/fertilisation part of the proposed shared equation remains
untested; do not insert probability as adult density. Assume Owen already
embeds pelagic mortality, so do not apply a second pelagic modifier.

Reconstruct the archived linear boundary and logged external flux before
running; require their established numerical tolerances. Compare control
and condition-proxy treatment under the same hydrodynamic seed. Record
per-reef first/second peak timing and heights, three-year-smoothed trough
ratio, frozen loss, coral, external/internal/local recruitment, failures,
and all input/output hashes. A useful treatment must meet the frozen
`<=0.10` trough screen on all three two-wave reefs with seven matches and
credible peak amplitudes/coral; a lower absolute trough accompanied by a
smaller first peak is not success. Reject the proxy as a calibrated source
model even if it passes: source condition needs independent GBR coral or
adult-density evidence. Roll back by using the archived linear boundary;
neither the core transition nor the default forcing changes.

The completed builder `20261002T160100_owen_source_condition_proxy` replayed
the archived linear boundary to `4.5e-16` and required one global late-wave
normaliser `1.58115`. It scarcely reduced 2004-09 incoming supply: Lizard
increased 0.229 to 0.236, MacGillivray remained 0.179, and North Direction
fell 0.379 to 0.359 pre-settlement recruits/ha/year. The paired one-seed
run `20261002T170000_owen_source_condition_screen` exactly replayed the
control and kept seven observed-peak matches, but its modeled trough ratios
were 0.505, 0.433, and 0.488 at Lizard, MacGillivray, and North Direction
(control 0.481, 0.405, 0.518). Mean frozen loss fell from 2.923 to 2.838,
illustrating why the distinct trough gate is needed. Reject this proxy as a
deep-trough mechanism; do not infer that real maternal condition is
irrelevant. The exact external density/condition production equation remains
untested until independent source histories or a full-GBR source model are
available. Do not fit a stronger proxy to the Lizard trough.

### Fecundity-by-consumption screen (pre-registered 2026-10-02)

Hypothesis: a larger first COTS cohort combined with faster per-adult coral
consumption depletes prey, reducing body condition and starvation survival
after the peak strongly enough to make a deep trough, while recovered coral
permits the second peak. This is a test of existing equations, not an assumed
result or a change in mortality/Allee meaning. Fertilisation is already
represented by `D^2/(D^2+3^2)` multiplying adult density in the local
source-production equation. Do not change the fixed `3 COTS/ha` threshold.

Use one seed (`20260930`) and a 2-by-2 factorial around the archived
joint-fourfold Owen treatment. Local `larval_fecundity` is either its archived
value or `1.5` times that value. Fast- and slow-coral consumption coefficients
`a_F` and `a_S` are either their archived values or both multiplied by
`1.5`; all other parameters, hydrodynamic selection, external boundary
production, settlement, area-weighted observation mapping, and frozen score
stay fixed. The 1.5 multiplier is a diagnostic perturbation, **not** an
evidence-derived bound. The historic optimiser allowed `a_F<=2.0` and
`a_S<=0.9`, but the current factor table lists a narrower `a_F<=1.0`; the
archived candidate already has `a_F>1.0`. Record this contract mismatch and
do not promote any candidate until bounds are reconciled with evidence.

Expected interaction: fecundity alone may amplify both peaks and troughs;
consumption alone may lower coral and COTS; their combination should deepen
the inter-peak minimum *relative to maintained peak heights* if prey-depletion
feedback is sufficient. Reject as a cycle explanation if all three two-wave
reefs do not meet the pre-frozen trough ratio `<=0.10` with seven observed
peak matches, or if timing, amplitude, coral, or frozen loss degrade
unacceptably. Log local production, retention, internal/external supply,
maturation, adult density, coral, selected hydrodynamic year, failures, and
source/output hashes. The control must replay archived adult and coral
trajectories to `1e-12`. No production optimisation or default change follows
from a one-seed pass; replication and evidence bounds remain promotion gates.

The screen was completed as
`20261002T180000_owen_fecundity_consumption_screen` (one seed, four
treatments, no failures). The control reproduced archived adult and coral
trajectories to `1e-12` and the exact selected hydrodynamic years. Although
all treatments retained seven observed-peak matches, the peak detector found
an extra early peak at MacGillivray under either single perturbation and at
all three two-wave reefs under both perturbations. The frozen mean loss
worsened from `2.923` (control) to `3.416` (fecundity), `3.194`
(consumption), and `3.556` (both). Where exactly two modeled peaks remained,
trough/smaller-peak ratios ranged from `0.405` to `0.540` on the three
two-wave reefs, still far above `0.10`. A ratio is deliberately undefined
when three modeled peaks appear; the extra peak is itself a qualitative
failure. Mean modeled 2004-09 COTS/tow across those reefs fell from `0.224`
to `0.148` in the combined treatment, but this is not a deep inter-wave
crash and comes with worse first-wave timing and/or height. The 2000-10
external pelagic supply averaged `0.477` pre-settlement recruits/ha/year in
all treatments, while the combined treatment reduced mean internal
immigration from `0.473` to `0.105`. Thus the tested prey-depletion feedback
can reduce local/internal supply, but cannot by itself synchronize a deep
trough with the unchanged external boundary. Coral in the combined treatment
is lower than the control on three of four reefs in 2015. Reject this
parameter-only explanation under the registered gate. Do not optimize this
factorial or promote its bounds; obtain independent source histories or a
source-state model before reconsidering the shared recruitment hypothesis.

### Low-food mortality trigger diagnostic (pre-registered 2026-10-02)

Hypothesis: the legacy food-survival multiplier is too weak **during the
inter-wave years relative to both outbreak waves**. Before adding a new core
equation, replay the archived joint-fourfold Owen control and reconstruct
site-level pre-feeding prey cover, body condition, food-survival multiplier,
and juvenile/adult food losses from existing state and flux logs. Check the
reconstruction against the legacy equation and require the archived adult,
coral, and hydrodynamic trajectory replay. No core or default change is
authorized by this diagnostic alone.

Use one fixed, uncalibrated linear stress screen for each of three proposed
inputs: pre-feeding prey cover (C) with threshold `0.10` cover fraction;
pre-update body condition (B) with threshold `0.50`; and prey cover per adult
`C/max(D,0.01)` with threshold `0.08 cover-fraction ha/COTS`. Stress is
`max(0, 1-x/threshold)` and a hypothetical additional survival penalty would
be `0.5*stress` on both post-baseline juvenile and adult survival. These
thresholds are diagnostic contrasts, **not evidence-derived biological
bounds**. Do not fit them to the observed trough. Compare adult-weighted
stress in 1994-97 (first wave), 2004-09 (trough), and 2012-15 (second wave).
The qualitative pre-gate requires trough stress at least twice the larger
wave stress, and at least `0.10` in absolute stress, on all three two-wave
reefs. Failure means this trigger cannot selectively deepen the trough in
the archived state trajectory; do not wire it into the core. Passing would
justify a separate opt-in dynamic implementation with unit/integration tests,
fixed `3 COTS/ha` Allee threshold, unchanged external boundary, frozen
scoring, and a rollback switch. The current model remains the control.

The completed site replay `20261002T201000_owen_food_trigger_diagnostic`
reconstructed legacy food survival to `4.4e-15`. Direct prey-cover stress
passed the temporal pre-gate on Lizard, MacGillivray, and North Direction:
trough stress was `0.200`, `0.217`, and `0.240`, respectively, at least
`2.84` times the larger wave stress. Condition and prey-per-adult stress
failed because they penalised the second or first wave too. Proceed with
**only** the direct-cover dynamic treatment below; the other two candidates
are not core-change candidates in this cycle.

Dynamic hypothesis: a second, low-cover food-survival factor, independently
controlled from `C_max`, may deepen the inter-wave crash without the broad
first-wave suppression caused by raising `C_max`. For juvenile and adult
post-baseline survival, multiply the legacy `f(C)` by
`1 - lambda*max(0,1-C/0.10)`, with pre-feeding prey cover `C=F+S` as a
dimensionless fraction. Compare the archived control (`switch=false`) to
`lambda=0.5` and `lambda=1.0` (`switch=true`), one recorded seed (`20260930`),
and no other parameter changes. These strengths and the `0.10` threshold are
bounded diagnostic perturbations, not independent biological estimates.
Keep production, external boundary, connectivity, source-survival assumption,
settlement, area-weighted observation mapping, frozen score, and `3 COTS/ha`
Allee threshold fixed. Expected: reduced 2004-09 adult survival and
maturation, a trough/smaller-peak ratio closer to `<=0.10`, both original
peaks retained in timing and amplitude, and no unacceptable coral decline.
Reject if the trough gate fails on any of the three two-wave reefs, seven
observed-peak matches are not retained, or first/second peaks or coral
degrade materially. Log local/internal/external fluxes and food-stress
diagnostics. The rollback is the default `switch=false`; no optimization or
promotion follows from a one-seed pass.

The completed dynamic screen `20261002T210000_owen_low_cover_mortality_screen`
replayed archived adult and coral trajectories to `1e-12` and hydrodynamic
years exactly; no scenario failed. At `lambda=0.5`, all seven observed peaks
remained matched, but the Lizard/MacGillivray/North Direction trough ratios
were `0.472/0.407/0.500` against control `0.481/0.405/0.518`, and mean
frozen loss worsened from `2.923` to `3.041`. At `lambda=1.0`, Lizard's first
formal peak disappeared, matches fell to six, the remaining MacGillivray
and North Direction trough ratios were `0.391/0.477`, and mean loss worsened
to `3.680`. Across the three two-wave reefs in 2000-10, mean effective
adult food survival only moved `0.850 -> 0.837 -> 0.833`, local fecundity
`1.401 -> 1.251 -> 1.152`, and internal immigration `0.473 -> 0.434 ->
0.401` pre-settlement recruits/ha/year. External pelagic supply remained
`0.477` in all cases. Total food-related deaths per hectare need not rise
monotonically as strength increases because fewer adults remain and coral
recovers; distinguish death rates from deaths. The offline timing signal
therefore did **not** survive coupled dynamics. Keep the mechanism off and
do not increase mortality further to fit the trough. The next gate remains
an independently defensible history of outside-reef production or another
joint internal/external recruitment explanation, not optimization of this
one-seed screen.

### Owen initialization and low-boundary diagnostic (2026-10-07)

Hypothesis: the inherited spatial COTS initialization, rather than external
larval supply, causes prematurely high early CPUE and forces a muted calibrated
growth response. The area-weighted, fixed-scale seed audit must precede this
test because the first model year of `cots_log` is zero-initialized and the
spatial seed includes juvenile and subadult cohorts as well as adults.

Use the archived joint-fourfold Owen control, 2018-22 annual global COTS
matrices sampled with seed `20260930`, `3 adult COTS/ha` Allee threshold,
unchanged area-weighted reef mapping and frozen three-year-smoothed score.
Change only the spatial seed multiplier (archived value versus one-quarter of
it) and the absolute outside-source pre-pelagic production (archived value,
10%, or zero), a 2 x 3 one-seed factorial. The 10% and zero settings test
whether local and within-domain production can sustain the waves; they do not
assert that outside reefs were historically absent. Log initial age-class
densities separately from modeled 1986 onward states, annual internal/external
fluxes, COTS peaks and coral cover. Require the full/full cell to reproduce
the archived joint4 control to `1e-12` for adult and coral trajectories.

Qualitative pass: all three two-wave reefs retain both observed waves, the
inter-peak smoothed minimum is at most 10% of the smaller adjacent peak, and
neither wave's height nor coral trajectory degrades materially. Loss of the
first wave or the second wave is a failure even if mean loss improves. This
screen is not a new calibration, does not alter the default, and does not
authorize optimization. The rollback is the unchanged archived control.
An adult seed near 3/ha is **not** part of this factorial: with the existing
6:3:1 juvenile:subadult:adult initializer it would simultaneously preload
roughly 18 juvenile and 9 subadult COTS/ha. A separate opt-in stage-profile
experiment with an evidence-based cohort structure would be needed to test
that hypothesis cleanly.

Completed runs: `20261007T101500_owen_initialization_audit` and
`20261007T103000_owen_seed_boundary_screen`. The initial audit matched all
23 recorded Owen source/output hashes. The full/full simulation exactly
replayed archived adult and coral trajectories and hydrodynamic years. Seeded
area-weighted adult densities were only `0.337`, `0.389`, `0.100`, and
`0.357 COTS/ha` at Lizard, MacGillivray, North Direction, and Eyrie;
no target site began at the fixed 3/ha adult Allee threshold. The corresponding
initial mean Allee multipliers were `0.0129`, `0.0166`, `0.0035`, and
`0.0143`. The archived 1985 COTS log zeros are unfilled output, not the seed.
The V2 seed field itself was sampled from
`COTS_prob_0.02_cpue_year2025_clean.tif`, so it is a 2025 spatial proxy in a
1985 initialization; 2,881 of 2,895 sites and all 113 reefs are seeded.
Historical initialization and spatial-support validity are therefore a
separate, prior gate. Do not infer an ecological growth-rate deficiency from
this particular 2025-seeded transient alone.

At full seed/full boundary, seven observed peaks were matched, but the three
two-wave trough ratios were `0.481/0.405/0.518`. Quarter seed/full boundary
retained six matches and ratios `0.387/0.366/0.454`, but shifted the first
modeled Lizard/MacGillivray peaks from 1997 to 2002. Full seed/10% boundary
and full seed/zero boundary each retained only three matches and no two-wave
target; quarter seed/10% retained five matches with ratios
`0.693/0.693/0.713`; quarter seed/zero retained none. No cell passed the
qualitative gate. In particular, the fixed parameter set cannot sustain the
observed two waves with negligible external production. This does **not**
establish structural impossibility across parameter space. Keep the existing
control and Allee threshold; do not promote any seed or boundary change.

Next discriminating gate: define a historically defensible 1985 seed location
and age-class profile, then test an opt-in adult-near-threshold seed with an
independently specified, low juvenile/subadult cohort (rather than raising
the current 6:3:1 seed multiplier), plus a bounded local-production response,
under zero/low external boundary. The experiment needs an explicit age-class
initialization contract, unit and integration tests, initial-state logs, and
the same frozen peak/coral score. A pass must recover both observed waves and
the trough without importing unsupported outside production. Otherwise,
local transmission/feedback or regional source dynamics remain limiting.

### Experimental 1991 observation-based initialization (2026-10-07)

Hypothesis: the anachronistic 2025 COTS probability seed and 2023-derived coral
cover obscure the early outbreak transient. Moving the experimental simulation
start to 1991 and initializing from nearby 1989-1991 LTMP observations may
shift the first peak without suppressing the second, while a historical coral
surface may alter predation and crash depth. This changes **initial conditions
and the time window only**, not any COTS transition equation or default domain.

`build_1991_initial_surfaces.py` pools raw COTS counts and tow effort by reef
over 1989-1991 (COTS/tow) and averages available annual 9 m LTMP manta hard-
coral medians (cover fraction). It uses audited reef-name/ID matches, including
an area-weighted centroid for the four Lizard Isles components. For each V2
site, COTS IDW weights are `tows / max(distance_km, 1)^2`; coral IDW weights
are `1 / max(distance_km, 1)^2`. Only donor reefs within 111 km count, with at
least three required. No missing estimate may be silently converted to zero.
The 1991 GBR outbreak prediction is keyed to the site parent reef and recorded
as `P(COTS/tow > 0.22)`; it is not multiplied by CPUE. The existing outside-
reef boundary still uses its prediction series independently. The builder
records source/output hashes, unmapped names, donor counts, dimensions, and
units without altering raw data.

`run_1991_initialization_test.jl` slices the actual 1985-2024 DHW cube to
1991-2024. It compares (1) a 1991 inherited-seed/inherited-coral control and
(2) a 2 x 2 factorial: adult-only IDW COTS seed under provisional observation
factors 0.015 or 0.150 **COTS/tow per adult COTS/ha**, crossed with inherited
or interpolated 1991 coral cover. The low factor is only a threshold-equivalence
sensitivity; the high factor brackets the fitted nuisance observation scales.
Neither is a validated manta-tow conversion. The opt-in site-keyed initial
state CSV supplies juvenile, subadult and adult COTS/ha separately; juvenile
and subadult densities are zero in this bounded adult-only test, deliberately
avoiding the inherited 6:3:1 preload. The interpolated 9 m manta cover is an
explicit **whole-reef proxy**, not a calibrated reef-wide observation model:
only total starting cover is changed, preserving inherited coral group/size
proportions. This habitat-support mismatch must remain visible in plots and
interpretation. AIMS observations used for initialization (1989-1991) are
excluded from scoring; scoring starts in 1992 using the frozen reef-specific
observation factors and objective. Keep the supported adult Allee threshold
at 3 COTS/ha and the Owen source/sink, annual hydrodynamic-year, and external
production settings fixed.

Expected observable consequences: later first peak relative to the inherited
1991 control, sensitivity of peak height and coral consumption to the 1991
cover proxy, and possibly a deeper trough if early coral depletion controls
survival. Log per-reef initial adult density, coral cover, hydrodynamic year,
adult and coral trajectories, local and external fluxes, matched peaks and
errors for every cell. No treatment is promoted if it loses either observed
wave, leaves the two-wave trough ratio above 0.1, or degrades coral dynamics
materially. A one-seed pass is a qualitative gate only; replication, held-out
years/reefs, source-survival sensitivity, and an independently supported
COTS/tow-to-density conversion remain required. Rollback is the unchanged
1985-2024 control and the opt-in switch unset.

```powershell
python -m unittest discover -s sandbox/calibration -p test_1991_initial_surfaces.py -v
python sandbox/calibration/build_1991_initial_surfaces.py --output sandbox/calibration/runs/<unused_surface_run_id>
$env:OWEN_1991_SURFACE_DIR = (Resolve-Path sandbox/calibration/runs/<unused_surface_run_id>).Path
julia --project=sandbox sandbox/calibration/test_1991_initial_state.jl
julia --project=sandbox sandbox/calibration/run_1991_initialization_test.jl
```

### Bounded 1991 cohort-history and inter-wave grazing diagnostic (2026-10-07)

Hypothesis: omitting immature 1991 cohorts delays the first wave; COTS grazing
after the first wave prevents coral recovery and weakens the second wave. Keep
the established 3 adult/ha Allee threshold, ecological parameters, Owen
sampled-year matrices, external boundary, coral initial cover, seed, frozen
reef observation mapping, and score unchanged. This is a diagnostic, not an
optimization or a proposal to remove COTS grazing from the default model.

First, compare IDW adult seeds at the provisional `0.015 COTS/tow per adult/ha`
factor with 0, 0.5, and 1.0 times the **model's own** 6:3:1
juvenile:subadult:adult initialization ratio; also check the full ratio at
`0.150`. The ratios stand in for unknown pre-1991 demographic history, not
measured cohorts. Control is the inherited 1991 seed. Expected result if
missing history matters: earlier first peak without losing its height or the
second wave. Failure is a later or collapsed first peak, an inter-wave trough
still above 0.1 of the smaller peak, or worse coral/held-out behavior.
Formally, for adult seed `A_i = CPUE_i / alpha`, set juvenile `J_i = 6kA_i`
and subadult `S_i = 3kA_i` for `k in {0, 0.5, 1}`. The diagnostic grazing bypass
evaluates the ordinary COTS predation/demography step on a copy of cover but
sets post-predation coral cover equal to pre-predation cover only in 2004-2012;
all other years use the ordinary transition. This preserves COTS demographic
updates in-window while deliberately removing direct coral consumption.

Second, compare the inherited 1991 trajectory to (a) a truly zero-COTS coral
control, which must have zero stage densities and no external/internal larval
supply, and (b) an opt-in 2004-2012 no-grazing diagnostic that still advances
COTS stages and recruitment. Treatment (a) bounds the *aggregate* COTS effect;
it cannot isolate remnant grazing because the entire earlier history changes.
Treatment (b) starts after the first modeled wave and isolates the direct
inter-wave predation term. A response in (b) supports the remnant-grazing
hypothesis but is not, by itself, a viable ecological mechanism. Record
population, coral, fluxes, hydrodynamic years, inputs, failures and scores.
Any temporary core-model bypass must be behind an explicit default-off switch
with an integration test proving that the default path replays. Rollback is
simply to leave the switch unset.

The COTS/tow-to-adult COTS/ha factor is an observation parameter, not an
ecological transition parameter. Before promotion, obtain either a published
calibration applicable to LTMP *COTS counts per two-minute tow* or paired
LTMP counts/effort and independently surveyed adult densities (COTS/ha), with
reef, year, survey area, habitat and detection information. Do not treat tow
search area alone or a scars-based translation as a count-to-density
calibration. Fit/propagate its uncertainty in the observation model and then
repeat the bounded screen under supported bounds. The 1991 cohort profile
likewise requires pre-1991 juvenile/subadult observations or an explicitly
model-derived warm-up state before being called historical initialization.

Outcome of the bounded screen: the full model-derived 6:3:1 cohort with IDW
adults and provisional factor `0.015` moves Lizard's first peak to 1997, but
the worst two-wave trough ratio rises to `0.594` (gate `<=0.10`). A true
zero-COTS 1991 control, including a zero external boundary, reaches `0.637`
Lizard coral in 2012 versus `0.130` with COTS. The 2004-2012 grazing bypass
leaves pre-2004 adult/coral trajectories unchanged and raises 2012 Lizard
coral to `0.482`; however, the trough ratio worsens from `0.286` to `0.451`
and the second peak falls from `0.610` to `0.546 COTS/tow`. Thus ongoing
grazing explains much of the constrained coral recovery *in this model*, but
coral release alone does not fix the COTS crash/second-wave shape. No treatment
passes the qualitative gate. The first no-COTS diagnostic was invalid because
the experimental external boundary was still on; its run is retained, not
used for inference. The corrected run asserts zero COTS in every state row.

The two user-supplied 2021-2026 workbooks were audited without fitting a
conversion. Six explicit reef crosswalks yield 14 SALAD–EOTR manta reef-years,
9 surveyed within 30 days, and 13 with culling records; only 6 of the 9
near-date groups have nonzero manta COTS counts. SALAD's `No. COTS / area`
reconciles to its all-size `Density(Ha)`, but this is not adult-only density.
Manta data in this workbook are EOTR, not the historical LTMP method used in
the calibration score. A chain through RosettaCOTS requires SALAD's
COTS-plus-scars definition, the cull-selection process, and uncertainty on
both links; it is not a replacement for a direct count-to-density study.
The size audit found 991 individually recorded COTS in 22 selected reef-years,
all with a size: 949 are at least 250 mm and 929 at least 260 mm. Only 13/22
reef-years reconcile exactly with the SALAD track counts; only 6/9 near-date
manta groups reconcile, of which 3 have positive manta COTS counts. Do not
fit a point conversion from these aggregates. Next resolve count mismatches,
perform a spatial/date track-to-tow match, agree the adult-size cutoff, and
fit a zero-aware observation model with transfer uncertainty for historical
LTMP. Hold the ecological parameters fixed meanwhile.

A coordinate/date proximity screen tightens the decision: within 1 km of
SALAD track midpoints and 30 days, 40 SALAD tracks record 106 COTS near 153
distinct EOTR manta tows that record zero COTS. Even at 0.5 km/7 days, 14
SALAD tracks record 51 COTS near 30 zero-count tows. These are not verified
shared search paths and do not yield a detectability estimate, but they rule
out treating reef-year ratios as a robust scalar conversion. Before fitting,
validate path intersections and survey protocols, resolve which size class
the ecological adult state represents, and decide how to model zero counts,
observer/detection variation, and the EOTR-to-LTMP transfer.
