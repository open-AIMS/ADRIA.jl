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
0.00085-0.00127. The next mechanism pilot must therefore resolve the units and
low-density range of `allee_threshold` and separately log fecundity, background
immigration, and dispersed larval supply.

Leave-one-reef-out runs are supported through `BBO_EXCLUDE_REEFS`, but should
not be promoted as cross-validation until a training fit passes the peak
contract. Regional validation requires a historical Moore/Cairns domain; the
available Moore package begins in 2025 and cannot test the historical tow
series.
