# Lizard reef observation contract for the V1/V2 gate

This contract is frozen for `compare_lizard_domain_gate.jl`. It does not change
the legacy optimizer objective or promote a new ecological parameterization.

1. The COTS adult state is a density in COTS/ha at each model site. Each target
   reef's annual state is the polygon-area-weighted mean of its sites:
   `sum(area_m2_i * adults_ha_i) / sum(area_m2_i)`. The result remains COTS/ha.
   Site membership is keyed by parent `UNIQUE_ID`, not by map row order. Site
   areas must be finite and positive, and within-reef weights sum to one.
2. The four AIMS tow-series names are mapped explicitly to Lizard Island,
   MacGillivray, North Direction, and Eyrie reefs. `reef_cots.csv` is unchanged.
   Only surveys through 2024 enter scoring because the model ends in 2024.
   The input has one survey record per target reef-year in that window, so no
   within-year averaging or tow weighting is needed.
3. The conversion from simulated COTS/ha to COTS/tow is one non-negative
   least-squares scale per reef. It is estimated *only* on the V1 legacy/closed
   treatment at seed `20260930` using raw survey CPUE, then held constant for
   every domain, transport, boundary, seed, and smoothing treatment. This is
   an exploratory full-window scale, not a temporal holdout fit.
4. The primary observation sensitivity is a symmetric three-calendar-year
   running mean: at a surveyed year, average available survey CPUE in that
   year and the adjacent year on either side. Missing years are not filled or
   interpolated; an isolated survey keeps its raw value. The raw series is
   scored and reported separately. Both scores use the same conversion scale
   and the same objective. Model peak detection retains its existing
   three-year smoothing. A five-year observed mean is **not** adopted: it
   shifted peaks further and created an Eyrie peak in the diagnostic.
5. The peak detector uses normalized height and prominence; a detected local
   peak need not exceed the independently defined outbreak threshold of
   `0.22 COTS/tow`. Report both detected peaks and their actual COTS/tow
   heights. The plot shows raw observations, the three-year observed mean,
   model trajectories, detected peaks, and the threshold separately.

The comparison outputs include the complete per-site weights, reef area
summary, observed raw/smoothed CPUE, per-reef scores, and trajectories. The
source hashes, script hashes, seeds, selected forcing mode, and known Owen
policies are recorded in each run's `metadata.toml`. A `cycle` run is
deterministic under changed seeds; `sample` mode is a hydrodynamic-year
sensitivity, not a historical hindcast. Do not pool V1 and V2 losses as if
they came from the same spatial/data domain.

## Coral-cover diagnostic mapping (not part of the COTS score)

The corrected `plot_lizard_ltmp_coral.py` uses the same four survey aliases as
item 2 above: `Lizard Island Reef` -> `Lizard Isles`, `MacGillivray Reef` ->
`Macgillivray Reef`, `North Direction Reef` -> `North Direction Island`, and
`Eyrie Reef` -> `Eyrie Reef`. The first coral plot searched for exact model
names and consequently omitted three available manta series; that plot is
superseded, not overwritten.

- Manta: `reef_manta.csv`, `data_type=manta`, `domain_category=reef`,
  `variable=HC`, `purpose=MANTA`, `project_code=LTMP`, and `depth=9`.
- Photo transect: `reef_photo_transect.csv`,
  `data_type=photo-transect`, `domain_category=reef`,
  `variable=HARD CORAL`, `purpose=GROUP_LEVEL`, `project_code=LTMP`, and
  `depth=9`. Species/family `COMPOSITION` rows are excluded so hard coral is
  not double-counted. No matching Eyrie photo series is present.
- Use `report_year` for alignment with the 1985-2024 model calendar years;
  retain fractional `date` and source `id` in the derived observation CSV.
  Assert at most one record per reef, method, and report year. Do not fill
  missing years. Show the source-reported `lower` and `upper` intervals without
  assuming they are a specified confidence level.
- Model coral cover is a polygon-area-weighted whole-reef fraction; both LTMP
  measurements represent 9 m field methods. Keep manta and photo separate,
  and treat model/survey differences as diagnostics until a spatial/depth
  observation model is justified. Do not add coral to the frozen COTS objective
  as a post-hoc weight.
