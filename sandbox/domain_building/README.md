# Lizard domain build and connectivity contract

`Lizard_Historical_v0.1` remains the reproducible calibration control. The V2
migration is additive and is written to `sandbox/data/Lizard_Historical_v0.2`.
No builder in this directory overwrites v0.1 unless that version is selected
explicitly.

## V2 source-to-domain mapping

The V2 source polygons contain 2,895 rows but only 1,436 unique supplied
`site_id` values because some sites are multipart fragments. The domain therefore
uses `indexV2`, asserted to be unique and exactly `1:2895`, as its authoritative
order. Canonical IDs are `FNG_V2_0001` through `FNG_V2_2895`; the supplied ID is
retained as `source_site_id`.

The spatial builder:

- reads `farNorthGbrSubreefsV2LL.shp` as EPSG:4326;
- repairs the six invalid polygon geometries with `make_valid`, retaining only
  polygonal components;
- calculates area in EPSG:32755 and stores it in square metres;
- samples the COTS suitability raster at each projected polygon centroid;
- joins parent reef names from the supplied global reef shapefile; and
- records all source and output hashes in `provenance.json`.

The interim workbook `conMatFngbrDist 1.xlsx` has no row or column labels. Its
`in` sheet is streamed without modifying the workbook, asserted to be 2,895 by
2,895, finite, non-negative, and symmetric, then labelled in strict `indexV2`
order. Rows are sources and columns are sinks. The values are interim
distance-based connectivity weights; they are not relabelled as empirically
estimated larval probabilities.

Historical DHW is rebuilt from `dhw_historical.csv` using string-normalized RME
reef IDs. The source contains two identical copies of every reef-year value;
the builder verifies equality before dropping the duplicate copy. This fixes the
all-zero v0.1 DHW artifact and is therefore a material domain-data change that
must not be mixed with v0.1 calibration results.

Initial cover reproduces the v0.1 ReefMod extraction semantics: mean 2023 cover
for species 2-6 is divided by the ReefMod reef-specific habitable fraction and
split into 35 group-size bins. It is not yet reconciled to the Lizard spatial
layer's fixed `k = 0.5`; that correction is intentionally deferred as a separate
ecological/state-definition change.

## COTS connectivity datasets

The loader recognizes `ADRIA_COTS_CONNECTIVITY_DATASET`:

- `legacy` (default) reads `cots_connectivity/` and retains the six ReefMod
  2010-2017 spawning-season matrices as the control;
- `owen_global_2018_2023` reads
  `cots_connectivity/datasets/owen_global_2018_2023/` and is experimental and
  opt-in.

The Owen source order is the row order of `reefShapefiles/gbrShapeLL.shp` and
contains 3,861 reefs. All 113 V2 parent reefs match. The builder keeps annual
matrices separate, uses source-row/sink-column orientation, expands the static
mean to sites by dividing each reef-to-reef value among sink sites, and checks
that local reef-to-reef mass is preserved.

The supplied Owen and ReefMod ancillary data do not yet form a complete
evidence-backed forcing contract:

- the first four Owen matrices contain wholly non-finite source rows;
- 56 Owen source reefs lack ReefMod larval-survival and habitable-area data; and
- the available survival file ends in 2017, before the 2018-2022 Owen
  hydrodynamic years.

Consequently, the Owen builder fails by default. The currently built
experimental artifact was generated only with three explicit provisional
policies: zero wholly non-finite rows, exclude the 56 unmapped external sources,
and reuse the latest available 2017 survival modifier. These choices are logged
in `provenance.json` and `annual_forcing/provenance.toml`; they are not model
defaults and have not passed the peak-timing/amplitude promotion gates.

The checked-in `audit_owen_forcing_gaps.py` quantifies, but does not fill, these
gaps. Its `20261001T115318_owen_coverage_audit` output finds that only ten
unmatched source-year rows have nonzero raw connection to Lizard after the
explicit zero-row policy. The missing share of raw external connection weight
averages 0.14% across sink-years but reaches 10.6% at one sink in 2022. Most
unmatched rows are wholly non-finite in 2018–21. For mapped sources, using
each source's 2010–17 mean survival instead of 2017 gives a median 1.84×
external-survival-weighted coefficient across sink-years. Neither comparison
identifies the actual missing reef areas or 2018–22 survival; those remain
data-provider/evidence gates.

## Reproduction commands

Run from the repository root:

```powershell
# 1. Build the versioned spatial domain, interim connectivity, cover, and DHW.
python sandbox\domain_building\build_lizard_domain_v2.py

# 2. Build the retained legacy COTS control for the V2 sites.
$env:LIZARD_DOMAIN_VERSION = 'Lizard_Historical_v0.2'
$env:COTS_WATER_QUALITY_SCENARIO = 'q3baseline'
julia --project=sandbox sandbox\domain_building\build_lizard_cots_connectivity.jl

# 3. Build the explicitly provisional Owen option.
python sandbox\domain_building\build_lizard_cots_connectivity_owen.py `
  --nonfinite-source-policy zero_full_rows `
  --missing-source-policy exclude `
  --survival-year-policy latest_available

# 4. Validate the legacy control and the full V2/Owen load contracts.
$env:LIZARD_DOMAIN_VERSION = 'Lizard_Historical_v0.2'
$env:ADRIA_COTS_CONNECTIVITY_DATASET = 'legacy'
julia --project=sandbox sandbox\domain_building\validate_lizard_cots_connectivity.jl
julia --project=sandbox sandbox\domain_building\validate_lizard_domain_v2.jl

# 5. Run one deterministic integration smoke test (seed 20260930).
julia --project=sandbox sandbox\domain_building\smoke_lizard_domain_v2.jl

# 6. Audit unmatched Owen sources and historical survival sensitivity.
python sandbox\domain_building\audit_owen_forcing_gaps.py
```

To select Owen matrices in a calibration or simulation, all of the following
must be recorded:

```powershell
$env:LIZARD_DOMAIN_VERSION = 'Lizard_Historical_v0.2'
$env:ADRIA_COTS_CONNECTIVITY_MODE = 'cots'
$env:ADRIA_COTS_CONNECTIVITY_DATASET = 'owen_global_2018_2023'
$env:COTS_CONNECTIVITY_TEMPORAL_MODE = 'cycle'  # or mean/sample by design
$env:COTS_CONNECTIVITY_SEED = '20260930'
```

The bounded same-seed V2 pilot is reproduced with:

```powershell
julia --project=sandbox sandbox\calibration\compare_lizard_v2_connectivity.jl
```

Its `20261001T103104_lizard_v2_connectivity_seed20260930` output compares
legacy/Owen with closed/evidence-boundary conditions. With Owen forcing, the
closed subset has no detected peaks at the four scored reefs; the provisional
boundary creates five detected peaks and lowers mean loss from 7.148 to 4.511.
The main peaks are late and below 0.22 COTS/tow, so neither Owen nor the
boundary is promoted. The BlackBoxOptim driver blocks V2 until an
area-weighted observation mapping and peak gate are frozen.

The subsequent three-seed paired V1/V2/Owen comparison uses the tested
area-weighted observation contract in
`sandbox/calibration/LIZARD_OBSERVATION_CONTRACT.md`. The `cycle` control
(`20261001T112847_lizard_domain_gate`) is seed-invariant by design. The
`sample` sensitivity (`20261001T113401_lizard_domain_gate`) varies selected
hydrodynamic years but is not a historical hindcast. With three-year-smoothed
observations, V1/V2 legacy closed losses are 5.779/5.770, while V2 Owen
closed/boundary losses are 7.122/4.489. The late Owen/boundary peaks remain
below 0.22 COTS/tow; see the calibration README and updated peak plot. The
V1/V2 comparison changes site geometry, total area, DHW, and raster sampling
together, so it cannot attribute effects to connectivity alone.

## Next promotion gates

1. Obtain or justify 2018-2022 source-specific pelagic-survival modifiers.
2. Resolve the 56 new source reefs with evidence-backed survival and habitable
   areas, then quantify their effect on every Lizard sink.
3. Confirm the non-finite rows with the data provider and replace the provisional
   zero-row policy.
4. The area-weighted observation aggregation is tested and frozen for bounded
   comparison, but the production optimizer still uses its legacy unweighted
   objective and remains gated. Reconcile the large V1/V2 area difference
   before interpreting ecological parameter shifts.
5. Paired V1/V2 and three-seed sampled-year screens are complete; repeat with
   evidence-backed Owen survival/area policies and a genuine temporal holdout
   before claiming peak timing, height, prominence, and period validation.
6. For a promotion-grade run, archive package versions and the working-tree
   state alongside revisions, seeds, input/output hashes, per-reef scores, and
   failed evaluations.

The Allee threshold remains fixed at its evidence-based value during these
tests. Juvenile maturation and other biological parameters are not adjusted to
compensate for unresolved domain or transport contracts.
