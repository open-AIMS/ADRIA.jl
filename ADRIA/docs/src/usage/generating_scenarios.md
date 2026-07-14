# Generating scenarios

Typical use of ADRIA is to generate a number of scenarios by sampling from combinations of
possible factors relating to environmental, intervention, and coral conditions.

A scenario is defined as a combination of all factors (i.e., all the model inputs).

```julia
# Load domain before generating scenarios
dom = ADRIA.load_domain("path to domain data package", "<RCP>")

# Generate 128 scenarios based on available environmental data layers and model parameters
scens = ADRIA.sample(dom, 128)
```

Here, the `scens` variable holds a DataFrame of scenarios of shape $N$ by $D$, where
$N$ is the number of scenarios (rows) and $D$ is the number of factors (columns).
Because it is a DataFrame, it can be modified after the fact.

The Sobol' method (Sobol' 1993, 2001) with Owen Scrambling is the default sampling approach.
The Sobol' sampling method is a deterministic low-discrepancy quasi-monte carlo method.
Samples are described as having low discrepancy if the samples are equi-distributed, and
thus guarantee an even exploration of factor space. Owen Scrambling is an approach to
introduce randomness to the quasi-monte carlo sequence, and belong to a class of sampling
methods known as randomized Quasi-Monte Carlo (rQMC). rQMC approaches offer a balance
between good space-filling properties and associated exploration of factor space with
improved convergence characteristics.

One limitation of the Sobol' method is that all factors are assumed to be independent.
For most factors represented in ADRIA, this assumption holds true. Specific factors relating
to intervention options may conditionally co-vary however, and this dependency is introduced
by adjusting the sample values _a posteriori_ to restrict sampled values to their plausible
combinations, and to map continuous values to their expected discrete factor
values (where necessary), as is in the case with categorical factors. The Sobol'
scheme is therefore disrupted due to the adjustment and so a Sobol' sensitivity
analysis may exhibit comparatively poor convergence. Subsequent assessment of
uncertainty and sensitivity is instead conducted with the distribution-based
PAWN method (Pianosi and Wagener 2015, 2018). See
[Conditional factor dependencies](@ref) below for what this adjustment does in practice,
and the [architecture overview](@ref "Conditional sampling dependencies") for how it is
implemented.

!!! note "Sobol' samples"
    The convergence properties of the Sobol' sequence is only valid if the number of
    samples is a power of 2.

Samples for factors with non-uniform distributions are transformed to their indicated
distributions using the Inverse Cumulative Distribution Function method.

Although the Sobol' method is the default, any sampler supported by the
[QuasiMonteCarlo.jl](https://github.com/SciML/QuasiMonteCarlo.jl) package may be used.
Below is an example using Latin Hypercube sampling.

```julia
import ADRIA.QuasiMonteCarlo: LatinHypercubeSample

scens = ADRIA.sample(dom, 100, LatinHypercubeSample())
```

## On model parameters and specifications

The current default values can be extracted with:

```julia
params = ADRIA.param_table(dom)
```

Again, `params` is a DataFrame of a single row and $D$ factors:
A single scenario with model factors set to their default values.

Running specific user-defined scenarios is as simple as modifying the DataFrame
(referred to as the "scenario specification"). A set of scenarios may be specified
simply by extending the number of rows of the DataFrame. Details of the ADRIA
model - parameter names, the default values, and their bounds - can be extracted
as well.

```julia
# Get model specification
model_spec = ADRIA.model_spec(dom)

# Sometimes it is useful to export the model specification to CSV
ADRIA.model_spec(dom, "model_spec.csv")
```

## Fixing factor sampling bounds

At times, it is necessary to create samples while holding some model factors constant.

Although a scenario set could be modified to make specific factors constant, doing so
runs the risk of creating (many) identical scenarios, thereby wasting computational
effort. A more efficient approach is to modify the model specification itself to treat
those factors as constants. These then get ignored for the purpose of scenario
generation.

```julia
dom = ADRIA.load_domain("path to domain data package", "<RCP>")

# Could keep a copy of the original model parameters/bounds
# to reset to later.
# orig_spec = DataFrame(dom.model)

# Make the assisted adaptation factor a constant
ADRIA.fix_factor!(dom, :a_adapt)

# Set the assisted adaptation factor to a given constant value
ADRIA.fix_factor!(dom, :a_adapt, 3.0)

# Pass in factor names and their constant values as named arguments
# to fix a set of factors.
ADRIA.fix_factor!(dom;
    N_seed_TA=Int64(5e5),
    N_seed_CA=Int64(5e5),
    SRM=0.0,  # Never shade
    fogging=0.0,  # Never fog
    a_adapt=3.0,  # only deploy +3 DHW enhanced corals
    seed_years=5,
    shade_years=0,
    seed_deployment_freq=0,
    seed_year_start=3,
    shade_year_start=3,
    seed_coral_cover=1.0
)
```

## Setting different sampling bounds

Samples can also be taken over a constrained range. For example, one can investigate only
scenarios with high fogging and seeding, and select a specific MCDA decision method:

```julia
dom = ADRIA.load_domain("path to domain data package", "<RCP>")

# Adjust seeding bounds. Note only lower and upper bounds are needed because the factors in
# question have a uniform distribution.
ADRIA.set_factor_bounds!(dom, :N_seed_TA, (500000.0, 1000000.0))
ADRIA.set_factor_bounds!(dom, :N_seed_CA, (500000.0, 1000000.0))
ADRIA.set_factor_bounds!(dom, :N_seed_SA, (500000.0, 1000000.0))

# Adjust fogging bounds. Note lower, upper and mode parameters are needed because it
# is a triangular distribution.
ADRIA.set_factor_bounds!(dom, :fogging, (0.2, 0.3, 0.1))

# Adjust multiple factors simultaneously (more efficient than setting these one at a time)
ADRIA.set_factor_bounds!(dom;
    seed_heat_stress=(0.3, 0.7),
    N_seed_TA=(500000.0, 1000000.0),
    N_seed_CA=(500000.0, 1000000.0))

# List of available MCDA decision methods
ADRIA.decision.mcda_method_names()

# Constrain guided/counterfactual regime, and separately select a specific MCDA decision method
ADRIA.set_factor_bounds!(dom, :guided, (0.0, 1.0))  # unguided or guided only, no counterfactual
ADRIA.set_factor_bounds!(dom, :mcda_method, ("COCOSO",))
```

## Conditional factor dependencies

Some factors are only meaningful under certain scenario regimes. For example,
seeding-related deployment factors have no effect on a counterfactual (no intervention)
scenario, and decision-strategy criteria weights have no effect unless a `guided` MCDA
approach is selected. Rather than sampling these factors unconditionally and discarding
the result, ADRIA resolves such dependencies automatically, either by:

- fixing the affected factors to a constant *before* sampling, when the governing setting
  (currently `guided`) is known for the whole call (e.g. `ADRIA.sample_cf`,
  `ADRIA.sample_guided`, `ADRIA.sample_unguided`), or
- zeroing/adjusting the affected factors *after* sampling, on a per-scenario (row) basis,
  when the governing setting is itself sampled (as with plain `ADRIA.sample(dom, n)`, where
  `guided` varies from row to row) or when the dependency is only evaluable from a drawn
  value (e.g. gating on whether a sampled `fogging` value is greater than zero).

See the [architecture overview](@ref "Conditional sampling dependencies") for how this is
implemented, if extending or debugging this behavior.

This means the returned scenario DataFrame will commonly contain columns set to `0.0` (or
another sentinel constant) for rows where the corresponding intervention/strategy is
inactive, rather than a value drawn from that factor's nominal distribution.

```julia
cf_scens = ADRIA.sample_cf(dom, 128)

# Intervention and criteria weight columns are fixed to a sentinel value for
# counterfactual scenarios rather than sampled, since they have no effect.
cf_scens[:, [:N_seed_TA, :fogging, :seed_heat_stress]]
```

Marine Cloud Brightening factors are handled the same way — see
[Marine Cloud Brightening (MCB) Scenarios](@ref) below for a concrete example of fixing
those factors to values available in a given dataset.

!!! warning "Effect on Sobol' sensitivity analysis"
    Fixing or adjusting sampled values after the fact disrupts the low-discrepancy
    structure of the underlying Sobol' sequence. Scenario sets produced this way should
    **not** be treated as valid Sobol' samples for the purpose of estimating Sobol'
    indices — such samples exist only to explore plausible factor combinations. Where
    Sobol' indices are required, PAWN is used instead (see
    [PAWN sensitivity (heatmap overview)](@ref)), which does not require this structure to
    be preserved.

!!! note "Duplicate scenarios"
    Because multiple distinct pre-adjustment samples can collapse onto the same fixed
    values (e.g. every counterfactual sample ends up with `fogging = 0.0`), it is expected
    for a scenario set to contain some duplicate rows. ADRIA will emit a warning
    reporting the proportion of duplicates when this occurs.

## Marine Cloud Brightening (MCB) Scenarios

When a domain is loaded with a 5D DHW dataset (containing `mcb_durations` and `albedo` dimensions), ADRIA automatically populates MCB-specific intervention factors. These are prefixed with `mcb_` and their sampling distributions are derived from the NetCDF axis labels.

The primary MCB factors are:
- `mcb_albedo`: The reflectiveness level to apply.
- `mcb_duration`: The yearly duration (in days) of MCB deployment.
- `mcb_deployment_freq`: How often to deploy (e.g., every 1 year, every 2 years).

Note that `mcb_start_year` is currently hardcoded to **2035**.

Because these factors are tied to specific levels available in the provided NetCDF, it is often necessary to fix them to a specific value or adjust their bounds to match the dataset's constraints.

```julia
dom = ADRIA.load_domain("path/to/5d/domain", "45")

# Fix MCB to a specific duration and albedo level available in the NetCDF
ADRIA.fix_factor!(dom, :mcb_duration, 50.0)
ADRIA.fix_factor!(dom, :mcb_albedo, 0.3)

# Or allow them to vary across their available categorical range
# (Default behavior if not fixed)
scens = ADRIA.sample(dom, 128)
```

## Sampling by intervention regime

Alongside plain `ADRIA.sample(dom, n)` (which samples the `guided` factor together with
everything else, so each row may land in a different intervention regime), ADRIA provides
a family of convenience functions that fix the intervention regime for the whole call
up-front. This avoids wasting samples on factor combinations that are meaningless for the
regime under study (e.g. seeding-related factors when only counterfactual scenarios are
wanted), and lets the conditional factor dependencies described above be resolved *before*
sampling rather than row-by-row afterwards (see
[Conditional factor dependencies](@ref)).

| Function | Regime sampled | `n` refers to |
|----------|-----------------|---------------|
| `ADRIA.sample(dom, n)` | All regimes (`guided` itself is sampled) | total scenarios |
| `ADRIA.sample_cf(dom, n)` | Counterfactual only (no interventions) | total scenarios |
| `ADRIA.sample_unguided(dom, n)` | Unguided interventions only | total scenarios |
| `ADRIA.sample_guided(dom, n)` | Guided (MCDA-driven) interventions only | total scenarios |
| `ADRIA.sample_balanced(dom, n)` | Equal parts counterfactual, unguided and guided | scenarios *per regime* (returns `3n` rows) |
| `ADRIA.sample_selection(dom, n)` | Guided location-selection factors only; coral factors held at their defaults | total scenarios |
| `ADRIA.sample_paired(dom, n)` | Guided scenarios paired 1:1 with a derived counterfactual row sharing the same non-intervention factor values | guided scenarios drawn (returns `2n` rows) |
| `ADRIA.sample_paired_stratified(dom, n, ssps)` | As `sample_paired`, repeated independently for each RCP/SSP in `ssps` | guided scenarios drawn *per SSP stratum* |
| `ADRIA.sample_matched(dom, n)` | Guided scenarios paired 1:1 with derived rows for each regime in `regimes` (default: counterfactual **and** unguided) | guided scenarios drawn (returns `(1 + length(regimes)) * n` rows) |
| `ADRIA.sample_matched_stratified(dom, n, ssps)` | As `sample_matched`, repeated independently for each RCP/SSP in `ssps` | guided scenarios drawn *per SSP stratum* |

```julia
# Counterfactual scenarios only (no interventions)
cf_scens = ADRIA.sample_cf(dom, 1024)

# Unguided intervention scenarios only
ug_scens = ADRIA.sample_unguided(dom, 1024)

# Guided (MCDA-driven) intervention scenarios only
gd_scens = ADRIA.sample_guided(dom, 1024)

# 1024 scenarios of each regime (3072 rows total)
balanced_scens = ADRIA.sample_balanced(dom, 1024)
```

For causal comparisons between an intervention and "no intervention", `sample_paired`
draws guided scenarios and derives a matching counterfactual row for each one, holding all
non-intervention factors (environmental layers, etc.) fixed between the pair so that
outcome differences can be attributed to the intervention itself. `sample_paired_stratified`
extends this across multiple RCP/SSP scenarios, using an independent Owen scramble per
stratum:

```julia
# 128 guided/counterfactual pairs (256 rows), tagged by `:run_type`
paired_scens = ADRIA.sample_paired(dom, 128)

# As above, repeated for each listed RCP/SSP, tagged by `:ssp` and `:run_type`
strat_scens = ADRIA.sample_paired_stratified(dom, 128, ["45", "60", "85"])
```

`sample_paired` only pairs guided scenarios against a derived counterfactual row, so it
cannot support a 3-way cross-comparison against unguided scenarios. `sample_matched`
generalises this: it draws guided scenarios and derives a matched row for each regime
requested via `regimes` (any non-empty subset of `(:counterfactual, :unguided)`, default
both), all sharing the same non-intervention factor draws — `sample_paired` is just
`sample_matched` with `regimes=(:counterfactual,)`. Note that unguided rows keep the same
intervention deployment amounts as their guided counterpart (only the MCDA criteria
weights/`plan_horizon` are zeroed and the strategy is fixed to periodic) — unlike
counterfactual rows, where deployment amounts are zeroed too:

```julia
# 128 guided/counterfactual/unguided triples (384 rows), tagged by `:run_type`
matched_scens = ADRIA.sample_matched(dom, 128)

# Only guided + unguided (256 rows)
matched_ug_scens = ADRIA.sample_matched(dom, 128; regimes=(:unguided,))

# As `sample_matched`, repeated for each listed RCP/SSP
strat_matched_scens = ADRIA.sample_matched_stratified(dom, 128, ["45", "60", "85"])
```

This differs from `sample_balanced`, which draws each regime *independently* — those rows
are not matched row-for-row and cannot be used for this kind of differencing.

!!! note "Sobol' samples"
    As with plain `ADRIA.sample`, `n` (or `n_per_stratum`) should be a power of 2 when
    using the default Sobol' sampler.

## References

1. Sobol’, I. M. 1993.
   Sensitivity analysis for non-linear mathematical models.
   Mathematical Modelling and Computational Experiment 1:407–414.
   [Translated from Russian, accessible at: https://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=mm&paperid=2320&option_lang=eng]
2. Sobol′, I. M. 2001.
   Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates.
   Mathematics and Computers in Simulation 55:271–280.
   https://doi.org/10.1016/S0378-4754(00)00270-6
3. Pianosi, F., and T. Wagener. 2015.
   A simple and efficient method for global sensitivity analysis based on cumulative distribution functions.
   Environmental Modelling & Software 67:1–11.
   https://dx.doi.org/10.1016/j.envsoft.2015.01.004
4. Pianosi, F., and T. Wagener. 2018.
   Distribution-based sensitivity analysis from a generic input-output sample.
   Environmental Modelling & Software 108:197–207.
   https://dx.doi.org/10.1016/j.envsoft.2018.07.019
