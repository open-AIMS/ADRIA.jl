# Generating scenarios

Typical use of ADRIA is to generate a number of scenarios by sampling from combinations of
possible factors relating to environmental, intervention, and coral conditions.

A scenario is defined as a combination of all factors (i.e., all the model inputs).

```julia
# Load domain before generating scenarios
dom = ADRIA.load_domain("path to domain data package")

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
PAWN method (Pianosi and Wagener 2015, 2018).

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

## Constrained sampling

At times, it is necessary to create samples while holding some model factors constant.

Although a scenario set could be modified to make specific factors constant, doing so
runs the risk of creating (many) identical scenarios, thereby wasting computational
effort. A more efficient approach is to modify the model specification itself to treat
those factors as constants. These then get ignored for the purpose of scenario
generation.

```julia
dom = ADRIA.load_domain("path to domain data package")

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

Samples can also be taken over a constrained range. For example, if one wanted to investigate only
scenarios with high fogging and seeding, the following could be used:

```julia
dom = ADRIA.load_domain("path to domain data package")

# Adjust seeding bounds. Note only lower and upper bounds are needed because the factors in
# question have a uniform distribution.
dom = ADRIA.set_factor_bounds(dom, :N_seed_TA, (500000.0, 1000000.0))
dom = ADRIA.set_factor_bounds(dom, :N_seed_CA, (500000.0, 1000000.0))
dom = ADRIA.set_factor_bounds(dom, :N_seed_SA, (500000.0, 1000000.0))

# Adjust fogging bounds. Note lower, upper and mode parameters are needed because it
# is a triangular distribution.
dom = ADRIA.set_factor_bounds(dom, :fogging, (0.2, 0.3, 0.1))

# Adjust multiple factors simultaneously.
dom = ADRIA.set_factor_bounds(dom;
    seed_heat_stress=(0.3, 0.7),
    N_seed_TA=(500000.0, 1000000.0),
    N_seed_CA=(500000.0, 1000000.0))
```

## Sampling counterfactuals only

A convenience function to create scenarios with no interventions (counterfactuals).

```julia
cf_scens = ADRIA.sample_cf(dom, 1024)
```

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
