# Generating scenarios

Typical use of ADRIA is to generate a number of scenarios by sampling from combinations of
possible factors relating to environmental, intervention, and coral conditions.

A scenario is defined as a combination of all factors (i.e., all the model inputs).

```julia
# Generate 128 scenarios based on available environmental data layers and model parameters
scens = ADRIA.sample(dom, 128)
```

Here, the `scens` variable holds a DataFrame of scenarios of shape $N$ by $D$, where
$N$ is the number of scenarios (rows) and $D$ is the number of factors (columns).

The Sobol' method (Sobol' 1993, 2001) is the default sampling approach. It is a
deterministic low-discrepancy quasi-monte carlo sampler. Samples are described as
having low discrepancy if the samples are equi-distributed, and thus guarantee an
even exploration of parameter space. One limitation of the Sobol' method is that 
all factors are assumed to be independent. For most factors represented in ADRIA,
this assumption holds true. Specific factors relating to intervention options may
conditionally co-vary however, and this dependency is introduced by adjusting the
sample values _a posteriori_ to restrict sampled values to their plausible 
combinations, and to map continuous values to their expected discrete factor
values (where necessary), as is in the case with categorical factors. The Sobol'
scheme is therefore disrupted due to the adjustment and so a Sobol' sensitivity
analysis cannot be relied on. Subsequent assessment of uncertainty and sensitivity
is instead conducted with the distribution-based PAWN method (Pianosi and 
Wagener 2015, 2018).

!!! note "Sobol' samples"
    The convergence properties of the Sobol' sequence is only valid if the number of
    samples is a power of 2.

Samples for factors with non-uniform distributions are transformed to their indicated
distributions using the Inverse Cumulative Distribution Function method.

Although the Sobol' method is the default, any sampler supported by the
[Surrogates.jl](https://github.com/SciML/Surrogates.jl) package may be used.
Below is an example using Latin Hypercube sampling.

```julia
import Surrogates.QuasiMonteCarlo: LatinHypercubeSample

scens = ADRIA.sample(dom, 100, LatinHypercubeSample())
```


### On model parameters and specifications

The current default values can be extracted with:

```julia
params = ADRIA.param_table(scenario_domain)
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
model_spec = ADRIA.model_spec(scenario_domain)

# Sometimes it is useful to export the model specification to CSV
ADRIA.model_spec(scenario_domain, "model_spec.csv")
```


# References

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
