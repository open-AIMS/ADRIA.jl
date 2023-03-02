# Generating scenarios

Typical use of ADRIA is to generate a number of scenarios by sampling from combinations of
possible factors relating to environmental, intervention, and coral conditions.

A scenario is defined as a combination of all factors (i.e., all the model inputs).

```julia
# Generate 128 scenarios based on available environmental data layers and model parameters
scens = ADRIA.sample(dom, 128)
```

Here, the `scens` variable holds a DataFrame of scenarios of shape
``N`` by ``D``, where ``N`` is the number of scenarios (rows) and ``D``
is the number of factors (columns).

!!! tip "Alternate samplers"
    The sampling method is compatible with any sampler supported by the 
    [Surrogates.jl](https://github.com/SciML/Surrogates.jl) package.

    The default sampler is Sobol'. Below is an example using 
    Latin Hypercube sampling.

    ```julia
    import Surrogates.QuasiMonteCarlo: LatinHypercubeSample

    scens = ADRIA.sample(ex_domain, 100, LatinHypercubeSample())
    ```


### On model parameters and specifications

The current default values can be extracted with:

```julia
param_df = ADRIA.param_table(scenario_domain)
```

Again, `param_df` is a DataFrame of ``1`` by ``D``:
A single scenario with input values set to their default.

Running specific user-defined scenarios is as simple as modifying
the scenario DataFrame (referred to as the "scenario specification", 
or "scenario spec").

On a related note, details of the ADRIA model - parameter names, the 
default values, and their bounds - can be extracted as well.

```julia
# Get model specification
model_spec = ADRIA.model_spec(scenario_domain)

# Sometimes it is useful to export the model specification to CSV
ADRIA.model_spec(scenario_domain, "model_spec.csv")
```