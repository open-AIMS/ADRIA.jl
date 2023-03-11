# An overview of ADRIA

ADRIA consists of three overarching components:

1. a set of Multi-Criteria Decision Analysis methods used to emulate decision-makers
2. a simplified coral ecology model
3. a suite of analysis and assessment methods

Each component may be applied separately, or altogether to perform an end-to-end analysis.
This page documents the underlying architectural implementation of ADRIA, detailing how the above
components interact.


## General Structure

The primary purpose of ADRIA is to support reef restoration and adaptation through the development
of robust intervention strategies under deep uncertainty. Here, "robustness" refers to the ability
of an intervention strategy to meet desired outcomes under uncertain future conditions, which
themselves are unknown and may be unexpected. To do so, ADRIA employs an Exploratory Scenario 
Modelling framework to explore the range of possible futures and their outcomes.

Exploratory Scenario Modelling (ESM) itself leverages uncertainty and sensitivity analysis (UA/SA).
Uncertainty analysis quantifies the variability of uncertainties in a given system and its
expected range of outcomes. Sensitivity analysis examines the effect of a change in a model's
inputs to its outputs. Common workflows to such analyses involve a three-step process (as
discussed in Pianosi et al., 2016):

1. Input sampling
2. Model evaluation
3. Post-processing

When ADRIA is applied in its entirety, "input sampling" is analogous to scenario generation: all
model parameters (the inputs) are collated and are sampled through a quasi-monte carlo process 
and subsequently adjusted to produce plausible combinations of parameter values. The Sobol' sampling
method is adopted as the default, although any method provided by the [QuasiMonteCarlo.jl](https://github.com/SciML/QuasiMonteCarlo.jl)
package can be used. Sample adjustment is required to map sampled values (which are continuous) to
categorical or whole number values (e.g., Baroni and Tarantola, 2014) as may be expected by some 
parameters. Values are also adjusted to avoid implausible factor combinations, such as active
intervention parameters in the case of non-intervention scenarios.

Model evaluation is simply running the model with the generated scenario set.

Post-processing is the analysis and visualization step.


## Model Factors

Factors in ADRIA are defined across four sub-components:

1. Intervention
2. Criteria
3. EnvironmentalLayers
4. Coral

Each sub-component is represented by a struct with fields for each parameter. The `Intervention`
sub-component holds parameters that define a given adopted strategy: how many (and type of) corals
are to be seeded, the length of any deployment, the start/end years, and so on. The `Criteria`
sub-component relates to the preferences for the Multi-Criteria Decision Analysis methods, further
detailed in the [MCDA page](TODO). For the ADRIA ecology model, `EnviromentalLayers` relate to the
environmental scenarios available for a given simulation (a time series of DHW and Wave stress),
itself determined on the loading of data (see [this page](TODO)).

The `Coral` sub-component relates to the coral ecosystem model, currently defined as an ecosystem of
six coral "types": Tabular Acropora, Corymbose Acropora, Small Massives and Large Massives, as well 
as enhanced (e.g., via genetic modification) Tabular and Corymbose Acropora. ADRIA represents
these across six size classes, with six parameters for each coral "type" and size class.

These six parameters are further detailed in [the Coral documentation](TODO), however, it results 
in 216 unique parameters (6 groups x 6 size classes x 6 parameters). Instead of specifying all
216 parameters by hand, ADRIA instead auto-generates the sub-component using a common template
(see [`coral_spec()`](TODO) and [`create_coral_struct()`](TODO)). Through discussion with expert
stakeholders, parameter bounds were set to +/- 10% of their default values following a triangular
distribution, the peak of which is the default value.

The [ModelParameters.jl](https://github.com/rafaqz/ModelParameters.jl) package is used to provide
a simple table-like interface to model parameters. Using the `Param` type provided by 
ModelParameters.jl allows the default value, the parameter bounds, their expected distribution,
and other associated metadata (e.g., a human-readable name and description) to be attached. An
example from the `Intervention` sub-component is shown below.

```julia
guided::N = Param(0, ptype="integer", bounds=(-1, 3 + 1), dists="unif",
        name="Guided", description="Choice of MCDA approach.")
```

Note the `+ 1` in the example above when specifying the upper bound. Usual sampling approaches,
such as the Sobol' method mentioned above, provide a continous value. These need to be 
transformed to a whole number when working with categorical factors or parameters which expect a
whole number. A "flooring trick" (as it is referred to here) is adopted to handle this 
transformation. The process is described in
[Baroni and Tarantola (2014)](https://doi.org/10.1016/j.envsoft.2013.09.022) and is illustrated
here with the example below.

Let $x_i$ be a parameter that is expected to be a discrete value between 1 and 3 (inclusive).
In other words, there are 3 valid options: $x_i = \{1, 2, 3\}$.

1. The upper bound is increased by 1, such that $x_i = \{1, 2, 3, 4\}$
2. Sample from this extended range with the given sampler, which returns a value $1 \le v_i \lt 4$, where $v_i \in \mathbb{R}$
3. Take the `floor` of $v_i$. For example, if $v_i = 3.9$, then $\text{floor}(v_i) = 3$.

In this manner the expected probability of a possible value being selected is maintained.

These parameter definitions are collectively known as the model specification, and
can be collated as a DataFrame.

!!! note "Parameter Collation and Scenario Generation"
    See [this page](TODO) for an example how-to on collating model parameters and
    generating samples.

Combination of the realized parameter values then represent a "scenario".


# References

Baroni, G., and S. Tarantola. 2014. 
A General Probabilistic Framework for uncertainty and global sensitivity analysis of deterministic models: A hydrological case study. 
Environmental Modelling & Software 51:26–34.

Pianosi, F., K. Beven, J. Freer, J. W. Hall, J. Rougier, D. B. Stephenson, and T. Wagener. 2016. 
Sensitivity analysis of environmental models: A systematic review with practical workflow. 
Environmental Modelling & Software 79:214–232.
