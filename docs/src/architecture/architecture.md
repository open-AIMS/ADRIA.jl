# Architectural overview

ADRIA consists of three overarching components:

1. a set of Multi-Criteria Decision Analysis methods used to emulate decision-makers
2. a coral ecosystem model (referred to as ADRIAmod)
3. a suite of analysis and assessment methods

Each component may be applied separately, or altogether to perform an end-to-end analysis.
This page documents the underlying architectural implementation of ADRIA, detailing how the above
components interact.

## General Structure

The primary purpose of ADRIA is to support reef restoration and adaptation through the development
of robust intervention strategies under deep uncertainty. Here, "robustness" refers to the ability
of an intervention strategy to meet desired outcomes under uncertain future conditions, which
themselves are unknown and may be unexpected. To do so, ADRIA adopts an Exploratory Scenario
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
model factors (the inputs) are collated and are sampled through a quasi-monte carlo process.
The Sobol' sampling method is adopted as the default, although any method provided by the
[QuasiMonteCarlo.jl](https://github.com/SciML/QuasiMonteCarlo.jl) package can be used.
Sample adjustment is required to map sampled values (which are continuous) to categorical or
whole number values (e.g., Baroni and Tarantola, 2014) as may be expected by some
factors. Values are also adjusted to avoid implausible factor combinations, such as active
intervention factors in the case of non-intervention scenarios.

Model evaluation is simply running the model with the generated scenario set.

Post-processing is the analysis and visualization step.

## Model Factors

Factors in ADRIA are defined across four sub-components:

1. Intervention
2. CriteriaWeights
3. EnvironmentalLayer
4. Coral

Each sub-component is represented by a struct with fields for each parameter. The `Intervention`
sub-component holds parameters that define a given adopted intervention strategy/option: how
many (and type of) corals are to be seeded, the length of any deployment, the start/end years,
and so on. The `CriteriaWeights` sub-component relates to the preferences for the Multi-Criteria
Decision Analysis methods, further detailed in [Dynamic Multi-Criteria Decision Analysis](@ref). For the ADRIA ecosystem model
(ADRIAmod), `EnviromentalLayer` relate to the environmental scenarios available for a given
simulation (a time series of DHW and Wave stress), itself determined on the loading of data
(see [Running scenarios](@ref)).

The `Coral` sub-component relates to ADRIAmod, currently representing six coral species:

1. Arborescent Acropora
2. Tabular Acropora
3. Corymbose Acropora
4. Corymbose non-Acropora
5. Small massives and encrusting
6. Large massives

ADRIAmod represents these across six size classes, with six parameter sets for each coral
species and size class. These six parameter sets are further detailed in [Model Factors](@ref),
however, it results in a large number of unique factors (6 groups x 6 size classes x 6 parameters: 216 factors).
Instead of specifying all coral factors by hand, ADRIA instead auto-generates the sub-component
using a common template (see `coral_spec()` and `create_coral_struct()` in [General API](@ref)).
Through discussion with expert stakeholders, factor bounds were set to +/- 10% of their
default values following a triangular distribution, the peak of which is the default value.

The [ModelParameters.jl](https://github.com/rafaqz/ModelParameters.jl) package is used to provide
a simple table-like interface to model factors. ADRIA provides a wrapper to the `Param` type
provided by ModelParameters.jl called a `Factor`. These `Factor`s requires ADRIA-specific
metadata to be provided alongside the default value for a given model factor. These include
the assumed distribution of the factor, distribution factors (principally the lower and
upper bounds), and a human-readable name and description. An example from the `Intervention`
sub-component is shown below.

```julia
guided::N = Factor(0.0, ptype="ordered categorical", dist=(-1.0, 3.0), dist_params=DiscreteUniform,
        name="Guided", description="Choice of MCDA approach.")
```

Note that factor values are provided as floats - even where discrete values are expected -
to maintain type stability. Mixing floats with integers will lead to an error. Similarly,
values used to update the model should always be floats.

Combinations of the realized factor values then represent a "scenario".

!!! note "Parameter Collation and Scenario Generation"
    See [Cookbook examples](@ref) for an example how-to on collating model factors and
    generating samples.

# References

Pianosi, F., K. Beven, J. Freer, J. W. Hall, J. Rougier, D. B. Stephenson, and T. Wagener. 2016.
Sensitivity analysis of environmental models: A systematic review with practical workflow.
Environmental Modelling & Software 79:214–232.
