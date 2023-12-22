@doc """
Example script showcasing how to extract model details, such as
- the model parameter table
- model specification

and more specific information/data from the above.

This example file assumes Julia is started in the examples folder.
"""

using DataFrames
using ADRIA


@info "Loading data package"
data_pkg = "./Test_domain"
scenario_domain = ADRIA.load_domain(data_pkg, 45)

# Get current parameter table (fieldnames and their values)
param_df = ADRIA.param_table(scenario_domain)

# Get model specification with lower/upper bounds separated
model_spec = ADRIA.model_spec(scenario_domain)

# Export model specification to CSV
ADRIA.model_spec(scenario_domain, "model_spec.csv")


# Get parameter details

## Parameter names
p_names = scenario_domain.model[:fieldname]

## Current values
p_vals = scenario_domain.model[:val]

## ADRIA parameter types
p_types = scenario_domain.model[:ptype]

## Parameter bounds (for e.g., to pass into a sampler or optimizer)
## Note: ADRIA integer parameter bounds are set such that ℓ ≤ x ≤ u+1,
## where ℓ is the lower bound and u is the upper bound.
## This is because `floor(x)` is assigned with `update_params!()`.
## Instances where ℓ := x := u indicate uncertain parameters that
## are nevertheless assumed to be constant.
p_bounds = scenario_domain.model[:dist_params]

## Component groups
p_groups = scenario_domain.model[:component]

## All of above as a DataFrame
model_spec = DataFrame(scenario_domain.model)


# Get DataFrame of parameter information for a specific sub-component (Intervention, Criteria, Coral)
ADRIA.component_params(scenario_domain.model, Intervention)
