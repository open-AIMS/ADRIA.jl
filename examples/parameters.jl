using DataFrames
using ADRIA


@info "Loading data package"
data_pkg = "C:/development/ADRIA_brick/Brick"
scenario_domain = ADRIA.load_domain(data_pkg, 45)

# Get current parameter table
param_df = ADRIA.param_table(scenario_domain)

# Update a domain spec with new values from row of DataFrame
update_params!(scenario_domain, param_df[1, :])


# Get parameter details

## Parameter names
p_names = scenario_domain.model[:fieldname]

## Current values
p_vals = scenario_domain.model[:val]

## ADRIA parameter types
p_types = scenario_domain.model[:ptype]

## Parameter bounds (for e.g., to pass into a sampler or optimizer)
## Note: ADRIA integer parameter bounds are set to (ℓ-1, u) as the `floor(x+1)` is assigned with `update_params!()`.
p_bounds = scenario_domain.model[:bounds]

## Component groups
p_groups = scenario_domain.model[:component]

## All of above as a DataFrame
model_spec = DataFrame(scenario_domain.model)


# Get DataFrame of parameter information for a specific sub-component (Intervention, Criteria, Coral)
ADRIA.component_params(scenario_domain.model, Intervention)

