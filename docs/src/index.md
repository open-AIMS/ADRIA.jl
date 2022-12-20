# ADRIA.jl Documentation

ADRIA is a multi-criteria decision support platform for informing reef restoration and adaptation interventions.


## Quick start


```julia
# Import ADRIA package
using ADRIA


# Load input dataset ("Input Set") for a spatial domain
dom = ADRIA.load_domain("path to Input Set")

# Generate 100 scenarios based on available environmental data layers and model parameters
scens = ADRIA.sample(dom, 100)

# Run sampled scenarios for a given RCP
rs = ADRIA.run_scenarios(scens, dom, "45")

# ... or repeatedly run scenarios across several RCPs
rs = ADRIA.run_scenarios(scens, dom, ["45", "60", "85"])

# then extract metrics for analysis
tac = ADRIA.metrics.total_absolute_cover(rs)
```

- [Synopsis](@ref)