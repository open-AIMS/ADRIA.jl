# ADRIA.jl

ADRIA: Adaptive Dynamic Reef Intervention Algorithms.

[![Release](https://img.shields.io/github/v/release/open-AIMS/ADRIA.jl)](https://github.com/open-AIMS/ADRIA.jl/releases)  [![Documentation](https://img.shields.io/badge/docs-stable-blue)](https://open-aims.github.io/ADRIA.jl/stable/)  [![DOI](https://zenodo.org/badge/483052659.svg)](https://zenodo.org/badge/latestdoi/483052659)

ADRIA is a decision-support tool designed to help reef managers, modellers and decision-makers
address the challenges of adapting to climate change in coral reefs. It provides line of sight
to conservation solutions in complex settings where multiple objectives need to be considered,
and helps investors identify which options represent the highest likelihood of providing
returns on investment. ADRIA uses a set of dynamic Multi-Criteria Decision Analyses (dMCDA)
which simulates a reef decision maker to identify candidate locations for intervention
deployment which consider ecological, economic and social benefits.

ADRIA also includes a simplified coral ecosystem model to allow exploration of outcomes as a 
result of intervention decisions made across a wide range of possible future conditions.

Usage is demonstrated in the [documentation](https://open-aims.github.io/ADRIA.jl/stable/usage/getting_started/)


## Quick start

### Setup

To specify user-specific options, a `config.toml` file should be created with the following options (adjusted to suit your needs):

```toml
[operation]
num_cores = 2     # No. of cores to use. Values <= 0 will use all available cores.
threshold = 1e-6  # Result values below this will be set to 0.0
debug = false     # Disable multi-processing to allow error messages to be shown

[results]
output_dir = "./Outputs"  # Change this to point to where you want to store results
```

!!! tip "Performance"
    ADRIA uses an on-disk data store to hold results from model runs.
    Setting `output_dir` to a directory on an SSD (Solid State Drive)
    will maximize performance.


### Usage

```julia
# Import ADRIA package
using ADRIA


# Load data for a spatial domain
dom = ADRIA.load_domain("path to domain data package")

# Generate 128 scenarios based on available environmental data layers and model parameters
scens = ADRIA.sample(dom, 128)

# Run sampled scenarios for a given RCP
rs = ADRIA.run_scenarios(scens, dom, "45")

# ... or repeatedly run scenarios across several RCPs
rs = ADRIA.run_scenarios(scens, dom, ["45", "60", "85"])

# then extract metrics for analysis
tac = ADRIA.metrics.total_absolute_cover(rs)
```
