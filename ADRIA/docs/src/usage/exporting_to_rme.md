# Exporting to RME Format

The `ADRIA.RMEExport.export_to_rme()` utility allows ADRIA `ResultSet` objects to be exported into a format compatible with [ReefModEngine.jl](https://github.com/open-AIMS/ReefModEngine.jl). This is primarily used to link ADRIA ecological simulations with the **CREAM (Comprehensive Reef Economic Analysis Model)** suite via the `cost-eco-model-linker`.

## Domain Compatibility

The export targets **Great Barrier Reef (GBR) domains**. Both a standard ADRIA GBR
`Domain` and an `RMEDomain` (`load_domain(RMEDomain, ...)`) are supported. It requires:

- A `GBRMPA_ID` column in the domain spatial data (`rs.loc_data`), used to name the yearly
  reefsets in `scenario_info.json`. Current ADRIA GBR DataPackages and current RME data
  packages both provide it. Older packages without `GBRMPA_ID`, and non-GBR or synthetic
  domains, are not supported.
- The standard 5-group functional coral model and the standard seeding factors
  (`N_seed_*`, `seeding_devices_per_m2`). The evenness normalisation hard-codes 5 groups.
- An `RCP` in the results that also exists in the domain's DHW data (needed for the
  GCM-name lookup).

## Overview

Because ADRIA and RME use different internal standards, the export process performs several critical "translation" steps:

1.  **Unit Scaling**: Proportions (0.0 - 1.0) are scaled to percentages (0 - 100) where required by the economic model.
2.  **Coordinate Variables**: NetCDF outputs include human-readable coordinate variables for `timesteps` (years) and `locations` (Unique Reef IDs).
3.  **Intervention Realization**: ADRIA's spatial intervention logs (`seed_log`) are translated into absolute coral counts and realized intervention areas (km²) for logistical cost calculations.

## Usage Workflow

### 1. Generate Paired Scenarios
The economic model calculates benefits by comparing an intervention to its corresponding counterfactual. Use `sample_pairwise()` to ensure your scenario set contains the necessary pairs in the correct order.

```julia
using ADRIA

# Standard ADRIA GBR domain
dom = ADRIA.load_domain("path/to/domain", "45")

# ...or an RME data package
dom = ADRIA.load_domain(RMEDomain, "path/to/rme_data_package", "45")

# Generate 16 paired scenarios (16 Counterfactuals + 16 matching Interventions)
scens = ADRIA.sample_pairwise(dom, 16)
```

### 2. Run Simulations
Run the scenarios normally. Ensure you are using a domain compatible with the Great Barrier Reef (GBR) if you intend to use the standard CREAM cost models.

```julia
rs = ADRIA.run_scenarios(dom, scens, "45")
```

### 3. Export Results
Call `export_to_rme()` to generate the RME-compatible directory structure.

`export_to_rme()` takes the `Domain` used to run the scenarios, the `ResultSet`, and
the output directory.

```julia
using ADRIA.RMEExport

out_dir = "./economic_analysis_input"
export_to_rme(dom, rs, out_dir)
```

Optional keyword arguments: `map_counterfactuals`, `pdp_geopackage_path`, `scen_spec`.

## Exported Files

The utility produces the following files in the target directory:

*   **`results.nc`**: NetCDF containing `total_cover`, `relative_juveniles`, `relative_shelter_volume` (also written as `relative`), `evenness`, and zero-filled `rubble` / `cots` placeholders, plus `timesteps` / `locations` / `scenarios` coordinate variables.
*   **`iv_yearly_scenarios.csv`**: A year-by-year log of absolute coral counts, seeding densities, and deployment areas per scenario.
*   **`scenario_info.json`**: Metadata mapping scenario indices to counterfactual status and defining the unique "reefsets" used each year.
*   **`reef_information.csv`**: Spatial metadata for the reefs involved in the simulation.

## Technical Notes

### Shelter Volume Scaling
The downstream economic model applies a **9.33x scaling factor** to the exported `relative_shelter_volume`. This reconciles ADRIA's reference standard (relative to 1 m²) with the historical GBR standard (relative to a 95cm diameter circle).

### Counterfactual Detection
A scenario is flagged as a `counterfactual` in the export if all intervention factors (`N_seed`, `fogging`, `SRM`, etc.) are set to zero. This bitmask is stored in `scenario_info.json`.
