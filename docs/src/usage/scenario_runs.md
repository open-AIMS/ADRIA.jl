# Running scenarios

```julia
# Run sampled scenarios for a given RCP
rs = ADRIA.run_scenarios(dom, scens, "45")

# ... or repeat scenario runs across multiple RCPs
rs = ADRIA.run_scenarios(dom, scens, ["45", "60", "85"])

# The location of the outputs stored on disk
@info ADRIA.store_name(rs)
# "Example_domain__RCPs45__2022-10-19_12_01_26_965"

@info ADRIA.store_location(rs)
# "[some location]/Example_domain__RCPs45__2022-10-19_12_01_26_965"
```

The `rs` variable is an `ResultSet` object which acts as an interface to the stored results.

The `ResultSet` provides:

- An overview of scenarios run
- Access to results from key ADRIA metrics
- Seeding/Shading/Fogging logs
- domain spatial data

```julia
print(rs)
```

!!! note "on-disk data store"
    ADRIA uses an on-disk data store (in Zarr format) to reduce memory use.
    The primary location for these is defined in the project's `config.toml` file
    (see instructions in [Getting Started](@ref)).

!!! tip "Reloading results"
    Pre-existing results can also be reloaded by providing the path to the data store.

    ```julia
    rs = ADRIA.load_results("path to result set")
    ```
