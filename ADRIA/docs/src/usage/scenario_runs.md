```@meta
EditURL = "scenario_runs.jl"
```

# Running scenarios

```julia
# Run sampled scenarios for a given RCP
rs = ADRIA.run_scenarios(dom, scens, "45")

# Repeat scenario runs across multiple RCPs
rs = ADRIA.run_scenarios(dom, scens, ["45", "60", "85"])

# The location of the outputs stored on disk
@info ADRIA.store_name(rs)
# "Example_domain__RCPs45__2022-10-19_12_01_26_965"

@info ADRIA.result_location(rs)
# "[some location]/Example_domain__RCPs45__2022-10-19_12_01_26_965"
```

The `rs` variable is an `ADRIAResultSet` that acts as an interface to the on-disk store.
It provides access to scenario inputs, outcome arrays, intervention logs, and spatial data.
A summary can be printed with:

```julia
print(rs)
```

Commonly accessed fields include `rs.inputs` (the scenario DataFrame), `rs.outcomes`
(named outcome arrays), `rs.ranks` (location ranking logs), `rs.seed_log`,
`rs.shading_log`, and `rs.loc_data` (spatial attributes).
See [Loading Results](@ref) for a complete field reference.

!!! note "on-disk data store"
    ADRIA uses an on-disk data store (in Zarr format) to reduce memory use.
    The primary location for these is defined in the project's `config.toml` file
    (see instructions in [Getting Started](@ref)).

!!! tip "Reloading results"
    Pre-existing results can be reloaded by providing the path to the data store.

    ```julia
    rs = ADRIA.load_results("path to result set")
    ```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

