# Loading a Domain

ADRIA is designed to work with `Domain` data packages.
In short, these are pre-packaged data sets that hold all the necessary data to run
simulations for a given spatial domain.

See [Architectural overview](@ref) for more information.

A `Domain` may be loaded with the `load_domain` function.
By convention we assign the `Domain` to `dom`, although this variable can be named anything.

```julia
dom = ADRIA.load_domain("path to domain data package")
```

ReefMod Engine datasets can also be used to run ADRIAmod simulations for the Great Barrier
Reef.

```julia
dom = ADRIA.load_domain(RMEDomain, "path to ReefModEngine dataset", "45")
```

Note that at the moment the target RCP has to be specified.
