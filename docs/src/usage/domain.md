# Loading a Domain

ADRIA is designed to work with `Domain` data packages.
In short, these are pre-packaged data sets that hold the all necessary data to run
simulations for a given spatial domain.

See [Architecture](@ref) for more information.

A `Domain` may be loaded with the `load_domain` function.
By convention we assign the `Domain` to `dom`, although this variable can be named anything.

```julia
dom = ADRIA.load_domain("path to Input Set")
```