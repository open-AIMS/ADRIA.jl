using ADRIAinterfaces: ParamSpec

# `ParamSpec` kind => `ptype` of the parameter table
const _KIND_PTYPE = (
    real="continuous", integer="ordered discrete", bool="unordered categorical",
    categorical="unordered categorical"
)

# `ParamSpec` dist symbol => (`dist` of the parameter table, number of `dist_params`)
function _spec_dist(spec::ParamSpec)
    d = spec.dist
    d === :uniform && return Uniform, 2
    d === :triang && return TriangularDist, 3
    d === :discrete_ordered && return DiscreteOrderedUniformDist, 3
    d === :discrete_uniform && return DiscreteUniform, 2
    d === :categorical && return CategoricalDistribution, length(spec.levels)
    throw(ArgumentError("$(spec.name): unknown dist :$d"))
end

"""
    spec_rows(
        specs::AbstractVector{ParamSpec}; overrides=Dict{String,Float64}()
    )::Vector{Pair{Symbol,NamedTuple}}

Turn `ParamSpec`s (declarations by submodels, see `ADRIAinterfaces.ParamSpec`) into the
factor specs (in order) that the parameter table (`SpecStore`) is built from, as
`name => NamedTuple` with the `val`, `ptype`, `dist`, `dist_params`, `default_dist_params`,
`name` and `description` of `_factor_spec`. Lives in ADRIA because the table columns hold
ADRIA's distribution types, which the plain-data `ParamSpec` must not know about.
`overrides` (keyed by factor name) replaces defaults only, bounds are unaffected.

# Mapping to table columns

| `ParamSpec`                  | column             | value                                    |
|:-----------------------------|:-------------------|:-----------------------------------------|
| `name`                       | `fieldname`        | as is                                    |
| `label`                      | `name`             | as is                                    |
| `description`                | `description`      | as is                                    |
| `default` (or override)      | `val`              | as is                                    |
| `kind` `:real`               | `ptype`            | `"continuous"`                           |
| `kind` `:integer`            | `ptype`            | `"ordered discrete"`                     |
| `kind` `:bool`/`:categorical`| `ptype`            | `"unordered categorical"`                |
| `dist` `:uniform`            | `dist`             | `Uniform`, params (lo, hi)               |
| `dist` `:triang`             | `dist`             | `TriangularDist`, params (lo, hi, mode)  |
| `dist` `:discrete_ordered`   | `dist`             | `DiscreteOrderedUniformDist`, (lo, hi, step) |
| `dist` `:discrete_uniform`   | `dist`             | `DiscreteUniform`, params (lo, hi)       |
| `dist` `:categorical`        | `dist`             | `CategoricalDistribution`, params `levels` |

`dist_params` of the `ParamSpec` (`levels` for `:categorical`) is the factor's full range.
The constant-by-default rule: with `constant=true` (the default) the table's `dist_params`
collapse to the default value (`fill(default, n)`, as `fix_factor!` does) so `sample()` does
not vary the factor, and the full range is kept as `default_dist_params` (what
`set_factor_bounds!` can widen to). With `constant=false`, `dist_params` equals
`default_dist_params`, the full range. The number of parameters must match the dist
(2, 3, 3, 2, `length(levels)` as above), and a `kind` of `:bool`/`:categorical` needs the
`:categorical` dist (and vice versa).
"""
function spec_rows(
    specs::AbstractVector{ParamSpec};
    overrides::AbstractDict{String,Float64}=Dict{String,Float64}()
)::Vector{Pair{Symbol,NamedTuple}}
    rows = Pair{Symbol,NamedTuple}[]
    for s in specs
        dist, n = _spec_dist(s)
        (s.kind in (:bool, :categorical)) == (s.dist === :categorical) || throw(
            ArgumentError("$(s.name): kind :$(s.kind) does not match dist :$(s.dist)")
        )
        full = s.dist === :categorical ? s.levels : s.dist_params
        length(full) == n || throw(
            ArgumentError(
                "$(s.name): dist :$(s.dist) needs $n dist_params, got $(length(full))"
            )
        )
        val = get(overrides, String(s.name), s.default)
        push!(
            rows,
            s.name => _factor_spec(
                val;
                ptype=getproperty(_KIND_PTYPE, s.kind),
                dist=dist,
                dist_params=s.constant ? fill(Float64(val), n) : full,
                default_dist_params=full,
                name=s.label,
                description=s.description
            )
        )
    end
    return rows
end
