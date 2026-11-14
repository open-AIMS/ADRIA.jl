"""
Flatten.jl / ModelParameters concrete-method overrides for large `Param`-only structs.

## Problem

`Flatten._flatten`, `_reconstruct`, and `metaflatten` are `@generated` functions that
build a per-field inline expression for every struct they encounter.  For structs with
many homogeneous fields (`Coral` has 331, `GrowthAcceleration` has 36) the expression
size grows as O(N), but Julia must *type-infer* each intermediate tuple concatenation step,
so the total compile cost is O(N²).  For the domain model tuple assembled
in `_assemble_domain_model`, the outer `@generated _reconstruct` additionally tries to
inline all ~440 fields at once, giving O(440²) compile cost (the regression this
caused was minutes of cold JIT; 357 s measured by the precompile budget test).

## Solution

Julia dispatches to the *most specific* applicable method first, so adding a concrete
(non-generated) method for a specific type completely bypasses the `@generated` path with
no correctness cost.  This file defines five categories of override:

- **Fix A** — `_flatten(::Tuple{Vararg{AbstractParam}}, ...)`: identity with
  `@nospecialize`, so the return type is always `Tuple{Vararg{AbstractParam}}` (not the
  concrete N-element form), preventing O(N²) specialisation of downstream callers.

- **Fix B** — `params(::Coral)` / `params(::GrowthAcceleration)`: direct accessors that
  pin the statically-inferred return type to `Tuple{Vararg{Param}}` (via
  `Tuple(::Vector{Param})`), ensuring `_expandkeys` compiles for the variadic form.

- **Fix C** — `_expandkeys(::Tuple{Vararg{AbstractParam}})` with `@nospecialize`: one
  compilation covers all concrete subtypes; the fast path detects homogeneous key sets
  (the common case) and returns early.

- **Fix D** — `_reconstruct(::Tuple, ...)` with `@nospecialize`: replaces the
  `@generated` outer-tuple reconstruct with a simple loop that dispatches `_reconstruct`
  per component, so compile cost is O(1) per component rather than O(total_fields²).

- **Struct-specific overrides** for `Coral` and `GrowthAcceleration`: concrete
  `_flatten`, `_reconstruct`, and `metaflatten` methods that use a `Vector` accumulator
  loop instead of the generated field-by-field expression.

## Placement

This file is included from `ADRIA.jl` immediately after `GrowthAcceleration.jl`, so both
`Coral` and `GrowthAcceleration` are fully defined before these methods are compiled.

## Note on `Base.keys` / `Base.get`

We do NOT redefine `Base.keys` / `Base.get` on `AbstractParam` here.
`ModelParameters` already provides those definitions (`param.jl:45`).  Redefining them in
ADRIA triggers "Method overwriting is not permitted during Module precompilation", which
prevents ADRIA from caching its precompile artefacts and adds ~55 s to every
`using ADRIA` session.
"""

import Flatten as _Flatten
import ModelParameters: AbstractParam
import ModelParameters as _MP

# ---------------------------------------------------------------------------
# Type-erased tuple -> Vector{Any}.  `collect(Any, t)` specialises on the concrete
# (~440-element, heterogeneous) tuple type, ~2.5 s cold per distinct `Param` mix.
# `@nospecialize` compiles this once; `t[i]` with a runtime index is a plain
# dynamic lookup.
# ---------------------------------------------------------------------------
function _tuple_to_vec(@nospecialize(t::Tuple))
    n = length(t)
    v = Vector{Any}(undef, n)
    for i = 1:n
        v[i] = t[i]
    end
    return v
end

# ---------------------------------------------------------------------------
# Fix A: @nospecialize prevents Julia from specialising this identity function
# on each concrete subtype of Tuple{Vararg{AbstractParam}}.  The return type
# is therefore Tuple{Vararg{AbstractParam}} regardless of input, which makes
# map() inside _expandkeys compile as a generic loop (O(1)) rather than a
# 331-element unrolled specialisation (O(N²)).
# ---------------------------------------------------------------------------
function _Flatten._flatten(
    @nospecialize(params::Tuple{Vararg{AbstractParam}}), ::Function, ::Type, ::Type
)
    return params   # already flat -- identity function, no generated expression needed
end

# ---------------------------------------------------------------------------
# Fix C: _expandkeys override with @nospecialize.
#
# Even with Fix A+B, _expandkeys is dispatched at RUNTIME on the concrete
# NTuple{N,Param{...}} type (since that is what Tuple(v::Vector{Param})
# produces).  Without @nospecialize, Julia compiles a specialisation of
# _expandkeys for NTuple{N,Param{...}} which unrolls map() N times -> O(N²).
#
# With @nospecialize, Julia compiles ONE version for all
# Tuple{Vararg{AbstractParam}} subtypes.  Inside that single body, pars is
# typed as Tuple{Vararg{AbstractParam}} (unknown length), so map() and all()
# compile as generic loops -> O(1).
#
# Correctness: _expandkeys gives every Param the union of all keys across all
# params (filling absent keys with nothing).  For Coral and GrowthAcceleration
# every Param already carries the same 7 keys, so _expandkeys is a semantic
# no-op.  The fast path detects this and returns pars unchanged.  The slow
# path handles structs with heterogeneous Params.
# ---------------------------------------------------------------------------
function _MP._expandkeys(@nospecialize(pars::Tuple{Vararg{AbstractParam}}))
    isempty(pars) && return pars
    # Iterate a type-erased `Vector`, never the tuple itself: `enumerate(pars)` /
    # `for p in pars` over the runtime-concrete ~440-element tuple type were each
    # JIT-specialised on it (~20 s cold).
    pv = _tuple_to_vec(pars)
    k1 = keys(pv[1])
    # Fast path: all params already share the same key set.
    uniform = true
    for i in eachindex(pv)
        if keys(pv[i]) != k1
            uniform = false
            break
        end
    end
    uniform && return pars
    # Slow path: heterogeneous key sets (e.g. only some params carry `:direction`).
    # Plain loops over `Vector`s, not splatted `union(map(keys, pars)...)` / `map` over the
    # tuple, which compile a ~440-way specialisation (>10 s cold).  Key order matches the
    # original `union` (first-seen order).
    allkeys_v = Symbol[]
    for i in eachindex(pv)
        for k in keys(pv[i])
            k in allkeys_v || push!(allkeys_v, k)
        end
    end
    allkeys = Tuple(allkeys_v)
    out = Vector{Any}(undef, length(pv))
    for i in eachindex(pv)
        par = pv[i]
        vals = Tuple(get(par, key, nothing) for key in allkeys_v)
        out[i] = _MP.rebuild(par, NamedTuple{allkeys}(vals))
    end
    return Tuple(out)
end

# ---------------------------------------------------------------------------
# Coral-specific overrides
# ---------------------------------------------------------------------------

# Fix B (Coral): pins the statically-inferred return type of params(coral) to
# Tuple{Vararg{Param}} so _expandkeys compiles for the variadic form (O(1))
# rather than the 331-element concrete NTuple form (O(N²)).
function _MP.params(coral::Coral)
    return _Flatten._flatten(coral, _Flatten.flattenable, _MP.SELECT, _MP.IGNORE)
end

function _Flatten._flatten(coral::Coral, ::Function, ::Type, ::Type)
    fnames = fieldnames(Coral)
    N = length(fnames)
    v = Vector{Param}(undef, N)
    @inbounds for i = 1:N
        v[i] = getfield(coral, fnames[i])
    end
    return Tuple(v)   # Tuple{Vararg{Param}} — avoids @generated 331-element expression
end

function _Flatten._reconstruct(::Coral, data, ::Function, ::Type, ::Type, n::Int)
    N = fieldcount(Coral)
    v = Vector{Any}(undef, N)
    @inbounds for i = 1:N
        v[i] = data[n + i - 1]
    end
    return (Coral(v...), n + N)   # positional _apply path — O(1) compile cost
end

function _Flatten.metaflatten(
    ::Coral, func::Function, ::Function, ::Type, ::Type, P, ::Any
)
    fnames = fieldnames(Coral)
    return Tuple(func(Coral, Val{fnames[i]}) for i in eachindex(fnames))
end

# ---------------------------------------------------------------------------
# GrowthAcceleration-specific overrides
# ---------------------------------------------------------------------------

# Fix B (GrowthAcceleration): same rationale as params(coral::Coral) above.
function _MP.params(ga::GrowthAcceleration)
    return _Flatten._flatten(ga, _Flatten.flattenable, _MP.SELECT, _MP.IGNORE)
end

function _Flatten._flatten(ga::GrowthAcceleration, ::Function, ::Type, ::Type)
    fnames = fieldnames(GrowthAcceleration)
    N = length(fnames)
    v = Vector{Param}(undef, N)
    @inbounds for i = 1:N
        v[i] = getfield(ga, fnames[i])
    end
    return Tuple(v)
end

function _Flatten._reconstruct(
    ::GrowthAcceleration, data, ::Function, ::Type, ::Type, n::Int
)
    N = fieldcount(GrowthAcceleration)
    v = Vector{Any}(undef, N)
    @inbounds for i = 1:N
        v[i] = data[n + i - 1]
    end
    return (GrowthAcceleration(v...), n + N)   # positional _apply path — O(1) compile cost
end

function _Flatten.metaflatten(
    ::GrowthAcceleration, func::Function, ::Function, ::Type, ::Type, P, ::Any
)
    fnames = fieldnames(GrowthAcceleration)
    return Tuple(func(GrowthAcceleration, Val{fnames[i]}) for i in eachindex(fnames))
end

# ---------------------------------------------------------------------------
# Fix D: Override Flatten._reconstruct for outer Tuple types.
#
# The @generated _reconstruct inlines code for ALL fields of ALL components in
# the tuple (the domain model components together hold ~450 params) ->
# O(440²) compile time (minutes cold).  This override loops over tuple elements and
# dispatches _reconstruct per component, so each component is compiled
# independently in O(1) time via its own concrete override above.
# ---------------------------------------------------------------------------
function _Flatten._reconstruct(
    @nospecialize(obj::Tuple),
    @nospecialize(data),
    ft::Function,
    ::Type{U},
    ::Type{I},
    n::Int
) where {U,I}
    len = length(obj)
    results = Vector{Any}(undef, len)
    # Type-erase `data`.  Components without a concrete override above (`Intervention`,
    # the criteria weights, `EnvironmentalLayer`, ...) still go through the `@generated`
    # `_reconstruct`, which specialises on the type of `data`.  Handing it the full
    # heterogeneous ~440-element tuple made each of those compile against that giant type
    # (~170 s cold, mostly `Intervention`).  A `Vector{Any}` gives one cheap, shared
    # specialisation.
    # Only convert when needed; anything else is passed through unchanged (the
    # `@generated` path accepts it exactly as before).
    vdata = if data isa Vector{Any}
        data
    elseif data isa Tuple
        _tuple_to_vec(data)
    elseif data isa AbstractVector
        Vector{Any}(data)
    else
        data
    end
    @inbounds for i = 1:len
        results[i], n = _Flatten._reconstruct(obj[i], vdata, ft, U, I, n)
    end
    return (Tuple(results), n)
end
