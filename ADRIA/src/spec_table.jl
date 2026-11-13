"""
Table-backed parameter store of a Domain.

A Domain holds one `SpecStore`: a plain `DataFrame` with one row per model factor, in
component order (`EnvironmentalLayer`, `Intervention`, `CAqCriteriaWeights`,
`FogCriteriaWeights`, `LvMCriteriaWeights`, `DepthThresholds`, `Coral`,
`GrowthAcceleration`, `DepthAttenuation`, `LowDHWBleaching`, `CotsParams`).

The table is the source of truth for factor values and bounds. It replaces the
`ModelParameters.Model` the Domain used to hold: every `Model` operation (`DataFrame(model)`,
`model[:col] = ...`) specialises on the ~440-element heterogeneous `Param` tuple and cost
tens of seconds of JIT per Julia session, whereas operations on a `DataFrame` compile once.

# Where the rows come from
Rows are derived from the component definitions, never copied into a static data file, so a
changed default cannot silently diverge from the table:

- `Coral` and `GrowthAcceleration` from `_build_coral_fields` / `_build_growth_acceleration_fields`
  (the same builders that generate their structs),
- `DepthAttenuation`, `LowDHWBleaching` and `CotsParams` from their `_*_specs` builders (the
  same factor definitions their struct defaults use),
- the remaining components by reading the `Param` fields of a default instance (`Intervention()`,
  `CAqCriteriaWeights()`, ..., and `EnvironmentalLayer(dhw, wave, cyclone)` /
  `Intervention(; iv_MCB_...)` whose factor ranges depend on the domain data).

Calibration overrides (`CalibOverrides`) are applied by the builders, so calibrated rows are
built directly, without constructing component instances.

# Columns
`component` (String), `fieldname` (Symbol), `val` (Float64), `ptype` (String), `dist`
(Type), `dist_params` and `default_dist_params` (`Vector{Float64}`), `direction`
(`Union{Nothing,Function}`; criteria weights only), `name` and `description` (String).
`model_spec` adds the derived `lower_bound`, `upper_bound` and `is_constant` columns.
"""
mutable struct SpecStore
    tbl::DataFrame
    row_of::Dict{Symbol,Int}
end

"""
    copy(s::SpecStore)::SpecStore

Independent copy of the table. `update!`/`update_params!` write into the table in place, so a
shallow copy of a Domain (`copy`/`Setfield` on the struct) shares its table with the original
and writes through either one are seen by both; use `deepcopy(dom)` (which copies the table)
to get an independent Domain.
"""
Base.copy(s::SpecStore)::SpecStore = SpecStore(deepcopy(s.tbl), copy(s.row_of))

"""
    spec_components()::NTuple{11,DataType}

Components of the parameter table, in row order. The single definition of which components
a domain's table holds (`build_spec_store` and the tests use it). A function, not a constant,
because some of the component types are defined after this file is included.
"""
spec_components() = (
    EnvironmentalLayer, Intervention, CAqCriteriaWeights, FogCriteriaWeights,
    LvMCriteriaWeights, DepthThresholds, Coral, GrowthAcceleration, DepthAttenuation,
    LowDHWBleaching, CotsParams
)

_component_name(component)::String = replace(string(component), "ADRIA." => "")

"""Factor spec (see `_factor_spec`) held by a `Param`."""
function _param_spec(p::Param)::NamedTuple
    return (; (k => getproperty(p, k) for k in keys(p))...)
end

"""Factor specs of the `Param` fields of a component instance, in field order."""
function _instance_specs(inst)::Vector{Pair{Symbol,NamedTuple}}
    return Pair{Symbol,NamedTuple}[
        fn => _param_spec(getfield(inst, fn)) for fn in fieldnames(typeof(inst))
    ]
end

_dict_specs(d::AbstractDict)::Vector{Pair{Symbol,NamedTuple}} =
    Pair{Symbol,NamedTuple}[Symbol(k) => v for (k, v) in d]

function _component_table(component, specs::Vector{Pair{Symbol,NamedTuple}})::DataFrame
    n = length(specs)
    return DataFrame(
        :component => fill(_component_name(component), n),
        :fieldname => Symbol[fn for (fn, _) in specs],
        :val => Float64[Float64(sp.val) for (_, sp) in specs],
        :ptype => String[sp.ptype for (_, sp) in specs],
        :dist => Type[sp.dist for (_, sp) in specs],
        :dist_params => Vector{Float64}[collect(Float64, sp.dist_params) for (_, sp) in specs],
        :default_dist_params => Vector{Float64}[
            collect(Float64, sp.default_dist_params) for (_, sp) in specs
        ],
        :direction => Union{Nothing,Function}[get(sp, :direction, nothing) for (_, sp) in specs],
        :name => String[sp.name for (_, sp) in specs],
        :description => String[sp.description for (_, sp) in specs];
        copycols=false
    )
end

"""
    build_spec_store(
        env::EnvironmentalLayer, intervention::Intervention;
        calib::CalibOverrides=CalibOverrides()
    )::SpecStore

Assemble the parameter table of a domain from the component definitions (see `SpecStore`).
`env` and `intervention` carry the data dependent factor ranges (number of DHW, wave and
cyclone scenarios; MCB albedo and duration levels). Calibration `calib` overrides default
`Coral`, `GrowthAcceleration`, `DepthAttenuation`, `LowDHWBleaching` and `CotsParams` values.

The growth limit caps (`_cap_growth_bounds!`) are applied.
"""
function build_spec_store(
    env::EnvironmentalLayer,
    intervention::Intervention;
    calib::CalibOverrides=CalibOverrides()
)::SpecStore
    specs = Dict{DataType,Vector{Pair{Symbol,NamedTuple}}}(
        EnvironmentalLayer => _instance_specs(env),
        Intervention => _instance_specs(intervention),
        CAqCriteriaWeights => _instance_specs(CAqCriteriaWeights()),
        FogCriteriaWeights => _instance_specs(FogCriteriaWeights()),
        LvMCriteriaWeights => _instance_specs(LvMCriteriaWeights()),
        DepthThresholds => _instance_specs(DepthThresholds()),
        Coral => _dict_specs(_build_coral_fields(; overrides=calib.coral)),
        GrowthAcceleration => _dict_specs(
            _build_growth_acceleration_fields(; overrides=calib.growth_accel)
        ),
        DepthAttenuation => _eff_dhw_specs(; overrides=calib.depth_atten),
        LowDHWBleaching => _low_dhw_bleach_specs(; overrides=calib.low_dhw),
        CotsParams => _cots_specs(; overrides=calib.cots)
    )
    # Every component must be listed in `spec_components()` and vice versa
    Set(keys(specs)) == Set(spec_components()) || error(
        "`spec_components()` and the components built in `build_spec_store` differ: " *
        "$(symdiff(Set(keys(specs)), Set(spec_components())))"
    )
    parts = DataFrame[_component_table(c, specs[c]) for c in spec_components()]
    tbl = reduce(vcat, parts)
    row_of = Dict{Symbol,Int}(f => i for (i, f) in enumerate(tbl.fieldname))
    length(row_of) == nrow(tbl) ||
        error("Factor names must be unique across the components of a domain.")

    store = SpecStore(tbl, row_of)
    _cap_growth_bounds!(store)
    return store
end

# Same post-processing as the former `model_spec(::Model)`, applied to the table (or a
# slice of it). Tuple columns use the eltype joined over the *full* table's lengths, so a
# component slice has the same eltype as the same rows taken from the full spec.
function _tuple_col(v::Vector{Vector{Float64}}, ref::Vector{Vector{Float64}})
    T = mapreduce(n -> NTuple{n,Float64}, typejoin, unique(length.(ref)))
    return T[Tuple(x) for x in v]
end

function _spec_shape(t::DataFrame, full::DataFrame=t)::DataFrame
    spec = DataFrame(t; copycols=true)
    spec[!, :dist_params] = _tuple_col(spec.dist_params, full.dist_params)
    spec[!, :default_dist_params] = _tuple_col(
        spec.default_dist_params, full.default_dist_params
    )
    DataFrames.hcat!(
        spec,
        DataFrame(
            :lower_bound => distribution_lower_bound.(spec.dist, spec.dist_params),
            :upper_bound => distribution_upper_bound.(spec.dist, spec.dist_params)
        )
    )
    spec[!, :is_constant] .= spec[!, :lower_bound] .== spec[!, :upper_bound]

    # Reorder so name/description appears at end
    # makes viewing as CSV a little nicer given description can be very long
    select!(spec, Not([:name, :description]), [:name, :description])

    return spec
end

"""
    model_spec(s::SpecStore)::DataFrame

Model specification as a DataFrame with lower and upper bounds (a copy: modify it freely and
write it back with `update!`).
"""
model_spec(s::SpecStore)::DataFrame = _spec_shape(s.tbl)

"""
    component_params(s::SpecStore, component)::DataFrame
    component_params(s::SpecStore, components::Vector{T})::DataFrame

Rows of one or several components, in the form of `model_spec`. Only the requested rows are
shaped, so this is cheap enough to call per scenario.
"""
function component_params(s::SpecStore, component)::DataFrame
    return _spec_shape(s.tbl[s.tbl.component .== _component_name(component), :], s.tbl)
end
function component_params(s::SpecStore, components::Vector{T})::DataFrame where {T}
    names = _component_name.(components)
    return _spec_shape(s.tbl[s.tbl.component .∈ [names], :], s.tbl)
end

"""
    param_table(s::SpecStore)::DataFrame

One-row DataFrame of the nominal factor values (`fieldname => val`).
"""
function param_table(s::SpecStore)::DataFrame
    return DataFrame(OrderedDict(string(k) => v for (k, v) in zip(s.tbl.fieldname, s.tbl.val)))
end

function _rows_of(s::SpecStore, fieldnames)::Vector{Int}
    return Int[
        get(s.row_of, Symbol(f)) do
            throw(ArgumentError("Unknown factor `$f`"))
        end for f in fieldnames
    ]
end

"""
    update!(s::SpecStore, spec::DataFrame)::Nothing

Write `val` (as `Float64`) and `dist_params` (as `Vector{Float64}`) of `spec` into the table.
Rows are matched by `fieldname` (so `spec` may hold a subset of the factors); `spec` must
have that column.

Throws an error, leaving the table untouched, if `spec` holds an unknown factor or if its
bounds would let a coral block skip a size class more than the current ones do (see
`_validate_growth_limits`).
"""
function update!(s::SpecStore, spec::DataFrame)::Nothing
    hasproperty(spec, :fieldname) ||
        throw(ArgumentError("`spec` needs a `fieldname` column to match factors by name"))
    rows = _rows_of(s, spec.fieldname)
    _validate_growth_limits(s.tbl, spec)

    vals = Float64.(spec.val)
    dist_params = Vector{Float64}[collect(Float64, dp) for dp in spec.dist_params]
    for (k, i) in enumerate(rows)
        s.tbl.val[i] = vals[k]
        s.tbl.dist_params[i] = dist_params[k]
    end

    return nothing
end

"""
    update_params!(s::SpecStore, params::DataFrameRow)::Nothing

Set the `val` of the factors named by the columns of `params` (e.g. a row of a scenario
set). Values are matched by name; a column named `RCP` is ignored and any other name that is
not a factor of the table is an error (nothing is written in that case).
"""
function update_params!(s::SpecStore, params::DataFrameRow)::Nothing
    nms = String[nm for nm in names(params) if nm != "RCP"]
    rows = _rows_of(s, nms)
    for (nm, i) in zip(nms, rows)
        s.tbl.val[i] = Float64(params[nm])
    end

    return nothing
end

"""
    _cap_growth_bounds!(s::SpecStore)::Nothing

Lower the upper sampling bound of each `linear_extension_scale` factor (and its default
bound) so the size class growth limit holds in the worst case. A bound is never lowered
below the factor's nominal value, in which case a warning is issued.
"""
function _cap_growth_bounds!(s::SpecStore)::Nothing
    tbl = s.tbl
    caps = _scale_ub_caps(_growth_limit_terms(tbl))
    isempty(caps) && return nothing

    fg_names = string.(functional_group_names())
    n_capped = 0
    above_nominal = Symbol[]
    for ((cb, g), cap) in caps
        fn = _scale_factor_name(cb, fg_names[g])
        i = s.row_of[fn]
        val = tbl.val[i]
        cap < val && push!(above_nominal, fn)

        capped_here = false
        for dp in (tbl.dist_params[i], tbl.default_dist_params[i])
            new_ub = max(min(dp[2], cap), val)
            if new_ub < dp[2]
                dp[2] = new_ub
                capped_here = true
            end
        end
        n_capped += capped_here
    end

    if n_capped > 0
        @info "Capped $(n_capped) `linear_extension_scale` upper bounds so growth cannot " *
            "exceed size class widths."
    end

    if !isempty(above_nominal)
        sort!(above_nominal)
        @warn "The worst-case growth limit is below the nominal value of " *
            "$(length(above_nominal)) `linear_extension_scale` factors, e.g. " *
            "$(first(above_nominal)). Their upper bounds were set to the nominal value; " *
            "the linear extension, scale or growth acceleration height values may need " *
            "recalibrating for the current size class edges."
    end

    return nothing
end
