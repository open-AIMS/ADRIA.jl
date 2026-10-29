# A coral block must never skip a size class within a timestep, so the effective linear
# extension of size classes 2 onwards cannot exceed the class width. Class 1 is exempt as
# blocks leaving it move wholesale into class 2, and the terminal class does not grow.
#
# Effective linear extension is `linear_extension * scale * growth_constraint`. The growth
# constraint is at most `max(1, 1 + height)`: growth acceleration lies in (1, 1 + height)
# and the spatial competition factor is <= 1. The worst case therefore combines the upper
# bound of the linear extension, of its scale factor and of the growth acceleration height.

const _GROWTH_LIMIT_TOL = 1e-12

_le_factor_name(fg::AbstractString, fg_idx::Int, size_class::Int) =
    Symbol("$(fg)_$(fg_idx)_$(size_class)_linear_extension")
_scale_factor_name(cb_group::Int, fg::AbstractString) =
    Symbol("linear_extension_scale_cb_group_$(cb_group)_$(fg)")
_height_factor_name(cb_group::Int) =
    Symbol("growth_acceleration_cb_group_$(cb_group)_height")

const _GrowthLimitTerms = NamedTuple{(:le_ub, :scale_ub, :gmax, :width),NTuple{4,Float64}}

_growth_excess(t::_GrowthLimitTerms)::Float64 = t.le_ub * t.scale_ub * t.gmax - t.width

"""
    _growth_limit_terms(spec::DataFrame)::Dict{Tuple{Int,Int,Int},_GrowthLimitTerms}

Worst-case terms of the size class growth limit for each (cb group, functional group, size
class) with size class in `2:(n_sizes - 1)`. Upper bounds are read from `dist_params`.
Combinations whose factors are not in `spec` are skipped.
"""
function _growth_limit_terms(spec::DataFrame)::Dict{Tuple{Int,Int,Int},_GrowthLimitTerms}
    upper = Dict{Symbol,Float64}()
    cb_groups = Int[]
    for (fn, dp) in zip(spec.fieldname, spec.dist_params)
        name = String(fn)
        height_match = match(r"^growth_acceleration_cb_group_(\d+)_height$", name)
        isnothing(height_match) || push!(cb_groups, parse(Int, height_match[1]))

        if endswith(name, "_linear_extension") ||
            startswith(name, "linear_extension_scale_cb_group_") ||
            !isnothing(height_match)
            upper[fn] = Float64(length(dp) >= 2 ? dp[2] : dp[1])
        end
    end
    unique!(cb_groups)

    widths = diff(bin_edges(; unit=:m).data; dims=2)
    fg_names = string.(functional_group_names())

    terms = Dict{Tuple{Int,Int,Int},_GrowthLimitTerms}()
    for cb in cb_groups
        height_ub = get(upper, _height_factor_name(cb), nothing)
        isnothing(height_ub) && continue
        gmax = max(1.0, 1.0 + height_ub)

        for (g, fg) in enumerate(fg_names)
            scale_ub = get(upper, _scale_factor_name(cb, fg), nothing)
            isnothing(scale_ub) && continue

            for s = 2:(size(widths, 2) - 1)
                le_ub = get(upper, _le_factor_name(fg, g, s), nothing)
                isnothing(le_ub) && continue
                terms[(cb, g, s)] = (
                    le_ub=le_ub, scale_ub=scale_ub, gmax=gmax, width=widths[g, s]
                )
            end
        end
    end

    return terms
end

"""
    _scale_ub_caps(terms)::Dict{Tuple{Int,Int},Float64}

Largest upper bound of each `linear_extension_scale` factor, keyed by (cb group, functional
group), that keeps every size class within its width given the other upper bounds.
"""
function _scale_ub_caps(
    terms::Dict{Tuple{Int,Int,Int},_GrowthLimitTerms}
)::Dict{Tuple{Int,Int},Float64}
    caps = Dict{Tuple{Int,Int},Float64}()
    for ((cb, g, _), t) in terms
        t.le_ub > 0.0 || continue
        caps[(cb, g)] = min(get(caps, (cb, g), Inf), t.width / (t.le_ub * t.gmax))
    end

    return caps
end

"""
    _cap_growth_bounds!(model::Model)::Nothing

Lower the upper sampling bound of each `linear_extension_scale` factor (and its default
bound) so the size class growth limit holds in the worst case. A bound is never lowered
below the factor's nominal value, in which case a warning is issued.
"""
function _cap_growth_bounds!(model::Model)::Nothing
    spec = model_spec(model)
    caps = _scale_ub_caps(_growth_limit_terms(spec))
    isempty(caps) && return nothing

    fg_names = string.(functional_group_names())
    row_of = Dict(fn => i for (i, fn) in enumerate(spec.fieldname))
    dist_params = Vector{Float64}[collect(Float64, dp) for dp in spec.dist_params]
    default_dist_params = Vector{Float64}[
        collect(Float64, dp) for dp in spec.default_dist_params
    ]

    n_capped = 0
    above_nominal = Symbol[]
    for ((cb, g), cap) in caps
        fn = _scale_factor_name(cb, fg_names[g])
        i = row_of[fn]
        val = Float64(spec.val[i])
        cap < val && push!(above_nominal, fn)

        capped_here = false
        for dp in (dist_params[i], default_dist_params[i])
            new_ub = max(min(dp[2], cap), val)
            if new_ub < dp[2]
                dp[2] = new_ub
                capped_here = true
            end
        end
        n_capped += capped_here
    end

    if n_capped > 0
        model[:val] = Float64.(spec.val)
        model[:dist_params] = dist_params
        model[:default_dist_params] = default_dist_params
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

"""
    _validate_growth_limits(old_spec::DataFrame, new_spec::DataFrame)::Nothing

Throw an error if `new_spec` breaks the size class growth limit by more than `old_spec`
did. Violations already present in `old_spec` are tolerated so unrelated edits to a model
that was created with infeasible nominal values still go through.
"""
function _validate_growth_limits(old_spec::DataFrame, new_spec::DataFrame)::Nothing
    new_terms = _growth_limit_terms(new_spec)
    isempty(new_terms) && return nothing
    old_terms = _growth_limit_terms(old_spec)

    worse = Tuple{Int,Int,Int}[]
    for (key, t) in new_terms
        excess = _growth_excess(t)
        excess > _GROWTH_LIMIT_TOL || continue
        old_excess = haskey(old_terms, key) ? _growth_excess(old_terms[key]) : 0.0
        excess > old_excess + _GROWTH_LIMIT_TOL && push!(worse, key)
    end
    isempty(worse) && return nothing
    sort!(worse)

    fg_names = string.(functional_group_names())
    caps = _scale_ub_caps(new_terms)
    details = String[]
    for (cb, g, s) in first(worse, 3)
        t = new_terms[(cb, g, s)]
        push!(
            details,
            "  - $(fg_names[g]), cb group $(cb), size class $(s): " *
            "$(round(t.le_ub * 100; digits=2)) cm linear extension x " *
            "$(round(t.scale_ub; digits=3)) scale x " *
            "$(round(t.gmax; digits=3)) growth acceleration = " *
            "$(round(t.le_ub * t.scale_ub * t.gmax * 100; digits=2)) cm " *
            "> class width of $(round(t.width * 100; digits=2)) cm. " *
            "Upper bound of `$(_scale_factor_name(cb, fg_names[g]))` must be at most " *
            "$(round(caps[(cb, g)]; digits=3))."
        )
    end
    more = length(worse) > 3 ? "\n  ... and $(length(worse) - 3) more." : ""

    return error(
        "Sampling bounds would let coral blocks skip a size class " *
        "($(length(worse)) case(s)). The worst-case effective linear extension " *
        "(upper bounds of the linear extension, its scale factor and growth acceleration " *
        "height) must not exceed the size class width.\n" *
        join(details, "\n") * more * "\n" *
        "Lower the upper bound of the scale factor, or of the linear extension or " *
        "growth acceleration height factors."
    )
end
