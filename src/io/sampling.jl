using Printf
using DataFrames, Distributions, LinearAlgebra

using ADRIA
using ADRIA: model_spec, component_params
using ADRIA.decision: mcda_normalize
import Distributions: sample
import QuasiMonteCarlo as QMC
import QuasiMonteCarlo: SobolSample, OwenScramble

import OrderedCollections: OrderedDict

const DISCRETE_FACTOR_TYPES = [
    "ordered categorical", "unordered categorical", "ordered discrete"
]

"""
    adjust_samples(d::Domain, df::DataFrame)::DataFrame
    adjust_samples!(spec::DataFrame, df::DataFrame)::DataFrame

Adjust given samples to ensure parameter value combinations for unguided
scenarios are plausible.
"""
function adjust_samples(d::Domain, df::DataFrame)::DataFrame
    return adjust_samples(model_spec(d), df)
end
function adjust_samples(spec::DataFrame, df::DataFrame)::DataFrame
    interv = component_params(spec, Intervention)
    seed_weights = component_params(spec, SeedCriteriaWeights)
    fog_weights = component_params(spec, FogCriteriaWeights)
    depth_offsets = component_params(spec, DepthThresholds)

    # If counterfactual, set all intervention options to 0.0
    df[df.guided .== -1.0, filter(x -> x ∉ [:guided, :heritability], interv.fieldname)] .=
        0.0

    # If unguided/counterfactual, set all preference criteria, except those related to depth, to 0.
    non_depth_names = vcat(
        seed_weights.fieldname,
        fog_weights.fieldname
    )
    df[df.guided .== 0.0, non_depth_names] .= 0.0  # Turn off weights for unguided
    df[df.guided .== -1.0, non_depth_names] .= 0.0  # Turn off weights for cf
    df[df.guided .== -1.0, depth_offsets.fieldname] .= 0.0  # No depth offsets for cf

    # If unguided, set planning horizon to 0.
    df[df.guided .== 0.0, :plan_horizon] .= 0.0

    # If no seeding is to occur, set related variables to 0
    not_seeded = (df.N_seed_TA .== 0) .& (df.N_seed_CA .== 0) .& (df.N_seed_SM .== 0)
    df[not_seeded, contains.(names(df), "seed_")] .= 0.0
    df[not_seeded, :a_adapt] .= 0.0

    # Same for fogging/shading
    not_fogged = (df.fogging .== 0)
    df[not_fogged, contains.(names(df), "fog_")] .= 0.0

    not_shaded = (df.SRM .== 0)
    df[not_shaded, contains.(names(df), "shade_")] .= 0.0

    # Normalize MCDA weights for fogging scenarios
    guided_fogged = (df.fogging .> 0.0) .& (df.guided .> 0)
    df[guided_fogged, fog_weights.fieldname] .= mcda_normalize(
        df[guided_fogged, fog_weights.fieldname]
    )
    # Normalize MCDA weights for seeding scenarios
    guided_seeded = .!(not_seeded) .& (df.guided .> 0)
    df[guided_seeded, seed_weights.fieldname] .= mcda_normalize(
        df[guided_seeded, seed_weights.fieldname]
    )

    if nrow(unique(df)) < nrow(df)
        perc = "$(@sprintf("%.3f", (1.0 - (nrow(unique(df)) / nrow(df))) * 100.0))%"
        @warn "Non-unique samples created: $perc of the samples are duplicates."
    end

    return df
end

"""
    sample(dom::Domain, n::Int, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Create samples and rescale to distribution defined in the model spec.

Notes:
- assumes all parameters are independent.

# Arguments
- `dom` : Domain
- `n` : Int
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample(
    dom::Domain, n::Int; sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    n > 0 ? n : throw(DomainError(n, "`n` must be > 0"))
    ms::DataFrame = model_spec(dom)
    cb_calib_groups::Vector{Int64} = dom.loc_data.CB_CALIB_GROUPS
    return sample(_filtered_model_spec(ms, cb_calib_groups), n, sample_method)
end

"""
    _filtered_model_spec(model_spec::DataFrame, cb_calib_groups::Vector{Int64})

Filters `cb_calib_group` specific factors correspondent to `cb_calib_groups` not present in
loaded domain spatial data.
"""
function _filtered_model_spec(model_spec::DataFrame, cb_calib_groups::Vector{Int64})
    factor_names = string.(model_spec.fieldname)

    # The extra "_" garantees that we don't select group "12", for example, when trying to
    # select just group "1"
    cb_calib_group_ids = string.(unique(cb_calib_groups))
    groups_to_keep = "cb_group_" .* cb_calib_group_ids .* "_"

    # Factors to remove
    group_mask =
        .!reduce.(
            (x, y) -> x .|| y, [occursin.(groups_to_keep, fname) for fname in factor_names]
        )
    growth_acc_mask = occursin.(r"growth_acceleration", factor_names) .&& group_mask
    mb_rate_mask = occursin.(r"mb_rate_scale", factor_names) .&& group_mask
    linear_extension_mask =
        occursin.(r"linear_extension_scale", factor_names) .&& group_mask

    # Factors to keep
    keep_mask = .!growth_acc_mask .&& .!mb_rate_mask .&& .!linear_extension_mask

    return keepat!(model_spec, keep_mask)
end

"""
    sample(dom::Domain, n::Int, component::Type, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Create samples and rescale to distribution defined in the model spec.

Notes:
- assumes all parameters are independent.

# Arguments
- `dom` : Domain
- `n` : Number of samples to create
- `component` : Component type, e.g. Intervention, Coral, etc.
- `sample_method` : sample_method to use

# Returns
Scenario specification
"""
function sample(
    dom::Domain,
    n::Int,
    component::Type,
    sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    n > 0 ? n : throw(DomainError(n, "`n` must be > 0"))

    spec = component_params(dom.model, component)
    return sample(spec, n, sample_method)
end

"""
    sample(spec::DataFrame, n::Int, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Create samples and rescale to distribution defined in the model spec.

# Arguments
- `spec` : DataFrame containing model parameter specifications.
- `n` : number of samples to generate.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample(
    spec::DataFrame,
    n::Int64,
    sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    if contains(string(sample_method), "SobolSample") && !ispow2(n)
        throw(DomainError(n, "`n` must be a power of 2 when using the Sobol' sampler"))
    end

    # Select non-constant params
    vary_vars = spec[spec.is_constant .== false, [:fieldname, :dist, :dist_params]]
    n_vary_params = size(vary_vars, 1)
    if n_vary_params == 0
        throw(DomainError(n_vary_params, "Number of parameters to perturb must be > 0"))
    end

    # Create distribution types
    vary_dists = try
        map(x -> x.dist(x.dist_params...), eachrow(vary_vars))
    catch err
        err_msg = "Some dist_params could not be converted to integer"

        isa(err, InexactError) ? error(err_msg) : rethrow(err)
    end

    # Create uniformly distributed samples for uncertain parameters
    samples = QMC.sample(n, zeros(n_vary_params), ones(n_vary_params), sample_method)

    # Scale uniform samples to indicated distributions using the inverse CDF method
    samples = Matrix(Distributions.quantile.(vary_dists, samples)')

    # Combine varying and constant values (constant params use their indicated default vals)
    full_df = hcat(fill.(spec.val, n)...)
    full_df[:, spec.is_constant .== false] .= samples

    # Ensure unguided scenarios do not have superfluous factor combinations
    return adjust_samples(spec, DataFrame(full_df, spec.fieldname))
end

"""
    sample_selection(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Create guided samples of factors relevant to location selection.
Coral factors are set to their default values and are not perturbed or sampled.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_selection(
    d::Domain, n::Int64, sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    subset_spec = component_params(
        d.model,
        [
            EnvironmentalLayer,
            Intervention,
            SeedCriteriaWeights,
            FogCriteriaWeights,
            DepthThresholds
        ]
    )

    # Only sample guided intervention scenarios
    _adjust_guided_lower_bound!(subset_spec, 1)

    # Create and fill scenario spec
    # Only Intervention, EnvironmentalLayer and CriteriaWeights factors are perturbed,
    # all other factors are fixed to their default values
    scens = repeat(param_table(d), n)
    scens[:, subset_spec.fieldname] .= sample(subset_spec, n, sample_method)

    return scens
end

"""
    sample_cf(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Generate only counterfactual scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_cf(
    d::Domain, n::Int64, sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    cb_calib_groups::Vector{Int64} = d.loc_data.CB_CALIB_GROUPS
    spec_df = _filtered_model_spec(model_spec(d), cb_calib_groups)

    # Unguided scenarios only
    guided_col = spec_df.fieldname .== :guided
    spec_df[guided_col, [:val, :lower_bound, :upper_bound, :dist_params, :is_constant]] .=
        [-1 -1 -1 (-1.0, -1.0) true]

    # Remove intervention scenarios as an option
    deactivate_interventions(spec_df)

    return sample(spec_df, n, sample_method)
end

"""
    _adjust_guided_lower_bound!(guided_spec::DataFrame, lower::Int64)::DataFrame

Adjust lower bound of guided parameter spec to alter sampling range.
"""
function _adjust_guided_lower_bound!(spec_df::DataFrame, lower::Int64)::DataFrame
    guided_col = spec_df.fieldname .== :guided
    g_upper = Float64(spec_df[guided_col, :upper_bound][1])

    # Update entries, standardizing values for bounds as floats
    spec_df[guided_col, [:val, :lower_bound, :dist_params]] .=
        [lower Float64(lower) (Float64(lower), g_upper)]

    return spec_df
end

"""
    sample_guided(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Generate only guided scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_guided(
    d::Domain, n::Int64, sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    cb_calib_groups::Vector{Int64} = d.loc_data.CB_CALIB_GROUPS
    spec_df = _filtered_model_spec(model_spec(d), cb_calib_groups)

    # Remove unguided scenarios as an option
    # Sample without unguided (i.e., values >= 1), then revert back to original model spec
    if !(spec_df[spec_df.fieldname .== :guided, :is_constant][1])
        _adjust_guided_lower_bound!(spec_df, 1)
        spec_df[!, :is_constant] .= spec_df[!, :lower_bound] .== spec_df[!, :upper_bound]
    end

    return sample(spec_df, n, sample_method)
end

"""
    sample_unguided(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Generate only unguided scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_unguided(
    d::Domain, n::Int64, sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    cb_calib_groups::Vector{Int64} = d.loc_data.CB_CALIB_GROUPS
    spec_df = _filtered_model_spec(model_spec(d), cb_calib_groups)

    # Fix guided factor to 0 (i.e., unguided scenarios only)
    guided_col = spec_df.fieldname .== :guided
    spec_df[guided_col, [:val, :lower_bound, :upper_bound, :dist_params, :is_constant]] .=
        [0 0 0 (0.0, 0.0) true]

    return sample(spec_df, n, sample_method)
end

"""
    deactivate_interventions(to_update::DataFrame)::Nothing
    deactivate_interventions!(dom::Domain)::Nothing

Deactivate all intervention factors (excluding `guided`) by settings these to 0.0.

# Arguments
- `to_update` : model specification to modify/update
- `dom` : Domain

# Returns
Nothing
"""
function deactivate_interventions(to_update::DataFrame)::Nothing
    intervs = component_params(to_update, Intervention)
    cols = Symbol[fn for fn in intervs.fieldname if fn != :guided]
    for c in cols
        _row = to_update.fieldname .== c
        _dparams =
            length(to_update[_row, :dist_params][1]) == 2 ? (0.0, 0.0) : (0.0, 0.0, 0.0)

        dval = _is_discrete_factor(to_update[_row, :ptype][1]) ? 0 : 0.0
        to_update[_row, [:val, :lower_bound, :upper_bound, :dist_params, :is_constant]] .=
            [dval 0.0 0.0 _dparams true]
    end

    return nothing
end
function deactivate_interventions!(dom::Domain)::Nothing
    iv_factors = component_params(model_spec(dom), Intervention).fieldname
    criteria_factors =
        component_params(
            model_spec(dom),
            [SeedCriteriaWeights, FogCriteriaWeights, decision.DepthThresholds]
        ).fieldname
    n_iv_factors = length(iv_factors)
    n_criteria_factors = length(criteria_factors)
    ADRIA.fix_factor!(dom; NamedTuple{Tuple(iv_factors)}(zeros(Float64, n_iv_factors))...)
    ADRIA.fix_factor!(
        dom; NamedTuple{Tuple(criteria_factors)}(zeros(Float64, n_criteria_factors))...
    )

    return nothing
end

"""
    fix_factor!(d::Domain, factor::Symbol)::Nothing
    fix_factor!(d::Domain, factor::Symbol, val::Real)::Nothing
    fix_factor!(d::Domain, factors::Vector{Symbol})::Nothing
    fix_factor!(d::Domain; factors...)::Nothing

Fix a factor so it gets ignored for the purpose of constructing samples.
If no value is provided, the default is used.

Note: Changes are permanent. To reset, either specify the original value(s)
      or reload the Domain.

# Examples
```julia
# Fix `guided` to default value
fix_factor!(dom, :guided)

# Fix `guided` to specified value
fix_factor!(dom, :guided, 3)

# Fix a set of factors to their default values
fix_factor!(dom, [:guided, :N_seed_TA])

# Fix specified factors to provided values
fix_factor!(dom; guided=3, N_seed_TA=1e6)
```
"""
function fix_factor!(d::Domain, factor::Symbol)::Nothing
    params = DataFrame(d.model)
    default_val = params[params.fieldname .== factor, :val][1]

    dist_params = params[params.fieldname .== factor, :dist_params][1]
    new_params = Tuple(fill(default_val, length(dist_params)))
    params[params.fieldname .== factor, :dist_params] .= [new_params]

    update!(d, params)
    return nothing
end
function fix_factor!(d::Domain, factor::Symbol, val::Real)::Nothing
    params = DataFrame(d.model)
    params[params.fieldname .== factor, :val] .= val

    dist_params = params[params.fieldname .== factor, :dist_params][1]
    new_dist_params = Tuple(fill(val, length(dist_params)))
    params[params.fieldname .== factor, :dist_params] .= [new_dist_params]

    update!(d, params)
    return nothing
end
function fix_factor!(d::Domain, factors::Vector{Symbol})::Nothing
    params = DataFrame(d.model)
    factor_rows = findall(in(factors), params.fieldname)

    # Get current values and dist_params lengths
    vals = params[factor_rows, :val]
    dist_lens = length.(params[factor_rows, :dist_params])

    # Create new dist_params tuples
    new_params = [Tuple(fill(v, len)) for (v, len) in zip(vals, dist_lens)]
    params[factor_rows, :dist_params] .= new_params

    update!(d, params)
    return nothing
end
function fix_factor!(d::Domain; factors...)::Nothing
    factor_names = keys(factors)
    factor_vals = collect(values(factors))

    params = DataFrame(d.model)

    target_order = vcat([findall(x -> x == y, params.fieldname) for y in factor_names]...)
    params[target_order, :val] .= factor_vals

    dist_params = params[target_order, :dist_params]
    new_dist_params = [
        Tuple(fill(v, length(d))) for (v, d) in zip(factor_vals, dist_params)
    ]
    params[target_order, :dist_params] .= new_dist_params

    update!(d, params)

    return nothing
end

"""
    _is_discrete_factor(dom::Domain, factor::Symbol)::Bool
    _is_discrete_factor(p_type::String)::Bool

Check `ptype` attribute whether the factor is a discrete variable type or not.
Returns `true` if discrete, `false` otherwise.

# Arguments
- `dom` : Domain
- `factor` : Name of model factor
- `ptype` : String representing variable type
"""
function _is_discrete_factor(dom::Domain, factor::Symbol)::Bool
    return _is_discrete_factor(get_attr(dom, factor, :ptype))
end
function _is_discrete_factor(p_type::String)::Bool
    return p_type ∈ DISCRETE_FACTOR_TYPES
end

"""
    get_bounds(dom::Domain, factor::Symbol)::NTuple{2,Float64}
    get_bounds(param::Param)::NTuple{2,Float64}

Get factor lower and upper bounds. If the factor has a triangular distribution, it returns
a 2-element tuple (without the peak value). Note that, for discrete factors, the actual
upper bound corresponds to the upper bound saved at the Domain's model_spec minus 1.0.

# Arguments
- `dom` : Domain
- `factor` : Factor name
- `param` : Parameter

# Returns
Minimum and maximum bounds associated with the parameter distribution.
"""
function get_bounds(dom::Domain, factor::Symbol)::NTuple{2,Float64}
    return get_attr(dom, factor, :dist_params)[1:2]
end
function get_bounds(param::Param)::NTuple{2,Float64}
    return param.dist_params[1:2]
end

"""
    get_default_bounds(dom::Domain, factor::Symbol)::NTuple{2,Float64}

"""
function get_default_bounds(dom::Domain, factor::Symbol)::NTuple{2,Float64}
    return get_attr(dom, factor, :default_dist_params)[1:2]
end

"""
    get_attr(dom::Domain, factor::Symbol, attr::Symbol)
"""
function get_attr(dom::Domain, factor::Symbol, attr::Symbol)
    ms = model_spec(dom)
    return ms[ms.fieldname .== factor, attr][1]
end

"""
    _update_decision_method!(dom, new_dist_params::Tuple)::Domain

Update the model spec with the tuple of MCDA methods to use.
"""
function _update_decision_method!(dom, new_dist_params::Tuple)::Domain
    new_method_names::Vector{String} = collect(new_dist_params)

    new_method_idxs = decision.decision_method_encoding.(new_method_names)

    ms = model_spec(dom)
    guided_row = findfirst(ms.fieldname .== :guided)
    @assert !isnothing(guided_row) "Guided variable not found in model spec."

    ms[guided_row, :dist_params] = Tuple(new_method_idxs)
    ms[guided_row, :val] = first(new_method_idxs)
    ms[guided_row, :is_constant] = length(new_method_idxs) == 1
    update!(dom, ms)

    return dom
end

"""
    set_factor_bounds!(dom::Domain, factor::Symbol, new_bounds::Tuple)::Nothing
    set_factor_bounds!(dom::Domain; factors...)::Nothing

Set new bound values for a given parameter. Sampled values for a parameter will lie
within the range `lower_bound ≤ s ≤ upper_bound`, for every sample value `s`.

Note: Changes are permanent. To reset, either specify the original value(s) or reload the
Domain.

# Arguments
- `dom` : Domain
- `factor` : Parameter whose bounds will be change to a new value
- `new_bounds` : Tuple bounds to be set as the new bounds of the respective factor. When
factor has a triangular distribution, `new_bounds` must be is a 3-element Tuple, with
`(new_min, new_max, new_peak)` values; when it has a uniform distribution, `new_bounds`
must be a 2-element Tuple, with `(new_lower, new_upper)` values.


# Examples
```julia
set_factor_bounds(dom, :wave_stress, (0.1, 0.2))
```
"""
function set_factor_bounds(dom::Domain, factor::Symbol, new_dist_params::Tuple)::Domain
    Base.@warn "set_factor_bounds is deprecated, use set_factor_bounds! instead"

    set_factor_bounds!(dom, factor, new_dist_params)
    return dom
end
function set_factor_bounds!(dom::Domain, factor::Symbol, new_dist_params::Tuple)::Domain
    if factor == :guided
        return _update_decision_method!(dom, new_dist_params)
    end

    # Set new base value
    old_val = get_attr(dom, factor, :val)
    new_val = mean(new_dist_params[1:2])
    (old_val isa Int) && (new_val = round(new_val))

    ms = model_spec(dom)
    ms[ms.fieldname .== factor, :dist_params] .= [new_dist_params]
    ms[ms.fieldname .== factor, :val] .= oftype(old_val, new_val)
    ms[!, :is_constant] .= (ms[!, :lower_bound] .== ms[!, :upper_bound])

    update!(dom, ms)
    return dom
end

function set_factor_bounds(dom::Domain; factors...)::Domain
    Base.@warn "set_factor_bounds is deprecated, use set_factor_bounds! instead"

    set_factor_bounds!(dom; factors...)
    return dom
end
function set_factor_bounds!(dom::Domain; factors...)::Domain
    factors = OrderedDict(factors)
    # Extract factor names and values
    factor_symbols = collect(keys(factors))
    new_params = collect(values(factors))

    # Handle categorical guided factor separately
    if :guided in factor_symbols
        factor_idx::Int64 = findfirst(factor_symbols .== :guided)
        dom = _update_decision_method!(dom, new_params[factor_idx])

        keep_idx = 1:length(factor_symbols) .!= factor_idx
        factor_symbols = factor_symbols[keep_idx]
        new_params = new_params[keep_idx]
    end

    ms = model_spec(dom)
    for (i, fn) in enumerate(factor_symbols)
        idx = findfirst(ms.fieldname .== fn)
        ms[idx, :dist_params] .= new_params[i]

        lb = new_params[i][1]
        ub = new_params[i][2]

        # Calculate new values preserving types
        old_val = ms[idx, :val]
        new_val = mean([lb, ub])
        ms[idx, :val] = old_val isa Int ? Int64(new_val) : new_val

        ms[idx, :lower_bound] = lb
        ms[idx, :upper_bound] = ub
        ms[idx, :is_constant] = lb == ub
    end

    update!(dom, ms)
    return dom
end

"""
    _validate_new_bounds(dom::Domain, factor::Symbol, new_dist_params::Tuple)::Nothing

Validate if new distribution bounds are within default bounds;
"""
function _validate_new_bounds(dom::Domain, factor::Symbol, new_bounds::Tuple)::Nothing
    err_msg = "Error setting new bound $new_bounds to $factor: bounds outside default range."
    default_lb, default_ub = get_default_bounds(dom, factor)
    (new_bounds[1] >= default_lb && new_bounds[2] <= default_ub) || error(err_msg)
    return nothing
end

_unzip(a) = map(x -> getfield.(a, x), fieldnames(eltype(a)))
_offdiag_iter(A) = collect(ι for ι in CartesianIndices(A) if ι[1] ≠ ι[2])

"""
    offdiag(A::AbstractArray)

Get off-diagonal values of matrix `A`.
"""
offdiag(A::AbstractArray) = A[_offdiag_iter(A)]

"""
    max_offdiag(A::AbstractArray)

Get the maximum off-diagonal values in matrix `A`.
"""
max_offdiag(A::AbstractArray) = maximum(offdiag(A))
function max_offdiag(df::DataFrame)
    return max_offdiag(cov(Matrix(df)))
end

"""
    max_maindiag(A::AbstractArray)

Get the maximum diagonal values in matrix `A`.
"""
max_maindiag(A::AbstractArray) = maximum(Matrix(I, size(A)...) .* A)
function max_maindiag(df::DataFrame)
    return max_maindiag(cov(Matrix(df)))
end
