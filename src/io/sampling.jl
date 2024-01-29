using Printf
using DataFrames, Distributions, LinearAlgebra
using ADRIA
using ADRIA: model_spec, component_params
using ADRIA.decision: mcda_normalize
import Surrogates: sample
import Surrogates.QuasiMonteCarlo as QMC
import Surrogates.QuasiMonteCarlo: SobolSample, OwenScramble

const DISCRETE_FACTOR_TYPES = ["ordered categorical", "unordered categorical", "discrete"]

"""
    _is_discrete_factor(p_type::String)::Bool
    _is_discrete_factor(dom, fieldname::Symbol)::Bool

Check ptype for discrete variable types. Returns true if discrete, false otherwise.

# Arguments
- `ptype` : String representing variable type
- `dom` : Domain
- `fieldname` : Name of model factor
"""
function _is_discrete_factor(p_type::String)::Bool
    return p_type ∈ DISCRETE_FACTOR_TYPES
end
function _is_discrete_factor(dom::Domain, fieldname::Symbol)::Bool
    model::Model = dom.model
    param_idx::Int64 = findfirst(x -> x == fieldname, model[:fieldname])
    ptype::String = model[:ptype][param_idx]
    return _is_discrete_factor(ptype)
end

"""
    _distribution_type(dom::Domain, fieldname::Symbol)
"""
function _distribution_type(dom::Domain, fieldname::Symbol)
    model::Model = dom.model
    param_idx::Int64 = findfirst(x -> x == fieldname, model[:fieldname])
    return model[:dist][param_idx]
end

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
    crit = component_params(spec, CriteriaWeights)
    interv = component_params(spec, Intervention)
    weights_seed_crit = criteria_params(crit, (:seed, :weight))
    weights_fog_crit = criteria_params(crit, (:fog, :weight))

    # If counterfactual, set all intervention options to 0.0
    df[df.guided .== -1.0, filter(x -> x ∉ [:guided, :heritability], interv.fieldname)] .=
        0.0

    # If unguided/counterfactual, set all preference criteria, except those related to depth, to 0.
    non_depth = filter(x -> x ∉ [:depth_min, :depth_offset], crit.fieldname)
    df[df.guided .== 0.0, non_depth] .= 0.0
    df[df.guided .== -1.0, non_depth] .= 0.0

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
    df[guided_fogged, weights_fog_crit.fieldname] .= mcda_normalize(
        df[guided_fogged, weights_fog_crit.fieldname]
    )
    # Normalize MCDA weights for seeding scenarios
    guided_seeded = .!(not_seeded) .& (df.guided .> 0)
    df[guided_seeded, weights_seed_crit.fieldname] .= mcda_normalize(
        df[guided_seeded, weights_seed_crit.fieldname]
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
function sample(dom::Domain, n::Int, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame
    n > 0 ? n : throw(DomainError(n, "`n` must be > 0"))
    return sample(model_spec(dom), n, sample_method)
end

"""
    sample(dom::Domain, n::Int, component::Type, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Create samples and rescale to distribution defined in the model spec.

Notes:
- assumes all parameters are independent.

# Arguments
- `dom` : Domain
- `n` : Number of samples to create
- `component` : Component type, e.g. CriteriaWeights
- `sample_method` : sample_method to use

# Returns
Scenario specification
"""
function sample(dom::Domain, n::Int, component::Type, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame
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
    sample_method=SobolSample(R=OwenScramble(base=2, pad=32))
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
    vary_dists = map(x -> x.dist(x.dist_params...), eachrow(vary_vars))

    # Create uniformly distributed samples for uncertain parameters
    samples = QMC.sample(n, zeros(n_vary_params), ones(n_vary_params), sample_method)

    # Scale uniform samples to indicated distributions using the inverse CDF method
    samples = Matrix(quantile.(vary_dists, samples)')

    # Combine varying and constant values (constant params use their indicated default vals)
    full_df = hcat(fill.(spec.val, n)...)
    full_df[:, spec.is_constant .== false] .= samples

    # Ensure unguided scenarios do not have superfluous factor combinations
    return adjust_samples(spec, DataFrame(full_df, spec.fieldname))
end

"""
    sample_site_selection(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Create guided samples of parameters relevant to site selection (EnvironmentalLayers, Intervention, CriteriaWeights).
All other parameters are set to their default values.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_site_selection(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame
    subset_spec = component_params(
        d.model, [EnvironmentalLayer, Intervention, CriteriaWeights]
    )

    # Only sample guided intervention scenarios
    _adjust_guided_lower_bound!(subset_spec, 1)

    # Create and fill scenario spec
    # Only Intervention, EnvironmentalLayer and CriteriaWeights factors are perturbed,
    # all other factors are fixed to their default values
    scens = repeat(param_table(d), n)
    select!(scens, Not(:RCP))  # remove RCP column added by param_table()
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
function sample_cf(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame
    spec_df = model_spec(d)

    # Unguided scenarios only
    guided_col = spec_df.fieldname .== :guided
    spec_df[guided_col, [:val, :lower_bound, :upper_bound, :dist_params, :is_constant]] .=
        [-1 -1 -1 (-1.0, -1.0) true]

    # Remove intervention scenarios as an option
    _deactivate_interventions(spec_df)

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
function sample_guided(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame
    spec_df = model_spec(d)

    # Remove unguided scenarios as an option
    # Sample without unguided (i.e., values >= 1), then revert back to original model spec
    _adjust_guided_lower_bound!(spec_df, 1)

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
function sample_unguided(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame
    spec_df = model_spec(d)

    # Fix guided factor to 0 (i.e., unguided scenarios only)
    guided_col = spec_df.fieldname .== :guided
    spec_df[guided_col, [:val, :lower_bound, :upper_bound, :dist_params, :is_constant]] .=
        [0 0 0 (0.0, 0.0) true]

    return sample(spec_df, n, sample_method)
end

"""
    _deactivate_interventions(to_update::DataFrame)::Nothing

Deactivate all intervention factors (excluding `guided`) by settings these to 0.0

# Arguments
- `to_update` : model specification to modify/update

# Returns
Scenario specification
"""
function _deactivate_interventions(to_update::DataFrame)::Nothing
    intervs = component_params(to_update, Intervention)
    cols = Symbol[fn for fn in intervs.fieldname if fn != :guided]
    for c in cols
        _row = to_update.fieldname .== c
        _dparams = length(to_update[_row, :dist_params][1]) == 2 ? (0.0, 0.0) : (0.0, 0.0, 0.0)

        dval = _is_discrete_factor(to_update[_row, :ptype][1]) ? 0 : 0.0
        to_update[_row, [:val, :lower_bound, :upper_bound, :dist_params, :is_constant]] .=
            [dval 0.0 0.0 _dparams true]
    end

    return nothing
end

"""
    fix_factor!(d::Domain, factor::Symbol)
    fix_factor!(d::Domain, factor::Symbol, val::Real)
    fix_factor!(d::Domain, factors...)

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
function fix_factor!(d::Domain; factors...)::Nothing
    for (factor, val) in factors
        try
            fix_factor!(d, factor, val)
        catch err
            if !(err isa MethodError)
                rethrow(err)
            end

            # Try setting value as an integer
            fix_factor!(d, factor, Int64(val))
        end
    end
    return nothing
end

"""
    get_bounds(dom::Domain, factor::Symbol)::Tuple
    get_bounds(param::Param)::Tuple

Get factor lower and upper bounds. If the factor has a triangular distribution, it returns
a 2-element tuple (without the peak value). Note that, for discrete factors, the actual
upper bound corresponds to the upper bound saved at the Domain's model_spec minus 1.0.

# Arguments
- `dom` : Domain
- `factor` : Name of the factor to get the bounds from
- `param` : Parameter

# Returns
Minimum and maximum bounds associated with the parameter distribution.
"""
function get_bounds(dom::Domain, factor::Symbol)::Tuple
    factor_filter::BitVector = collect(dom.model[:fieldname]) .== factor
    bounds::Tuple = dom.model[:dist_params][factor_filter][1]

    return (bounds[1], bounds[2])
end
function get_bounds(param::Param)::Tuple
    return param.dist_params[1:2]
end

"""
    lower_bound(param::Param)::Union{Int64, Float64}

Retrieve the lower bound of the parameter distribution
"""
function lower_bound(param::Param)::Union{Int64, Float64}
    return param.dist_params[1]
end

"""
    upper_bound(param::Param)::Union{Int64, Float64}

Retrieve the upper bound of the parameter distribution
"""
function upper_bound(param::Param)::Union{Int64, Float64}
    return param.dist_params[2]
end

"""
    get_default_bounds(dom::Domain, factor::Symbol)::Tuple

Get default distribution parameters for a factor.
Refer to `get_bounds` for more details of how the bounds work.

# Arguments
- `dom` : Domain
- `factor` : Name of the factor to get the bounds from
"""
function get_default_dist_params(dom::Domain, factor::Symbol)::Tuple
    model::Model = dom.model
    factor_filter::BitVector = collect(model[:fieldname]) .== factor
    default_params::Tuple = model[:default_dist_params][factor_filter][1]

    return default_params
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
set_factor_bounds!(dom, :wave_stress, (0.1, 0.2))
```
"""
function set_factor_bounds!(dom::Domain, factor::Symbol, new_bounds::Tuple)::Nothing
    _check_bounds_range(dom, factor, new_bounds)

    new_bounds, new_val = if _is_discrete_factor(dom, factor)
        _discrete_bounds(dom, factor, new_bounds)
    else
        _continuous_bounds(dom, factor, new_bounds)
    end

    params = model_spec(dom)
    params[params.fieldname .== factor, :dist_params] .= [new_bounds]
    params[params.fieldname .== factor, :val] .= new_val

    update!(dom, params)

    return nothing
end
function set_factor_bounds!(d::Domain; factors...)::Nothing
    for (factor, bounds) in factors
        set_factor_bounds!(d, factor, bounds)
    end

    return nothing
end

function _continuous_bounds(dom::Domain, factor::Symbol, new_bounds::Tuple)::Tuple
    new_lower, new_upper = new_bounds[1:2]
    new_val::Float64 = new_lower + 0.5 * (new_upper - new_lower)
    return (new_bounds, new_val)
end

function _discrete_bounds(dom::Domain, factor::Symbol, new_bounds::Tuple)::Tuple
    new_lower, new_upper = new_bounds[1:2]
    (new_lower % 1.0 == 0.0 && new_upper % 1.0 == 0.0) ||
        @warn "Upper and/or lower bounds for discrete variables should be integers."

    new_lower = round(new_lower)
    new_upper = min(ceil(new_upper), get_default_dist_params(dom, factor)[2])
    new_bounds = (new_lower, new_upper)

    new_val::Int64 = floor(new_lower + 0.5 * (new_upper - new_lower))

    return (new_bounds, new_val)
end

"""
    _check_bounds_range(dom::Domain, factor::Symbol, new_bounds::Tuple)::Nothing

Check new parameter bounds are within default parameter bounds
"""
function _check_bounds_range(dom::Domain, factor::Symbol, new_bounds::Tuple)::Nothing
    default_lower, default_upper = get_default_dist_params(dom, factor)
    new_lower, new_upper = new_bounds

    # Check if new bounds are within the default range
    out_of_bounds::Bool = (new_lower < default_lower) || (new_upper > default_upper)
    if out_of_bounds
        error(
            "Bounds should be within ($default_lower, $default_upper), received: ($new_lower, $new_upper).",
        )
    end

    if (_distribution_type(dom, factor) == "triang") && (length(new_bounds) !== 3)
        error("Triangular dist requires three parameters (minimum, maximum, peak).")
    elseif (_distribution_type(dom, factor) == "unif") && (length(new_bounds) !== 2)
        error("Uniform dist requires two parameters (minimum, maximum).")
    end

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
