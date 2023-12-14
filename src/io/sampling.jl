using Printf
using DataFrames, Distributions, LinearAlgebra
using ADRIA
using ADRIA: model_spec, _process_inputs!, component_params
using ADRIA.decision: mcda_normalize
import Surrogates: sample
import Surrogates.QuasiMonteCarlo: SobolSample


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
    # Map sampled values back to their discrete if necessary
    _process_inputs!(spec, df)

    crit = component_params(spec, CriteriaWeights)
    interv = component_params(spec, Intervention)
    weights_seed_crit = criteria_params(crit, (:seed, :weight))
    weights_fog_crit = criteria_params(crit, (:fog, :weight))

    # If counterfactual, set all intervention options to 0.0
    df[df.guided.==-1.0, filter(x -> x ∉ [:guided, :heritability], interv.fieldname)] .= 0.0

    # If unguided/counterfactual, set all preference criteria, except those related to depth, to 0.
    non_depth = filter(x -> x ∉ [:depth_min, :depth_offset], crit.fieldname)
    df[df.guided.==0.0, non_depth] .= 0.0
    df[df.guided.==-1.0, non_depth] .= 0.0

    # If unguided, set planning horizon to 0.
    df[df.guided.==0.0, :plan_horizon] .= 0.0

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
        df[guided_fogged, weights_fog_crit.fieldname],
    )
    # Normalize MCDA weights for seeding scenarios
    guided_seeded = .!(not_seeded) .& (df.guided .> 0)
    df[guided_seeded, weights_seed_crit.fieldname] .= mcda_normalize(
        df[guided_seeded, weights_seed_crit.fieldname],
    )

    # If use of distance threshold is off, set `dist_thresh` to 0.0
    df[df.use_dist.==0, :dist_thresh] .= 0.0

    if nrow(unique(df)) < nrow(df)
        perc = "$(@sprintf("%.3f", (1.0 - (nrow(unique(df)) / nrow(df))) * 100.0))%"
        @warn "Non-unique samples created: $perc of the samples are duplicates."
    end

    return df
end

"""
    sample(dom::Domain, n::Int, sampler=SobolSample())::DataFrame

Create samples and rescale to distribution defined in the model spec.

Notes:
- assumes all parameters are independent.

# Arguments
- `dom` : Domain
- `n` : Int
- `sampler` : type of sampler to use.

# Returns
Scenario specification
"""
function sample(dom::Domain, n::Int, sampler=SobolSample())::DataFrame
    n > 0 ? n : throw(DomainError(n, "`n` must be > 0"))
    return sample(model_spec(dom), n, sampler)
end

"""
    sample(dom::Domain, n::Int, component::Type)::DataFrame

Create samples and rescale to distribution defined in the model spec.

Notes:
- assumes all parameters are independent.

# Arguments
- `dom` : Domain
- `n` : Int
- `component` : Type, e.g. CriteriaWeights
- `sampler` : type of sampler to use.

# Returns
Scenario specification
"""
function sample(dom::Domain, n::Int, component::Type, sampler=SobolSample())::DataFrame
    n > 0 ? n : throw(DomainError(n, "`n` must be > 0"))

    spec = component_params(dom.model, component)
    return sample(spec, n, sampler)
end

"""
    sample(spec::DataFrame, n::Int, sampler=SobolSample())::DataFrame

Create samples and rescale to distribution defined in the model spec.

# Arguments
- `spec` : DataFrame containing model parameter specifications.
- `n` : number of samples to generate.
- `sampler` : type of sampler to use.

# Returns
Scenario specification
"""
function sample(spec::DataFrame, n::Int64, sampler=SobolSample(); supported_dists=Dict(
    "triang" => TriangularDist,
    "norm" => TruncatedNormal,
    "unif" => Uniform
))::DataFrame

    if Symbol(sampler) == Symbol("QuasiMonteCarlo.SobolSample()")
        ispow2(n) ? n : throw(DomainError(n, "`n` must be a power of 2 when using the Sobol' sampler"))
    end

    # Select non-constant params
    vary_vars = spec[spec.is_constant.==false, ["dists", "bounds"]]

    # Update range
    triang_params = vary_vars[vary_vars.dists.=="triang", "bounds"]
    vary_vars[vary_vars.dists.=="triang", "bounds"] .= map(x -> (x[1], x[2], (x[2] - x[1]) * x[3] + x[1]), triang_params)
    vary_dists = map((x) -> supported_dists[x.dists](x.bounds...), eachrow(vary_vars))

    # Create sample for uncertain parameters
    n_vary_params = size(vary_vars, 1)
    n_vary_params > 0 ? n_vary_params : throw(DomainError(n_vary_params, "Number of parameters to perturb must be > 0"))
    samples = sample(n, zeros(n_vary_params), ones(n_vary_params), sampler)

    # Convert vector of tuples to matrix
    samples = permutedims(hcat([collect(s) for s in samples]...))

    # Scale values to indicated distributions
    samples .= permutedims(hcat(map(ix -> quantile.(vary_dists[ix], samples[:, ix]), 1:size(samples, 2))...))'

    # Combine varying and constant values (constant params use their indicated default vals)
    full_df = hcat(fill.(spec.val, n)...)
    full_df[:, spec.is_constant.==false] .= samples

    # Adjust samples for discrete values using flooring trick
    # Ensure unguided scenarios do not have superfluous parameter values
    return adjust_samples(spec, DataFrame(full_df, spec.fieldname))
end

"""
    sample_site_selection(d::Domain, n::Int64, sampler=SobolSample())::DataFrame

Create guided samples of parameters relevant to site selection (EnvironmentalLayers, Intervention, CriteriaWeights).
All other parameters are set to their default values.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sampler` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_site_selection(d::Domain, n::Int64, sampler=SobolSample())::DataFrame
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
    scens[:, subset_spec.fieldname] .= sample(subset_spec, n, sampler)

    return scens
end

"""
    sample_cf(d::Domain, n::Int64, sampler=SobolSample())::DataFrame

Generate only counterfactual scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sampler` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_cf(d::Domain, n::Int64, sampler=SobolSample())::DataFrame
    spec_df = model_spec(d)

    # Unguided scenarios only
    guided_col = spec_df.fieldname .== :guided
    spec_df[guided_col, [:val, :lower_bound, :upper_bound, :bounds, :is_constant]] .= [-1 -1 -1 (-1.0, -1.0) true]

    # Remove intervention scenarios as an option
    _deactivate_interventions(spec_df)

    return sample(spec_df, n, sampler)
end

"""
    _adjust_guided_lower_bound!(guided_spec::DataFrame, lower::Int64)::DataFrame

Adjust lower bound of guided parameter spec to alter sampling range.
"""
function _adjust_guided_lower_bound!(spec_df::DataFrame, lower::Int64)::DataFrame
    guided_col = spec_df.fieldname .== :guided
    g_upper = Float64(spec_df[guided_col, :upper_bound][1])

    # Update entries, standardizing values for bounds as floats
    spec_df[guided_col, [:val, :lower_bound, :bounds]] .= [lower Float64(lower) (Float64(lower), g_upper)]
    return spec_df
end

"""
    sample_guided(d::Domain, n::Int64, sampler=SobolSample())::DataFrame

Generate only guided scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sampler` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_guided(d::Domain, n::Int64, sampler=SobolSample())::DataFrame
    spec_df = model_spec(d)

    # Remove unguided scenarios as an option
    # Sample without unguided (i.e., values >= 1), then revert back to original model spec
    _adjust_guided_lower_bound!(spec_df, 1)

    return sample(spec_df, n, sampler)
end

"""
    sample_unguided(d::Domain, n::Int64, sampler=SobolSample())::DataFrame

Generate only unguided scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sampler` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_unguided(d::Domain, n::Int64, sampler=SobolSample())::DataFrame
    spec_df = model_spec(d)

    # Fix guided factor to 0 (i.e., unguided scenarios only)
    guided_col = spec_df.fieldname .== :guided
    spec_df[guided_col, [:val, :lower_bound, :upper_bound, :bounds, :is_constant]] .= [0 0 0 (0.0, 0.0) true]

    return sample(spec_df, n, sampler)
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
        _bnds = length(to_update[_row, :bounds][1]) == 2 ? (0.0, 0.0) : (0.0, 0.0, 0.0)

        dval = _check_discrete(to_update[_row, :ptype][1]) ? 0 : 0.0
        to_update[_row, [:val, :lower_bound, :upper_bound, :bounds, :is_constant]] .= [dval 0.0 0.0 _bnds true]
    end

    return nothing
end


"""
    fix_factor(d::Domain, factor::Symbol)
    fix_factor(d::Domain, factor::Symbol, val::Real)
    fix_factor(d::Domain, factors...)

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
    default_val = params[params.fieldname.==factor, :val][1]

    bnds = params[params.fieldname.==factor, :bounds][1]
    new_bnds = Tuple(fill(default_val, length(bnds)))
    params[params.fieldname.==factor, :bounds] .= [new_bnds]

    update!(d, params)
end
function fix_factor!(d::Domain, factor::Symbol, val::Real)::Nothing
    params = DataFrame(d.model)
    params[params.fieldname.==factor, :val] .= val

    bnds = params[params.fieldname.==factor, :bounds][1]
    new_bnds = Tuple(fill(val, length(bnds)))
    params[params.fieldname.==factor, :bounds] .= [new_bnds]

    update!(d, params)
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
end

"""
    get_bounds(dom::Domain, factor::Symbol)::Tuple

Get factor bounds. If the factor is continuous, it will have a triangular distribution,
then bounds will be a tuple with three elements (minimum, maximum, peak). If the factor is
discrete, it will have a uniform distribution, then bounds will be a tuple with two elements
(lower_bound, upper_bound). Note that, for discrete factors, the actual upper bound
corresponds to the upper bound saved at the Domain's model_spec minus 1.0 because of how
sampling works.

# Arguments
- `dom` : Domain
- `factor` : Name of the factor to get the bounds from
"""
function get_bounds(dom::Domain, factor::Symbol)::Tuple
    model::Model = dom.model
    factor_filter::BitVector = collect(model[:fieldname]) .== factor

    bounds = if _check_discrete(dom, factor)
        lower, upper = model[:bounds][factor_filter][1]
        (lower, upper - 1.0)
    else
        model[:bounds][factor_filter][1]
    end

    return bounds
end

"""
    get_default_bounds(dom::Domain, factor::Symbol)::Tuple

Get factor default_bounds. Refer to `get_bounds` for more details of how the bounds work.

# Arguments
- `dom` : Domain
- `factor` : Name of the factor to get the bounds from
"""
function get_default_bounds(dom::Domain, factor::Symbol)::Tuple
    model::Model = dom.model
    factor_filter::BitVector = collect(model[:fieldname]) .== factor

    default_bounds = if _check_discrete(dom, factor)
        default_lower, default_upper = model[:default_bounds][factor_filter][1]
        (default_lower, default_upper - 1.0)
    else
        model[:default_bounds][factor_filter][1]
    end

    return default_bounds
end

"""
    set_factor_bounds!(dom::Domain, factor::Symbol, new_bounds::Tuple)::Nothing
    set_factor_bounds!(dom::Domain; factors...)::Nothing

Set new lower and upper bounds for a parameter, `new_bounds = (lower_bound, upper_bound)`.
All sampled values `s_vals` for this parameter will lie within the range
`lower_bound ≤ s_vals ≤ upper_bound`. When used to set new bounds for a discrete parameter,
the upper bound will be set to `upper_bound + 1` to guarantee that the sample values will
always be less or equal `upper_bound`.

Note: Changes are permanent. To reset, either specify the original value(s) or reload the
Domain.

# Arguments
- `dom` : Domain
- `factor` : Parameter whose bounds will be change to a new value
- `new_bounds` : Tuple of lower and upper bounds to be set as the new bounds of the
respective parameter


# Examples
```julia
set_factor_bounds!(dom, :wave_stress, (0.1,0.2))
```
"""
function set_factor_bounds!(dom::Domain, factor::Symbol, new_bounds::Tuple)::Nothing
    params = DataFrame(dom.model)

    # Get upper and lower bounds for default and new distributions
    default_lower, default_upper = params[params.fieldname .== factor, :default_bounds][1]
    new_lower, new_upper = new_bounds[1], new_bounds[2]

    _check_bounds_range(new_bounds, (default_lower, default_upper), factor, dom)

    if (params[params.fieldname .== factor, :dists][1] == "triang") &&
        (length(new_bounds) !== 3)
        error("Triangular dist requires three parameters (minimum, maximum, peak).")
    elseif (params[params.fieldname .== factor, :dists][1] == "unif") &&
        (length(new_bounds) !== 2)
        error("Uniform dist requires two parameters (minimum, maximum).")
    end

    new_val = new_lower + 0.5 * (new_upper - new_lower)

    if _check_discrete(dom, factor)
        if new_lower % 1.0 != 0.0 || new_upper % 1.0 != 0.0
            @warn "Upper and/or lower bounds for discrete variables should be integer numbers."
        end

        new_val = floor(Int64, new_val)
        new_lower = round(new_lower)
        # upper bound should not be greater than default_upper
        new_upper = min(round(new_upper) + 1.0, default_upper)
        new_bounds = (new_lower, new_upper)
    end

    params[params.fieldname .== factor, :bounds] .= [new_bounds]
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

_unzip(a) = map(x -> getfield.(a, x), fieldnames(eltype(a)))

"""
    function _check_bounds_range(new_bounds::Tuple, factor::Symbol, dom::Domain)

    Check new parameter bounds are within default parameter bounds
"""
function _check_bounds_range(
    new_bounds::Tuple, default_bounds::Tuple, factor::Symbol, dom::Domain
)::Nothing
    out_of_bounds::Bool = false
    new_lower, new_upper = new_bounds
    default_lower, default_upper = default_bounds

    if _check_discrete(dom, factor)
        # If factor is discrete, sampled values must be within default_lower and (default_upper - 1)
        out_of_bounds = (new_lower < default_lower) || (new_upper > (default_upper - 1))
    else
        out_of_bounds = (new_lower < default_lower) || (new_upper > default_upper)
    end

    if out_of_bounds
        error(
            "Bounds should be within ($default_lower, $default_upper), received: ($new_lower, $new_upper).",
        )
    end

    return nothing
end

"""
Check specified bounds for validity.

Raises error if lower bound values are greater than upper bounds.

# Arguments
- `lower` : lower bounds
- `upper` : upper bound values
"""
function _check_bounds(lower, upper)
    if any(lower .> upper)
        error("Bounds are not legal (upper bound must be greater than lower bound)")
    end
end

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
