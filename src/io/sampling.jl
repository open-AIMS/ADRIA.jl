using Printf
using DataFrames, Distributions, LinearAlgebra
using ADRIA
using ADRIA: model_spec, _process_inputs!, component_params
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

    crit = component_params(spec, Criteria)
    interv = component_params(spec, Intervention)

    # If counterfactual, set all intervention options to 0.0
    df[df.guided.==-1.0, filter(x -> x ∉ [:guided, :n_adapt], interv.fieldname)] .= 0.0

    # If unguided/counterfactual, set all preference criteria, except those related to depth, to 0.
    non_depth = filter(x -> x ∉ [:depth_min, :depth_offset], crit.fieldname)
    df[df.guided.==0.0, non_depth] .= 0.0
    df[df.guided.==-1.0, non_depth] .= 0.0

    # If unguided, set planning horizon to 0.
    df[df.guided.==0.0, :plan_horizon] .= 0.0

    # If no seeding is to occur, set related variables to 0
    not_seeded = (df.seed_TA .== 0) .& (df.seed_CA .== 0)
    df[not_seeded, contains.(names(df), "seed_")] .= 0.0
    df[not_seeded, :a_adapt] .= 0.0

    # Same for fogging/shading
    not_fogged = (df.fogging .== 0) .& (df.SRM .== 0)
    df[not_fogged, contains.(names(df), "shade_")] .= 0.0

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
- `component` : Type, e.g. Criteria
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
function sample(spec::DataFrame, n::Int, sampler=SobolSample(); supported_dists=Dict(
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
    sample_site_selection(d::Domain, n::Int, sampler=SobolSample())::DataFrame

Create guided samples of parameters relevant to site selection (EnvironmentalLayers, Intervention, Criteria).
All other parameters are set to their default values.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sampler` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_site_selection(d::Domain, n::Int, sampler=SobolSample())::DataFrame
    subset_spec = component_params(d.model, [EnvironmentalLayer, Intervention, Criteria])

    # Only sample guided intervention scenarios
    _adjust_guided_lower_bound!(subset_spec, 1)

    # Create and fill scenario spec
    # Only Intervention, EnvironmentalLayer and Criteria factors are perturbed,
    # all other factors are fixed to their default values
    scens = repeat(param_table(d), n)
    scens[:, subset_spec.fieldname] .= sample(subset_spec, n, sampler)

    return scens
end

"""
    sample_cf(d::Domain, n::Int, sampler=SobolSample())::DataFrame
    
Generate only counterfactual scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sampler` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_cf(d::Domain, n::Int, sampler=SobolSample())::DataFrame
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
    sample_guided(d::Domain, n::Int, sampler=SobolSample())::DataFrame
    
Generate only guided scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sampler` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_guided(d::Domain, n::Int, sampler=SobolSample())::DataFrame
    spec_df = model_spec(d)

    # Remove unguided scenarios as an option
    # Sample without unguided (i.e., values >= 1), then revert back to original model spec
    _adjust_guided_lower_bound!(spec_df, 1)

    return sample(spec_df, n, sampler)
end

"""
    sample_unguided(d::Domain, n::Int, sampler=SobolSample())::DataFrame
    
Generate only unguided scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sampler` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_unguided(d::Domain, n::Int, sampler=SobolSample())::DataFrame
    spec_df = model_spec(d)

    # Fix guided factor to 0 (i.e., unguided scenarios only)
    guided_col = spec_df.fieldname .== :guided
    spec_df[guided_col, [:val, :lower_bound, :upper_bound, :bounds, :is_constant]] .= [0 0 0 (0.0, 0.0) true]

    return sample(spec_df, n, sampler)
end

"""
    _deactivate_interventions(to_update::DataFrame)::Nothing

Deactivate all intervention factors (excluding `guided` and `n_adapt`) by settings these to 0.0

# Arguments
- `to_update` : model specification to modify/update

# Returns
Scenario specification
"""
function _deactivate_interventions(to_update::DataFrame)::Nothing
    intervs = component_params(to_update, Intervention)
    cols = Symbol[fn for fn in intervs.fieldname if fn != :n_adapt && fn != :guided]
    for c in cols
        _row = to_update.fieldname .== c
        _bnds = length(to_update[_row, :bounds][1]) == 2 ? (0.0, 0.0) : (0.0, 0.0, 0.0)

        dval = to_update[_row, :ptype][1] == "integer" ? 0 : 0.0
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
fix_factor!(dom; guided=3, seed_TA=1e6)
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
function fix_factor!(d::Domain; factors...)
    for (factor, val) in factors
        fix_factor!(d, factor, val)
    end
end


_unzip(a) = map(x -> getfield.(a, x), fieldnames(eltype(a)))


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
