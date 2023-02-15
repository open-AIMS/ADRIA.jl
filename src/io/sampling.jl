using Printf
using DataFrames, Distributions, LinearAlgebra
using ADRIA
import ADRIA: model_spec, process_inputs!, component_params
import Surrogates: sample
import Surrogates.QuasiMonteCarlo: SobolSample


"""
    adjust_samples(d::Domain, df::DataFrame)::DataFrame
    adjust_samples!(d::Domain, spec::DataFrame, df::DataFrame)::DataFrame

Adjust given samples to ensure parameter value combinations for unguided
scenarios are plausible.
"""
function adjust_samples(d::Domain, df::DataFrame)::DataFrame
    return adjust_samples(d, model_spec(d), df)
end
function adjust_samples(d::Domain, spec::DataFrame, df::DataFrame)::DataFrame
    process_inputs!(d, df)
    crit = component_params(spec, Criteria)
    interv = component_params(spec, Intervention)

    # If counterfactual, set all intervention options to 0.0
    df[df.guided.==-1.0, filter(x -> x ∉ [:guided, :n_adapt], interv.fieldname)] .= 0.0

    # If unguided/counterfactual, set all preference criteria, except those related to depth, to 0.
    non_depth = filter(x -> x ∉ [:depth_min, :depth_offset], crit.fieldname)
    df[df.guided.==0.0, non_depth] .= 0.0
    df[df.guided.==-1.0, non_depth] .= 0.0

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

function adjust_cf_samples(d::Domain, spec::DataFrame, df::DataFrame)::DataFrame
    # Get intervention columns that aren't "n_adapt"
    # and make them == 0
    intervs = component_params(spec, Intervention)
    cols = collect(skipmissing([fn != :n_adapt ? fn : missing for fn in intervs.fieldname]))
    df[:, cols] .= 0.0
    df[:, [:guided]] .= -1.0  # Mark counterfactuals as -1

    return adjust_samples(d, df)
end
function adjust_cf_samples(d::Domain, df::DataFrame)::DataFrame
    return adjust_cf_samples(d, model_spec(d), df)
end


"""
    sample(dom::Domain, n::Int)::DataFrame

Create samples and rescale to distribution defined in the model spec.

Notes:
- assumes all parameters are independent.

# Arguments
- `dom` : Domain
- `n` : Int
"""
function sample(dom::Domain, n::Int, sampler=SobolSample())::DataFrame
    n > 0 ? n : throw(DomainError(n, "`n` must be > 0"))

    spec = model_spec(dom)
    df = sample(spec, n, sampler)

    # Adjust samples for discrete values using flooring trick
    # Ensure unguided scenarios do not have superfluous parameter values
    return adjust_samples(dom, spec, df)
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
"""
function sample(dom::Domain, n::Int, component::Type, sampler=SobolSample())::DataFrame
    n > 0 ? n : throw(DomainError(n, "`n` must be > 0"))

    spec = component_params(dom.model, component)
    df = sample(spec, n, sampler)
    # Adjust samples for discrete values using flooring trick
    # Ensure unguided scenarios do not have superfluous parameter values
    process_inputs!(spec, df)

    return df
end

"""
    sample(spec::DataFrame, n::Int, sampler=SobolSample())::DataFrame

Create samples and rescale to distribution defined in the model spec.

# Arguments
- `spec` : DataFrame containing model parameter specifications.
- `n` : number of samples to generate.
- `sampler` : type of sampler to use.
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
    vary_vars = spec[spec.is_constant.==false, ["dists", "full_bounds"]]

    # Update range
    triang_params = vary_vars[vary_vars.dists.=="triang", "full_bounds"]
    vary_vars[vary_vars.dists.=="triang", "full_bounds"] .= map(x -> (x[1], x[2], (x[2] - x[1]) * x[3] + x[1]), triang_params)
    vary_dists = map((x) -> supported_dists[x.dists](x.full_bounds...), eachrow(vary_vars))

    # Create sample for uncertain parameters
    n_vary_params = size(vary_vars, 1)
    n_vary_params > 0 ? n_vary_params : throw(DomainError(n_vary_params, "Number of parameters to perturb must be > 0"))
    samples = sample(n, zeros(n_vary_params), ones(n_vary_params), sampler)

    # Convert vector of tuples to matrix
    samples = permutedims(hcat([collect(s) for s in samples]...))

    # Scale values to indicated distributions
    samples .= permutedims(hcat(map(ix -> quantile.(vary_dists[ix], samples[:, ix]), 1:size(samples, 2))...))'

    # Combine varying and constant values
    full_df = zeros(n, size(spec, 1))
    full_df[:, findall(spec.is_constant .== false)] .= samples

    df = DataFrame(full_df, spec.fieldname)

    # Adjust samples for discrete values using flooring trick
    # Ensure unguided scenarios do not have superfluous parameter values
    return df
end

"""
    sample_site_selection(d::Domain, n::Int, sampler=SobolSample())::DataFrame

Create samples of only site selection parameters and rescale to distribution defined in the model spec.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sampler` : type of sampler to use.
"""
function sample_site_selection(d::Domain, n::Int, sampler=SobolSample())::DataFrame

    crit_spec = component_params(d.model, Criteria)
    env_spec = component_params(d.model, EnvironmentalLayer)

    int_spec = component_params(d.model, Intervention)
    insertcols!(int_spec, :val, :bounds => copy([int_spec[:, :full_bounds]...]))
    guided_spec = _adjust_guided_lower_bound(int_spec[int_spec[:, :fieldname].==:guided, :], 1)
    select!(guided_spec, Not(:bounds))

    sample_df = vcat(vcat(env_spec, guided_spec), crit_spec)
    site_selection_sample = sample(sample_df, n, sampler)

    process_inputs!(sample_df, site_selection_sample)
    return site_selection_sample
end

"""
    sample_cf(d::Domain, n::Int, sampler=SobolSample())::DataFrame
    
Generate only counterfactual scenarios using any sampler from QuasiMonteCarlo.jl
"""
function sample_cf(d::Domain, n::Int, sampler=SobolSample())::DataFrame
    df::DataFrame = sample(d, n, sampler)
    return adjust_cf_samples(d, df)
end

"""
    adjust_guided_bounds(guided_spec::DataFrame, lower::Int64)::DataFrame
    
Adjust lower bound of guided parameter spec to alter sampling range.
"""
function _adjust_guided_lower_bound(spec_df::DataFrame, lower::Int64)::DataFrame
    guided_col = spec_df.fieldname .== :guided
    g_upper = spec_df[guided_col, :upper_bound]
    spec_df[guided_col, [:val, :lower_bound, :bounds, :full_bounds]] .= reshape([lower, lower, (lower, g_upper), (lower, g_upper)], (1, 4))

    return spec_df
end

"""
    sample_guided(d::Domain, n::Int, sampler=SobolSample())::DataFrame
    
Generate only guided scenarios using any sampler from QuasiMonteCarlo.jl
"""
function sample_guided(d::Domain, n::Int, sampler=SobolSample())::DataFrame
    spec_df = model_spec(d)

    # Remove unguided scenarios as an option
    mod_df = copy(spec_df)

    insertcols!(mod_df, :val, :bounds => copy([d.model[:bounds]...]))
    mod_df = _adjust_guided_lower_bound(mod_df, 1)

    # Sample without unguided, then revert back to original model spec
    ADRIA.update!(d.model, mod_df)
    samples = adjust_samples(d, sample(d, n, sampler))

    # Note: updating with spec_df does not work.
    mod_df = _adjust_guided_lower_bound(mod_df, 0)

    ADRIA.update!(d.model, mod_df)

    return samples
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
