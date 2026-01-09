import OrderedCollections: OrderedDict
using
    DataFrames,
    Distributions,
    LinearAlgebra,
    Printf

import Distributions: sample
import QuasiMonteCarlo as QMC
import QuasiMonteCarlo: SobolSample, OwenScramble

using ADRIA: model_spec, component_params
using .decision: mcda_normalize

include("sampling_interface.jl")
include("sampling_cf.jl")
include("sampling_unguided.jl")
include("sampling_guided.jl")
include("sampling_interventions.jl")
include("sampling_balanced.jl")

const DISCRETE_FACTOR_TYPES = [
    "ordered categorical", "unordered categorical", "ordered discrete"
]

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
