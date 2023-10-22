module sensitivity

using Logging

using
    AxisKeys,
    NamedDims,
    StaticArrays
using DataFrames

using
    Bootstrap,
    Distributions,
    HypothesisTests,
    StatsBase
using HypothesisTests: ApproximateKSTest

using FLoops

using ADRIA: ResultSet
using ADRIA.analysis: col_normalize


"""
    ks_statistic(ks)

Calculate the Kolmogorov-Smirnov test statistic.
"""
function ks_statistic(ks::ApproximateKSTest)::Float64
    n::Float64 = (ks.n_x * ks.n_y) / (ks.n_x + ks.n_y)

    return sqrt(n) * ks.δ
end


"""
    pawn(rs::ResultSet, y::Union{NamedDimsArray,AbstractVector{<:Real}}; S::Int64=10)::NamedDimsArray
    pawn(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}, factor_names::Vector{String}; S::Int64=10)::NamedDimsArray
    pawn(X::DataFrame, y::AbstractVector{<:Real}; S::Int64=10)::NamedDimsArray
    pawn(X::NamedDimsArray, y::Union{NamedDimsArray,AbstractVector{<:Real}}; S::Int64=10)::NamedDimsArray
    pawn(X::Union{DataFrame,AbstractMatrix{<:Real}}, y::AbstractMatrix{<:Real}; S::Int64=10)::NamedDimsArray

Calculates the PAWN sensitivity index.

The PAWN method (by Pianosi and Wagener) is a moment-independent approach to Global
Sensitivity Analysis. Outputs are characterized by their Cumulative Distribution Function
(CDF), quantifying the variation in the output distribution after conditioning an input over
"slices" (\$S\$) - the conditioning intervals. If both distributions coincide at all slices
(i.e., the distributions are similar or identical), then the factor is deemed
non-influential.

This implementation applies the Kolmogorov-Smirnov test as the distance measure and returns
summary statistics (min, mean, median, max, std, and cv) over the slices.

# Arguments
- `rs` : ResultSet
- `X` : Model inputs
- `y` : Model outputs
- `factor_names` : Names of each factor represented by columns in `X`
- `S` : Number of slides (default: 10)

# Returns
NamedDimsArray, of min, mean, median, max, std, and cv summary statistics.

# Examples
```julia
dom = ADRIA.load_domain("example_domain")
scens = ADRIA.sample(dom, 128)
rs = ADRIA.run_scenarios(dom, scens, "45")

# Get mean coral cover over time and locations
μ_tac = mean(ADRIA.metrics.scenario_total_cover(rs), dims=:timesteps)

ADRIA.sensitivity.pawn(rs, μ_tac)
```

# References
1. Pianosi, F., Wagener, T., 2018.
   Distribution-based sensitivity analysis from a generic input-output sample.
   Environmental Modelling & Software 108, 197-207.
   https://doi.org/10.1016/j.envsoft.2018.07.019

2. Baroni, G., Francke, T., 2020.
   GSA-cvd
   Combining variance- and distribution-based global sensitivity analysis
   https://github.com/baronig/GSA-cvd

3. Puy, A., Lo Piano, S., & Saltelli, A. 2020.
   A sensitivity analysis of the PAWN sensitivity index.
   Environmental Modelling & Software, 127, 104679.
   https://doi.org/10.1016/j.envsoft.2020.104679

4. https://github.com/SAFEtoolbox/Miscellaneous/blob/main/Review_of_Puy_2020.pdf

# Extended help
Pianosi and Wagener have made public their review responding to a critique of their method
by Puy et al., (2020). A key criticism by Puy et al. was that the PAWN method is sensitive
to its tuning parameters and thus may produce biased results. The tuning parameters referred
to are the number of samples (\$N\$) and the number of conditioning points - \$n\$ in Puy et
al., but denoted as \$S\$ here.

Puy et al., found that the ratio of \$N\$ (number of samples) to \$S\$ has to be
sufficiently high (\$N/S > 80\$) to avoid biased results. Pianosi and Wagener point out this
requirement is not particularly difficult to meet. Using the recommended value
(\$S := 10\$), a sample of 1024 runs (small for purposes of Global Sensitivity Analysis)
meets this requirement (\$1024/10 = 102.4\$). Additionally, lower values of \$N/S\$ is more
an indication of faulty experimental design moreso than any deficiency of the PAWN method.
"""
function pawn(
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    factor_names::Vector{String};
    S::Int64=10
)::NamedDimsArray
    N, D = size(X)
    step = 1 / S
    seq = 0.0:step:1.0

    # Preallocate result structures
    X_q = @MVector zeros(S + 1)
    pawn_t = @MArray zeros(S, D)
    results = @MArray zeros(D, 6)

    # Hide warnings from HypothesisTests
    with_logger(NullLogger()) do
        @floop for d_i in 1:D
            X_di = @view(X[:, d_i])
            X_q .= quantile(X_di, seq)

            Y_sel = @view(y[X_q[1].<=X_di.<=X_q[2]])
            if length(Y_sel) > 0
                pawn_t[1, d_i] = ks_statistic(ApproximateTwoSampleKSTest(Y_sel, y))
            end

            for s in 2:S
                Y_sel = @view(y[X_q[s].<X_di.<=X_q[s+1]])
                if length(Y_sel) == 0
                    continue  # no available samples
                end

                pawn_t[s, d_i] = ks_statistic(ApproximateTwoSampleKSTest(Y_sel, y))
            end

            p_ind = @view(pawn_t[:, d_i])
            p_mean = mean(p_ind)
            p_sdv = std(p_ind)
            p_cv = p_sdv ./ p_mean
            results[d_i, :] .= (
                minimum(p_ind),
                p_mean,
                median(p_ind),
                maximum(p_ind),
                p_sdv,
                p_cv
            )
        end
    end

    replace!(results, NaN => 0.0, Inf => 0.0)

    col_names = [:min, :mean, :median, :max, :std, :cv]
    return NamedDimsArray(results; factors=Symbol.(factor_names), Si=col_names)
end
function pawn(X::DataFrame, y::AbstractVector{<:Real}; S::Int64=10)::NamedDimsArray
    return pawn(Matrix(X), y, names(X); S=S)
end
function pawn(
    X::NamedDimsArray,
    y::Union{NamedDimsArray,AbstractVector{<:Real}};
    S::Int64=10
)::NamedDimsArray
    return pawn(X, y, axiskeys(X, 2); S=S)
end
function pawn(
    X::Union{DataFrame,NamedDimsArray},
    y::AbstractMatrix{<:Real};
    S::Int64=10
)::NamedDimsArray
    N, D = size(y)
    if N > 1 && D > 1
        msg::String = string(
            "The current implementation of PAWN can only assess a single quantity",
            " of interest at a time."
        )
        throw(ArgumentError(msg))
    end

    # The wrapped call to `vec()` handles cases where matrix-like data type is passed in
    # (N x 1 or 1 x D) and so ensures a vector is passed along
    return pawn(X, vec(y); S=S)
end
function pawn(
    rs::ResultSet,
    y::Union{NamedDimsArray,AbstractVector{<:Real}};
    S::Int64=10
)::NamedDimsArray
    return pawn(rs.inputs, y; S=S)
end


"""
    tsa(X::DataFrame, y::AbstractMatrix)::NamedDimsArray

Perform Temporal (or time-varying) Sensitivity Analysis using the PAWN sensitivity index.

The sensitivity index value for time \$t\$ is inclusive of all time steps prior to \$t\$.
Alternate approaches use a moving window, or only data for time \$t\$.

# Examples
```julia
rs = ADRIA.load_results("a ResultSet of interest")

# Get scenario outcomes over time (shape: `time ⋅ scenarios`)
y_tac = ADRIA.metrics.scenario_total_cover(rs)

# Calculate sensitivity of outcome to factors for each time step
ADRIA.sensitivity.tsa(rs.inputs, y_tac)
```

# Arguments
- `X` : Scenario specification
- `y` : scenario outcomes over time

# Returns
NamedDimsArray, of shape \$D\$ ⋅ 6 ⋅ \$T\$, where
- \$D\$ is the number of dimensions/factors
- 6 corresponds to the min, mean, median, max, std, and cv of the PAWN indices
- \$T\$ is the number of time steps
"""
function tsa(X::DataFrame, y::AbstractMatrix{<:Real})::NamedDimsArray
    local ts
    try
        ts = axiskeys(y, 1)
    catch err
        if err isa MethodError
            ts = 1:size(y, 1)
        else
            rethrow(err)
        end
    end

    t_pawn_idx = NamedDimsArray(
        zeros(ncol(X), 6, size(y, 1));
        factors=Symbol.(names(X)),
        Si=[:min, :mean, :median, :max, :std, :cv],
        timesteps=ts
    )

    for t in axes(y, 1)
        # t_pawn_idx[:, :, t] .= col_normalize(
        #     pawn(X, vec(mean(y[1:t, :], dims=1)))
        # )
        t_pawn_idx[:, :, t] .= pawn_cvm(X, vec(mean(y[1:t, :], dims=1)))
    end

    return t_pawn_idx
end
function tsa(rs::ResultSet, y::AbstractMatrix{<:Real})::NamedDimsArray
    return tsa(rs.inputs, y)
end

"""
    rsa(X::DataFrame, y::Vector{<:Real}; S=10)::NamedDimsArray
    rsa(rs::ResultSet, y::AbstractArray{<:Real}; S::Int64=10)::NamedDimsArray

Perform Regional Sensitivity Analysis.

Regional Sensitivity Analysis is a Monte Carlo Filtering approach which aims to
identify which (group of) factors drive model outputs within or outside of a specified bound.
Outputs which fall inside the bounds are regarded as "behavioral", whereas those outside
are "non-behavioral". The distribution of behavioral/non-behavioral subsets are compared for
each factor. If the subsets are not similar, then the factor is influential. The sensitivity
index is simply the maximum distance between the two distributions, with larger values
indicating greater sensitivity.

The implemented approach slices factor space into \$S\$ bins and iteratively assesses
behavioral (samples within the bin) and non-behavioral (out of bin samples) subsets with the
non-parametric \$k\$-sample Anderson-Darling test. Larger values indicate greater
dissimilarity (thus, sensitivity). The Anderson-Darling test places more weight on the tails
compared to the Kolmogorov-Smirnov test.

RSA can indicate where in factor space model sensitivities may be, and contributes to a
Value-of-Information (VoI) analysis.

Increasing the value of \$S\$ increases the granularity of the analysis, but necessitates
larger sample sizes.

Note: Values of type `missing` indicate a lack of samples in the region.

# Arguments
- `rs` : ResultSet
- `X` : scenario specification
- `y` : scenario outcomes
- `S` : number of bins to slice factor space into (default: 10)

# Returns
NamedDimsArray, [bin values, factors]

# Examples
```julia
ADRIA.sensitivity.rsa(X, y; S=10)
```

# References
1. Pianosi, F., K. Beven, J. Freer, J. W. Hall, J. Rougier, D. B. Stephenson, and
   T. Wagener. 2016.
   Sensitivity analysis of environmental models:
   A systematic review with practical workflow.
   Environmental Modelling & Software 79:214-232.
   https://dx.doi.org/10.1016/j.envsoft.2016.02.008

2. Saltelli, A., M. Ratto, T. Andres, F. Campolongo, J. Cariboni, D. Gatelli,
   M. Saisana, and S. Tarantola. 2008.
   Global Sensitivity Analysis: The Primer.
   Wiley, West Sussex, U.K.
   https://dx.doi.org/10.1002/9780470725184
   Accessible at: http://www.andreasaltelli.eu/file/repository/Primer_Corrected_2022.pdf
"""
function rsa(X::DataFrame, y::AbstractVector{<:Real}; S::Int64=10)::NamedDimsArray
    N, D = size(X)
    seq = collect(0.0:(1/S):1.0)

    X_di = @MVector zeros(N)
    X_q = @MVector zeros(S + 1)
    r_s = zeros(Union{Missing,Float64}, S, D)
    sel = trues(N)

    for d_i in 1:D
        X_di .= X[:, d_i]
        X_q .= quantile(X_di, seq)

        sel .= X_q[1] .<= X_di .<= X_q[2]
        if count(sel) == 0 || length(y[Not(sel)]) == 0 || length(unique(y[sel])) == 1
            # not enough samples, or inactive area of factor space
            r_s[1, d_i] = missing
        else
            r_s[1, d_i] = KSampleADTest(y[sel], y[Not(sel)]).A²k
        end

        for s in 2:S
            sel .= X_q[s] .< X_di .<= X_q[s+1]
            if count(sel) == 0 || length(y[Not(sel)]) == 0 || length(unique(y[sel])) == 1
                # not enough samples, or inactive area of factor space
                r_s[s, d_i] = missing
                continue
            end

            # bs = bootstrap(mean, y[b], BalancedSampling(n_boot))
            # ci = confint(bs, PercentileConfInt(conf))[1]
            r_s[s, d_i] = KSampleADTest(y[sel], y[Not(sel)]).A²k
        end
    end

    return col_normalize(
        NamedDimsArray(r_s; bins=string.(seq[2:end]), factors=Symbol.(names(X)))
    )
end
function rsa(rs::ResultSet, y::AbstractVector{<:Real}; S::Int64=10)::NamedDimsArray
    return rsa(rs.inputs, vec(y); S=S)
end


"""
    outcome_map(X::DataFrame, y::AbstractVecOrMat, rule, target_factors::Vector; S::Int=20, n_boot::Int=100, conf::Float64=0.95)::NamedDimsArray

Map normalized outcomes (defined by `rule`) to factor values discretized into `S` bins.

Produces a matrix indicating the range of (normalized) outcomes across factor space for
each dimension (the model inputs). This is similar to a Regional Sensitivity Analysis,
except that the model outputs are examined directly as opposed to a measure of sensitivity.

Note:
- `y` is normalized on a per-column basis prior to the analysis
- Empty areas of factor space (those that do not have any desired outcomes)
  will be assigned `NaN`

# Arguments
- `X` : scenario specification
- `y` : Vector or Matrix of outcomes corresponding to scenarios in `X`
- `rule` : a callable defining a "desirable" scenario outcome
- `target_factors` : list of factors of interest to perform analyses on
- `S` : number of slices of factor space. Higher values equate to finer granularity
- `n_boot` : number of bootstraps (default: 100)
- `conf` : confidence interval (default: 0.95)

# Returns
3-dimensional NamedDimsArray, of shape \$S\$ ⋅ \$D\$ ⋅ 3, where:
- \$S\$ is the slices,
- \$D\$ is the number of dimensions, with
- boostrapped mean (dim 1) and the lower/upper 95% confidence interval (dims 2 and 3).

# Examples
```julia
# Get metric of interest
mu_tac = vec(mean(ADRIA.metrics.scenario_total_cover(rs), dims=:timesteps))

# Factors of interest
foi = [:SRM, :fogging, :a_adapt]

# Find scenarios where all metrics are above their median
rule = y -> all(y .> 0.5)

# Map input values where to their outcomes
ADRIA.sensitivity.outcome_map(X, y, rule, foi; S=20, n_boot=100, conf=0.95)
```
"""
function outcome_map(
    X::DataFrame,
    y::AbstractVecOrMat{<:Real},
    rule::Union{Function,BitVector,Vector{Int64}},
    target_factors::Vector{String},
    model_spec::DataFrame;
    S::Int64=20,
    n_boot::Int64=100,
    conf::Float64=0.95
)::NamedDimsArray
    if !all(target_factors .∈ [model_spec.fieldname])
        missing_factor = .!(target_factors .∈ [model_spec.fieldname])
        error("Invalid target factors: $(target_factors[missing_factor])")
    end

    factors_to_assess = model_spec.fieldname .∈ [target_factors]
    foi_attributes = Symbol[:fieldname, :ptype, :lower_bound, :upper_bound]
    foi_spec = model_spec[factors_to_assess, foi_attributes]

    foi_cat = (foi_spec.ptype .== "categorical")
    if any(foi_cat)
        max_bounds = maximum(
            foi_spec[foi_cat, :upper_bound] .-
            foi_spec[foi_cat, :lower_bound]
        )
        S = round(Int64, max(S, max_bounds))
    end

    step_size = 1 / S
    steps = collect(0:step_size:1.0)

    p_table = NamedDimsArray(
        zeros(Union{Missing,Float64}, length(steps) - 1, length(target_factors), 3);
        bins=string.(steps[2:end]),
        factors=Symbol.(target_factors),
        CI=[:mean, :lower, :upper]
    )

    all_p_rule = _map_outcomes(y, rule)
    if length(all_p_rule) == 0
        @warn "No results conform to specified rule."
        return p_table
    end

    # Identify behavioural
    n_scens = size(X, 1)
    behave::BitVector = falses(n_scens)
    behave[all_p_rule] .= true

    X_q = zeros(S + 1)
    for (j, fact_t) in enumerate(target_factors)
        ptype = model_spec.ptype[model_spec.fieldname .== fact_t][1]
        if ptype == "categorical"
            fact_idx = foi_spec.fieldname .== fact_t
            lb = foi_spec.lower_bound[fact_idx][1]
            ub = foi_spec.upper_bound[fact_idx][1]
            X_q .= round.(quantile(lb:ub, steps)) .- 1
        else
            X_q .= quantile(X[:, fact_t], steps)
        end

        for i in 1:length(X_q[1:end-1])
            local b::BitVector
            if i == 1
                b = (X_q[i] .<= X[:, fact_t] .<= X_q[i+1])
            else
                b = (X_q[i] .< X[:, fact_t] .<= X_q[i+1])
            end

            b = b .& behave

            if count(b) == 0
                # No data to bootstrap (empty region)
                p_table[i, j, [1, 2, 3]] .= missing
                continue
            end

            bs = bootstrap(mean, y[b], BalancedSampling(n_boot))
            ci = confint(bs, PercentileConfInt(conf))[1]

            p_table[i, j, [1, 2, 3]] .= ci
        end
    end

    return p_table
end
function outcome_map(
    X::DataFrame,
    y::AbstractVecOrMat{<:Real},
    rule::Union{Function,BitVector,Vector{Int64}};
    S::Int64=20,
    n_boot::Int64=100,
    conf::Float64=0.95
)::NamedDimsArray
    return outcome_map(X, y, rule, names(X); S, n_boot, conf)
end
function outcome_map(
    rs::ResultSet,
    y::AbstractArray{<:Real},
    rule::Union{Function,BitVector,Vector{Int64}},
    target_factors::Vector{String};
    S::Int64=20,
    n_boot::Int64=100,
    conf::Float64=0.95
)::NamedDimsArray
    return outcome_map(rs.inputs, y, rule, target_factors, rs.model_spec; S, n_boot, conf)
end
function outcome_map(
    rs::ResultSet,
    y::AbstractArray{<:Real},
    rule::Union{Function,BitVector,Vector{Int64}};
    S::Int64=20,
    n_boot::Int64=100,
    conf::Float64=0.95
)::NamedDimsArray
    return outcome_map(rs.inputs, y, rule, names(rs.inputs), rs.model_spec; S, n_boot, conf)
end

function _map_outcomes(
    y::AbstractArray{<:Real},
    rule::Union{BitVector,Vector{Int64}}
)::Union{BitVector,Vector{Int64}}
    return rule
end
function _map_outcomes(y::AbstractArray{<:Real}, rule::Function)::Vector{Int64}
    _y = col_normalize(y)
    all_p_rule = findall(rule, eachrow(_y))

    return all_p_rule
end


function cramervonmises_twosample(sample1::AbstractVector{T}, sample2::AbstractVector{T}) where {T}
    n1, n2 = length(sample1), length(sample2)
    combined_data = vcat(sample1, sample2)
    ranks = tiedrank(combined_data)

    ecdf1 = cumsum([count(x -> x <= x_i, sample1) / n1 for x_i in combined_data])
    ecdf2 = cumsum([count(x -> x <= x_i, sample2) / n2 for x_i in combined_data])

    W = (sum((ecdf1 - ecdf2) .^ 2 .* ranks) / (n1 * n2)) + (1 / (12 * (n1 + n2)))

    return W
end

"""
Regional sensitivity analysis

Apply RSA using CvM statistic across factor space comparing y_b | y
"""
function regional_cvm(X, y; S=10)
    cvm = zeros(S, size(X, 2))
    steps = collect(0.0:1/S:1.0)
    for d in axes(X, 2)
        qs = unique(quantile(X[:, d], steps))
        for (i, qi) in enumerate(qs[2:end])
            if i == 1
                xd = qs[i] .<= X[:, d] .< qi
            else
                xd = qs[i] .< X[:, d] .<= qi
            end

            cvm[i, d] = cramervonmises_twosample(y[xd], y)
            # cvm[i, d] = cramervonmises_twosample(y[xd], y[.!xd])
        end

        cvm[isnan.(cvm[:, d]), d] .= 0.0
    end

    return cvm
end


"""
PAWN-CvM

PAWN modified to use the Cramer-von Mises test instead of the usual Kolmogorov-Smirnoff
test.

Default behavior is to take the maximum difference for each factor.
"""
function pawn_cvm(X, y; S=10, stat=maximum)
    r_cvm = regional_cvm(X, y; S=S)
    cvm_i = stat.(eachcol(r_cvm))

    return (cvm_i ./ maximum(cvm_i))
end

end
