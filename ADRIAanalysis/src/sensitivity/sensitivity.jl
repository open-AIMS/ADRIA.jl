module sensitivity

import ..feature_set

using Logging

using
    StaticArrays,
    YAXArrays
using DataFrames

using
    Bootstrap,
    Distributions,
    HypothesisTests,
    Random
using HypothesisTests: ApproximateKSTest

using ADRIA: DataCube, ResultSet, model_spec, ZeroDataCube

using ADRIA.analysis: col_normalize, normalize!

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
summary statistics (min, lower bound, mean, median, upper bound, max, std, and cv) over the slices.

# Arguments
- `rs` : ResultSet
- `X` : Model inputs
- `y` : Model outputs
- `factor_names` : Names of each factor represented by columns in `X`
- `S` : Number of slides (default: 10)

# Returns
YAXArray, of min, mean, lower bound, median, upper bound, max, std, and cv summary statistics.

# Examples
```julia
dom = ADRIA.load_domain("example_domain", "<RCP>")
scens = ADRIA.sample(dom, 128)
rs = ADRIA.run_scenarios(dom, scens, "45")

# Get mean coral cover over time and locations
μ_tac = mean(ADRIA.metrics.scenario_total_cover(rs), dims=:timesteps)

ADRIAanalysis.sensitivity.pawn(rs, μ_tac)
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
)::YAXArray
    N, D = size(X)
    step = 1 / S
    seq = 0.0:step:1.0

    # Preallocate result structures
    X_q = @MVector zeros(S + 1)
    pawn_t = @MArray zeros(S, D)
    results = @MArray zeros(D, 8)
    q_stats = [0.025, 0.5, 0.975]

    # Hide warnings from HypothesisTests
    with_logger(NullLogger()) do
        for d_i = 1:D
            X_di = @view(X[:, d_i])
            X_q .= quantile(X_di, seq)

            Y_sel = @view(y[X_q[1] .<= X_di .<= X_q[2]])
            if length(Y_sel) > 0
                pawn_t[1, d_i] = ks_statistic(ApproximateTwoSampleKSTest(Y_sel, y))
            end

            for s = 2:S
                Y_sel = @view(y[X_q[s] .< X_di .<= X_q[s + 1]])
                if length(Y_sel) == 0
                    continue  # no available samples
                end

                pawn_t[s, d_i] = ks_statistic(ApproximateTwoSampleKSTest(Y_sel, y))
            end

            p_ind = @view(pawn_t[:, d_i])
            p_mean = mean(p_ind)
            p_sdv = std(p_ind)
            p_cv = p_sdv ./ p_mean
            p_lb, p_med, p_ub = quantile(p_ind, q_stats)
            results[d_i, :] .= (
                minimum(p_ind),
                p_lb,
                p_mean,
                p_med,
                p_ub,
                maximum(p_ind),
                p_sdv,
                p_cv
            )
        end
    end

    replace!(results, NaN => 0.0, Inf => 0.0)

    col_names = [:min, :lb, :mean, :median, :ub, :max, :std, :cv]
    row_names = Symbol.(factor_names)
    return DataCube(results; factors=row_names, Si=col_names)
end
function pawn(X::DataFrame, y::AbstractVector{<:Real}; S::Int64=10)::YAXArray
    return pawn(Matrix(X), y, names(X); S=S)
end
function pawn(
    X::YAXArray,
    y::AbstractVector{<:Real};
    S::Int64=10
)::YAXArray
    return pawn(X, y, collect(X.axes[2]); S=S)
end
function pawn(
    X::YAXArray,
    y::YAXArray;
    S::Int64=10
)::YAXArray
    # YAXrrays will raise an error if any of the masked boolean indexing in pawn is empty so
    # vec(y) is required
    return pawn(X, vec(y), collect(X.axes[2]); S=S)
end
function pawn(
    X::Union{DataFrame,YAXArray},
    y::AbstractMatrix{<:Real};
    S::Int64=10
)::YAXArray
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
    y::Union{YAXArray,AbstractVector{<:Real}};
    S::Int64=10
)::YAXArray
    return pawn(rs.inputs, y; S=S)
end

"""
    convergence(X::DataFrame, y::YAXArray, target_factors::Vector{Symbol}; n_steps::Int64=10)::YAXArray
    convergence(rs::ResultSet, X::DataFrame, y::YAXArray, components::Vector{String}; n_steps::Int64=10)::YAXArray

Calculates the PAWN sensitivity index for an increasing number of scenarios where the
maximum is the total number of scenarios in scens. Number of scenario subsets determined by
N_steps. Can be calculated for individual factors or aggregated over factors for specified
model components.

# Arguments
- `rs` : Result set (only needed if aggregating over model components).
- `X` : Model inputs
- `y` : Model outputs
- `target_factors` : Names of target factors represented by columns in `X`.
- `components` : Names of model components to aggregate over (e.g. [:Intervention, :Criteria]).
- `n_steps` : Number of steps to cut the total number of scenarios into.

# Returns
YAXArray, of min, lower bound, mean, median, upper bound, max, std, and cv summary
statistics for an increasing number of scenarios.
"""
function convergence(
    X::DataFrame,
    y::YAXArray,
    target_factors::Vector{Symbol};
    Si::Function=pawn,
    n_steps::Int64=10
)::YAXArray
    N = length(y.scenarios)
    step_size = floor(Int64, N / n_steps)
    N_it = step_size == 0 ? collect(1:N) : collect(step_size:step_size:N)

    pawn_store = ZeroDataCube(;
        T=Float64,
        factors=target_factors,
        Si=[:min, :lb, :mean, :median, :ub, :max, :std, :cv],
        n_scenarios=N_it
    )
    scens_idx = randperm(N)

    for nn in N_it
        pawn_store[n_scenarios = At(nn)] .= Array(
            Si(
                X[scens_idx[1:nn], :], Array(y[scens_idx[1:nn]])
            )[factors = At(target_factors)]
        )
    end

    return pawn_store
end
function convergence(
    rs::ResultSet,
    X::DataFrame,
    y::YAXArray,
    components::Vector{Symbol};
    Si::Function=pawn,
    n_steps::Int64=10
)::YAXArray
    ms = model_spec(rs, Symbol.(names(X)))

    target_factors = [
        ms[ms[:, "component"] .== cc, "fieldname"] for
        cc in string.(components)
    ]

    Si_n = convergence(X, y, Symbol.(vcat(target_factors...)); Si=Si, n_steps=n_steps)

    # Note: n_steps only applies if it is > number of scenarios.
    Si_grouped = ZeroDataCube(;
        T=Float64,
        factors=components,
        Si=collect(Si_n.Si),
        n_scenarios=collect(Si_n.n_scenarios)
    )

    for (cc, factors) in zip(components, target_factors)
        Si_grouped[factors = At(cc)] .= dropdims(
            mean(Si_n[factors = At(Symbol.(factors))]; dims=:factors);
            dims=:factors
        )
    end

    return Si_grouped
end

"""
    tsa(X::DataFrame, y::AbstractMatrix)::YAXArray

Perform Temporal (or time-varying) Sensitivity Analysis using the PAWN sensitivity index.

The sensitivity index value for time \$t\$ is inclusive of all time steps prior to \$t\$.
Alternate approaches use a moving window, or only data for time \$t\$.

# Examples
```julia
rs = ADRIA.load_results("a ResultSet of interest")

# Get scenario outcomes over time (shape: `time × scenarios`)
y_tac = ADRIA.metrics.scenario_total_cover(rs)

# Calculate sensitivity of outcome to factors for each time step
ADRIAanalysis.sensitivity.tsa(rs.inputs, y_tac)
```

# Arguments
- `X` : Scenario specification
- `y` : scenario outcomes over time

# Returns
YAXArray, of shape \$D\$ × 6 × \$T\$, where
- \$D\$ is the number of dimensions/factors
- 6 corresponds to the min, mean, median, max, std, and cv of the PAWN indices
- \$T\$ is the number of time steps
"""
function tsa(X::DataFrame, y::AbstractMatrix{<:Real})::YAXArray
    ts = hasproperty(y, :axes) ? collect(y.axes[1]) : (1:size(y, 1))

    t_pawn_idx = ZeroDataCube(;
        T=Float64,
        factors=Symbol.(names(X)),
        Si=[:min, :lb, :mean, :median, :ub, :max, :std, :cv],
        timesteps=ts
    )

    for t in axes(y, 1)
        t_pawn_idx[:, :, t] .= col_normalize(
            pawn(X, vec(mean(y[1:t, :]; dims=1))).data
        )
    end

    return t_pawn_idx
end
function tsa(rs::ResultSet, y::AbstractMatrix{<:Real})::YAXArray
    return tsa(rs.inputs, y)
end

include("rsa.jl")

end
