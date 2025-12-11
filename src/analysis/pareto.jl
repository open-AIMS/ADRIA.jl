"""
    find_pareto_optimal(rs::ResultSet, y::AbstractArray, rcps::Vector{Int} offset::Int=0)
    find_pareto_optimal(scens::DataFrame, y::AbstractArray, rcps::Vector{Int}; offset::Int=0)

Identify the scenarios (by position) that are determined to be pareto optimal.

# Arguments
- `rs` : ResultSet
- `y` : scenario outcomes
- `rcps` : RCP ids as integers (e.g., 45, 60, 85)
- `offset` : include scenarios that are some distance from pareto front, where 0 refers to
             the pareto front itself.

# Returns
NamedTuple, where each entry relates to an RCP of interest, e.g., `(RCP45=[... scenario ids ...], RCP60=[ ...scenario_ids ...])`

# Examples
```julia
using ADRIA, Statistics


rs = ADRIA.load_results("some result set")

tac = ADRIA.metrics.scenario_total_cover(rs)
rsv = ADRIA.metrics.scenario_rsv(rs)
r_juves = ADRIA.metrics.scenario_juveniles(rs)

# Create matrix of mean scenario outcomes
mean_tac = vec(mean(tac, dims=1))
mean_sv = vec(mean(rsv, dims=1))
mean_juves = vec(mean(r_juves, dims=1))
y = hcat(mean_tac, mean_sv, mean_juves)

# Find all pareto optimal scenarios
optimal = ADRIA.analysis.find_pareto_optimal(rs, y, [45, 60])
# (RCP45 = [13, 48, 54, 65, 95], RCP60 = [274, 315, 356, 430, 455])
```
"""
function find_pareto_optimal(
    scens::DataFrame, y::AbstractArray, rcps::Vector{Int}; offset::Int=0
)::NamedTuple
    x_idx = [scens.RCP .== rcp for rcp in rcps]
    r_rcp = [reduce(vcat, nds(y[rcp_idx, :], offset)) for rcp_idx in x_idx]

    scen_ids = [findall(rcp_idx)[r_idx] for (rcp_idx, r_idx) in zip(x_idx, r_rcp)]

    return NamedTuple{Tuple(Symbol.("RCP" .* string.(rcps)))}(scen_ids)
end
function find_pareto_optimal(
    rs::ResultSet, y::AbstractArray, rcps::Vector; offset::Int=0
)::NamedTuple
    return find_pareto_optimal(rs.inputs, y, rcps; offset=offset)
end

"""
    find_robust(rs::ResultSet, y::AbstractArray, rule, rcps::Vector{Int}; offset::Int=0)
    find_robust(scens::DataFrame, y::AbstractArray, rule, rcps::Vector{Int}; offset::Int=0)

Identify the scenarios (by position) that are determined to be robust and pareto optimal.

# Arguments
- `rs` : ResultSet
- `y` : scenario outcomes
- `rule` : a function
- `rcps` : RCP ids as integers (e.g., 45, 60, 85)
- `offset` : include scenarios that are some distance from pareto front, where 0 refers to
             the pareto front itself.

# Returns
NamedTuple, where each entry relates to an RCP of interest, e.g., `(RCP45=[... scenario ids ...], RCP60=[ ...scenario_ids ...])`

# Examples
```julia
using ADRIA, Statistics


rs = ADRIA.load_results("some result set")

tac = ADRIA.metrics.scenario_total_cover(rs)
rsv = ADRIA.metrics.scenario_rsv(rs)
r_juves = ADRIA.metrics.scenario_juveniles(rs)

# Create matrix of mean scenario outcomes
mean_tac = vec(mean(tac, dims=1))
mean_sv = vec(mean(rsv, dims=1))
mean_juves = vec(mean(r_juves, dims=1))
y = hcat(mean_tac, mean_sv, mean_juves)

# Find all pareto optimal scenarios where all metrics >= 0.9
rule_func = x -> all(x .>= 0.9)
robust = ADRIA.analysis.find_robust(rs, y, rule_func, [45, 60])
# (RCP45 = [13, 65], RCP60 = [274, 455])
```
"""
function find_robust(
    scens::DataFrame, y::AbstractArray, rule, rcps::Vector{Int64}; offset::Int64=0
)::NamedTuple
    y_star = col_normalize(copy(y))

    opt = find_pareto_optimal(scens, y_star, rcps; offset=offset)

    vals = Vector{Int64}[]
    sizehint!(vals, length(opt))
    for o in opt
        r = map(rule, eachrow(y_star[o, :]))
        if any(r .> 0)
            push!(vals, o[r])
            continue
        end

        push!(vals, Int64[])
    end

    return NamedTuple{keys(opt)}(vals)
end
function find_robust(
    rs::ResultSet, y::AbstractArray, rule, rcps::Vector{Int64}; offset::Int64=0
)::NamedTuple
    return find_robust(rs.inputs, y, rule, rcps; offset=offset)
end
