# ADRIAviz/test/plotly_fixtures.jl
#
# Synthetic data factories for the PlotlyBase backend test suite.
#
# CRITICAL: Do NOT import any Makie package (WGLMakie, GLMakie, CairoMakie,
# GeoMakie, GraphMakie) in this file or any file that includes it. Doing so
# would trigger ADRIAvizMakieExt and cause method-ambiguity failures.

using ADRIA
using ADRIA: DataCube, AnnotatedOutcomes
using OrderedCollections
using DataFrames
using Random
import ArchGDAL as AG

# ─────────────────────────────────────────────────────────────────────────────
# AnnotatedOutcomes fixtures
# ─────────────────────────────────────────────────────────────────────────────

"""
    _plotly_scenario_ao(; n_timesteps=10, n_scenarios=12) -> AnnotatedOutcomes

2-D AnnotatedOutcomes (timesteps × scenarios) with both type and RCP group
metadata. Mirrors `_scenario_ao()` in annotated_outcomes.jl but imports no
Makie package.
"""
function _plotly_scenario_ao(; n_timesteps=10, n_scenarios=12)
    data = DataCube(
        rand(n_timesteps, n_scenarios);
        timesteps=1:n_timesteps, scenarios=1:n_scenarios
    )
    sc_type_groups = OrderedDict{Symbol,BitVector}(
        :counterfactual => vcat(trues(4), falses(8)),
        :unguided => vcat(falses(4), trues(4), falses(4)),
        :guided => vcat(falses(8), trues(4))
    )
    metadata = Dict{Symbol,Any}(
        :scenario_type_groups => sc_type_groups,
        :scenario_rcp_groups => OrderedDict{Symbol,BitVector}(
            :rcp45 => vcat(trues(6), falses(6)),
            :rcp60 => vcat(falses(6), trues(6))
        )
    )
    return AnnotatedOutcomes(data, metadata)
end

"""
    _plotly_rme_ao(; n_timesteps=10, n_scenarios=8) -> AnnotatedOutcomes

AnnotatedOutcomes where `:scenario_rcp_groups` is `nothing` (RME-style result).
Used to verify that `by_RCP=true` raises an `ArgumentError`.
"""
function _plotly_rme_ao(; n_timesteps=10, n_scenarios=8)
    data = DataCube(
        rand(n_timesteps, n_scenarios);
        timesteps=1:n_timesteps, scenarios=1:n_scenarios
    )
    sc_type_groups = OrderedDict{Symbol,BitVector}(
        :counterfactual => vcat(trues(4), falses(4)),
        :guided => vcat(falses(4), trues(4))
    )
    return AnnotatedOutcomes(
        data,
        Dict{Symbol,Any}(
            :scenario_type_groups => sc_type_groups,
            :scenario_rcp_groups => nothing
        )
    )
end

"""
    _plotly_bare_ao(; n_timesteps=5, n_scenarios=4) -> AnnotatedOutcomes

AnnotatedOutcomes with an empty metadata Dict — no `:scenario_type_groups`.
Used to verify that `scenarios(ao)` raises an `ArgumentError` referencing
`attach_scenario_metadata`.
"""
function _plotly_bare_ao(; n_timesteps=5, n_scenarios=4)
    data = DataCube(
        rand(n_timesteps, n_scenarios);
        timesteps=1:n_timesteps, scenarios=1:n_scenarios
    )
    return AnnotatedOutcomes(data, Dict{Symbol,Any}())
end

# ─────────────────────────────────────────────────────────────────────────────
# Sensitivity analysis fixtures
# ─────────────────────────────────────────────────────────────────────────────

"""
    _plotly_pawn_si(; n_factors=4) -> YAXArray

2-D YAXArray (factors × Si_cols) matching the structure returned by
`ADRIAanalysis.sensitivity.pawn()`. Used for `pawn(Si)` tests.
"""
function _plotly_pawn_si(; n_factors=4)
    factor_names = Symbol.("Factor" .* string.(1:n_factors))
    Si_cols = [:min, :lb, :mean, :median, :ub, :max, :std, :cv]
    return DataCube(
        abs.(rand(n_factors, length(Si_cols)));
        factors=factor_names, Si=Si_cols
    )
end

"""
    _plotly_tsa_si(; n_factors=3, n_timesteps=10) -> YAXArray

3-D YAXArray (factors × Si_cols × timesteps) for Plotly `tsa(Si)` tests.
The Plotly backend accepts this signature without a `ResultSet` argument
(colour-by-component information from model_spec is not needed for web rendering).
"""
function _plotly_tsa_si(; n_factors=3, n_timesteps=10)
    factor_names = Symbol.("Factor" .* string.(1:n_factors))
    Si_cols = [:min, :lb, :mean, :median, :ub, :max]
    return DataCube(
        abs.(rand(n_factors, length(Si_cols), n_timesteps));
        factors=factor_names, Si=Si_cols, timesteps=1:n_timesteps
    )
end

"""
    _plotly_rsa_data(; n_factors=3, n_quantiles=10, n_scenarios=20)
    -> (Si::YAXArray, factor_values::Matrix{Float64})

2-D YAXArray (factors × si_quantile) and a scenarios × factors matrix of raw
factor values. Used for `rsa(Si, factor_values)` Plotly-specific signature
(no `ResultSet` needed).
"""
function _plotly_rsa_data(; n_factors=3, n_scenarios=20)
    factor_names = Symbol.("Factor" .* string.(1:n_factors))
    X = DataFrame(rand(n_scenarios, n_factors), string.(factor_names))
    y = rand(n_scenarios)
    foi = (factor_names[1], factor_names[2])
    return X, y, foi
end

"""
    _plotly_outcome_map_data(; n_factors=3, n_scenarios=20)
    -> (X::DataFrame, y::Vector{Float64}, factors::Vector{Symbol})

DataFrame of scenario inputs plus a scalar outcome per scenario and a list of
factor column names. Used for `outcome_map(X, y, factor)` and
`outcome_map(X, y, factors)` Plotly-specific signatures.
"""
function _plotly_outcome_map_data(; n_factors=3, n_scenarios=20)
    factor_names = Symbol.("Factor" .* string.(1:n_factors))
    X = DataFrame(rand(n_scenarios, n_factors), string.(factor_names))
    y = rand(n_scenarios)
    return X, y, factor_names
end

"""
    _plotly_convergence_si(; n_factors=3, n_sample_sizes=8)
    -> NamedTuple{(:Si, :factors), ...}

3-D YAXArray (factors × Si_cols × n_scenarios) for `convergence(Si_conv,
factors)`. Sample sizes are log-spaced between 50 and 500.
"""
function _plotly_convergence_si(; n_factors=3, n_sample_sizes=8)
    factor_names = Symbol.("Factor" .* string.(1:n_factors))
    Si_cols = [:lb, :median, :ub]
    sample_sizes = Int.(round.(exp.(range(log(50), log(500); length=n_sample_sizes))))
    Si = DataCube(
        abs.(rand(n_factors, length(Si_cols), n_sample_sizes));
        factors=factor_names, Si=Si_cols, n_scenarios=sample_sizes
    )
    return (Si=Si, factors=factor_names)
end

# ─────────────────────────────────────────────────────────────────────────────
# Clustering fixtures
# ─────────────────────────────────────────────────────────────────────────────

"""
    _plotly_cluster_data(; n_timesteps=10, n_scenarios=12, n_clusters=3)
    -> (matrix, clusters::Vector{Int})

Random outcomes matrix (timesteps × scenarios) and integer cluster assignments
for `scenarios(matrix, clusters)` and `clustered_scenarios(matrix, clusters)`.
"""
function _plotly_cluster_data(; n_timesteps=10, n_scenarios=12, n_clusters=3)
    outcomes = rand(n_timesteps, n_scenarios)
    clusters = repeat(1:n_clusters, cld(n_scenarios, n_clusters))[1:n_scenarios]
    return outcomes, clusters
end

"""
    _plotly_bitvec_cluster_data(; n_timesteps=10, n_scenarios=12)
    -> (matrix, clusters::BitVector)

BitVector version of cluster assignments (target / non-target binary split).
"""
function _plotly_bitvec_cluster_data(; n_timesteps=10, n_scenarios=12)
    outcomes = rand(n_timesteps, n_scenarios)
    half = n_scenarios ÷ 2
    clusters = vcat(trues(half), falses(n_scenarios - half))
    return outcomes, clusters
end

# ─────────────────────────────────────────────────────────────────────────────
# Data envelopment analysis fixture
# ─────────────────────────────────────────────────────────────────────────────

"""
    _plotly_dea_output(; n_scenarios=20) -> NamedTuple

Named tuple that mirrors the fields accessed by
`ADRIA.viz.data_envelopment_analysis`:
- `Y`         : n_scenarios × 2 output metric matrix
- `vrs_vals`  : VRS efficiency scores (∈ [0, 1])
- `crs_vals`  : CRS efficiency scores (≤ VRS)
- `vrs_peers` : `(J = [indices],)` of efficient scenarios under VRS
- `crs_peers` : `(J = [indices],)` of efficient scenarios under CRS
"""
function _plotly_dea_output(; n_scenarios=20)
    Y = hcat(rand(n_scenarios), rand(n_scenarios))
    vrs_vals = clamp.(rand(n_scenarios) .+ 0.3, 0.0, 1.0)
    crs_vals = clamp.(vrs_vals .- abs.(rand(n_scenarios) .* 0.2), 0.0, 1.0)
    efficient = findall(vrs_vals .>= 0.9)
    isempty(efficient) && (efficient = [argmax(vrs_vals)])
    return (
        Y=Y,
        vrs_vals=vrs_vals,
        crs_vals=crs_vals,
        vrs_peers=(J=efficient,),
        crs_peers=(J=efficient,)
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Rule extraction fixtures
# ─────────────────────────────────────────────────────────────────────────────

"""
    _plotly_taxonomy_df(; n_timesteps=10, n_groups=5, n_scenarios=12)
    -> (scens::DataFrame, rtc::YAXArray)

Synthetic scenarios DataFrame and relative_taxa_cover YAXArray for
`taxonomy(scenarios, relative_taxa_cover)` Plotly tests.
"""
function _plotly_taxonomy_df(; n_timesteps=10, n_groups=5, n_scenarios=12)
    scens = _plotly_scenarios_df(; n_scenarios=n_scenarios)
    rtc = DataCube(
        rand(n_timesteps, n_groups, n_scenarios);
        timesteps=1:n_timesteps, groups=1:n_groups, scenarios=1:n_scenarios
    )
    return scens, rtc
end

"""
    _plotly_scenarios_rs()

Attempts to load a minimal ResultSet from the ADRIA test domain.
Returns `nothing` if the test domain is not present (CI skip guard).
Used by the `scenarios(rs, outcomes)` tests.
"""
function _plotly_scenarios_rs()
    test_data_dir = joinpath(pkgdir(ADRIA), "test", "data")
    domain_path = joinpath(test_data_dir, "Test_domain")
    isdir(domain_path) || return nothing
    try
        dom = ADRIA.load_domain(domain_path, "45")
        scens = ADRIA.sample(dom, 4)
        return ADRIA.run_scenarios(dom, scens, "45")
    catch
        return nothing
    end
end

"""
    _plotly_taxonomy_rs()

Attempts to load a minimal ResultSet from the ADRIA test domain.
Returns `nothing` if the test domain is not present (CI skip guard).
"""
function _plotly_taxonomy_rs()
    test_data_dir = joinpath(pkgdir(ADRIA), "test", "data")
    domain_path = joinpath(test_data_dir, "Test_domain")
    isdir(domain_path) || return nothing
    try
        dom = ADRIA.load_domain(domain_path, "45")
        scens = ADRIA.sample(dom, 4)
        return ADRIA.run_scenarios(dom, scens, "45")
    catch
        return nothing
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Rule extraction fixtures
# ─────────────────────────────────────────────────────────────────────────────

"""
    _plotly_rules(; n_rules=4) -> Vector{<:NamedTuple}

Rules mirroring `ADRIAanalysis.Rule` (a `.condition` field holding a list of
2-element clauses, each a 3-element vector `[feature_name, operator, threshold]`).
Built as NamedTuples rather than `ADRIAanalysis.Rule` directly since ADRIAviz's
test environment does not depend on ADRIAanalysis.
"""
function _plotly_rules(; n_rules=4)
    features = ["N_seed_TA", "fogging"]
    return [
        (;
            condition=[
                [features[1], "<", 0.3 + 0.05 * i],
                [features[2], "<=", 0.2 + 0.05 * i]
            ]
        )
        for i = 1:n_rules
    ]
end

"""
    _plotly_scenarios_df(; n_scenarios=20) -> DataFrame

Minimal scenario DataFrame with the columns expected by `rules_scatter`.
"""
function _plotly_scenarios_df(; n_scenarios=20)
    return DataFrame(;
        N_seed_TA=rand(n_scenarios),
        fogging=rand(n_scenarios),
        guided=rand(0:3, n_scenarios),
        SRM=zeros(n_scenarios),
        N_mc_settlers=zeros(Int, n_scenarios),
        mcb_duration=zeros(Int, n_scenarios),
        RCP=rand([45, 60], n_scenarios)
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Spatial fixtures
# ─────────────────────────────────────────────────────────────────────────────

"""
    _plotly_spatial_gdf(; n_sites=5) -> DataFrame

Synthetic GeoDataFrame with small square WGS84 polygon "reef patches" around
GBR coordinates. Includes `:site_id`, `:k`, and `:geometry` columns.
ArchGDAL is safe to import here — it does not trigger any Makie extension.
"""
# ─────────────────────────────────────────────────────────────────────────────
# Environment fixtures
# ─────────────────────────────────────────────────────────────────────────────

"""
    _plotly_cyclone_scens(; n_timesteps=10, n_locs=5, n_species=3, n_scens=4)
    -> YAXArray{Float64,4}

Synthetic cyclone mortality scenarios: `(timesteps × locations × species × scenarios)`.
Values are in [0, 1] (fractional mortality); the viz function multiplies by 100.
"""
function _plotly_cyclone_scens(; n_timesteps=10, n_locs=5, n_species=3, n_scens=4)
    return DataCube(
        rand(n_timesteps, n_locs, n_species, n_scens);
        timesteps=1:n_timesteps,
        locations=1:n_locs,
        species=1:n_species,
        scenarios=1:n_scens
    )
end

"""
    _plotly_dhw_scens(; n_timesteps=10, n_sites=5, n_scens=4) -> YAXArray{Float64,3}

Synthetic DHW scenarios: `(timesteps × sites × scenarios)`.
"""
function _plotly_dhw_scens(; n_timesteps=10, n_sites=5, n_scens=4)
    return DataCube(
        rand(n_timesteps, n_sites, n_scens) .* 8.0;   # DHW values typically 0–8
        timesteps=1:n_timesteps,
        sites=1:n_sites,
        scenarios=1:n_scens
    )
end

function _plotly_spatial_gdf(; n_sites::Int=5)
    h = 0.01  # half-side ~1 km at GBR latitudes
    geoms = map(1:n_sites) do i
        cx = 147.0 + i * 0.05
        cy = -18.0 + i * 0.05
        wkt = (
            "POLYGON((" *
            "$(cx - h) $(cy - h), $(cx + h) $(cy - h), " *
            "$(cx + h) $(cy + h), $(cx - h) $(cy + h), " *
            "$(cx - h) $(cy - h)))"
        )
        return AG.fromWKT(wkt)
    end
    return DataFrame(;
        site_id=string.(1:n_sites),
        k=fill(0.3, n_sites),
        geometry=geoms
    )
end
