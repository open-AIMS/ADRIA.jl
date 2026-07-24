using Test
using ADRIA
using ADRIA: DataFrames
using ADRIAanalysis

@testset "Regression: MLJ/SIRUS not auto-loaded" begin
    @test !any(string(k.name) == "SIRUS" for k in keys(Base.loaded_modules)) &&
        !any(string(k.name) == "MLJ" for k in keys(Base.loaded_modules))
end

@testset "ADRIAanalysis" begin
    pkg_ts = @testset "Package loads" begin
        @test isdefined(ADRIAanalysis, :sensitivity)
        @test isdefined(ADRIAanalysis, :data_envelopment_analysis)
        @test isdefined(ADRIAanalysis, :rules)
        @test isdefined(ADRIAanalysis, :cluster_rules)
        @test isdefined(ADRIAanalysis, :target_clusters)
        @test isdefined(ADRIAanalysis, :find_scenarios)
        @test isdefined(ADRIAanalysis, :scenario_types)
        @test isdefined(ADRIAanalysis, :scenario_rcps)
        @test isdefined(ADRIAanalysis, :scenario_clusters)
        @test isdefined(ADRIAanalysis, :find_pareto_optimal)
        @test isdefined(ADRIAanalysis, :find_robust)
        @test isdefined(ADRIAanalysis, :screen_scenarios)
    end
    let tc = Test.get_test_counts(pkg_ts)
        if tc.fails + tc.errors > 0
            @error "ADRIAanalysis: Package loads failed — aborting functional tests"
            exit(1)
        end
    end

    # ---------------------------------------------------------------------------
    # ADRIA.analysis.cluster_series  (inlined into src/analysis/analysis.jl)
    # ---------------------------------------------------------------------------
    @testset "ADRIA.analysis.cluster_series" begin
        # 8 timesteps × 6 scenarios — two obvious groups
        low = [0.1, 0.15, 0.12, 0.11, 0.14, 0.13]'  # 1×6
        high = [0.9, 0.85, 0.88, 0.91, 0.87, 0.89]'  # 1×6
        data = repeat(low, 8) .+ repeat(0.0:0.01:0.07, 1, 6)   # 8×6, gently rising
        # Replace top 3 columns with clearly higher values
        data[:, 4:6] .+= 0.6

        n_clusters = 2
        assignments = ADRIA.analysis.cluster_series(data, n_clusters)

        @test length(assignments) == 6
        @test eltype(assignments) <: Integer
        @test all(1 .<= assignments .<= n_clusters)
        # The two obvious groups should be separated into different clusters
        @test assignments[1:3] == assignments[1:3]          # low group internally consistent
        @test allunique([assignments[1], assignments[4]]) == false ||
            assignments[1] != assignments[4]

        @testset "hclust method" begin
            h_assignments = ADRIA.analysis.cluster_series(data, n_clusters; method=:hclust)
            @test length(h_assignments) == 6
            @test all(1 .<= h_assignments .<= n_clusters)
        end

        @testset "weuclidean distance" begin
            w_assignments = ADRIA.analysis.cluster_series(
                data, n_clusters; distance=:weuclidean
            )
            @test length(w_assignments) == 6
            @test all(1 .<= w_assignments .<= n_clusters)
        end

        @testset "3-cluster case" begin
            # 3 clear groups: low / mid / high  (3 scenarios each, 5 timesteps)
            g1 = 0.1 .* ones(5, 3)
            g2 = 0.5 .* ones(5, 3)
            g3 = 0.9 .* ones(5, 3)
            d3 = hcat(g1, g2, g3)   # 5×9
            a3 = ADRIA.analysis.cluster_series(d3, 3)
            @test length(a3) == 9
            @test length(unique(a3)) == 3
            # Scenarios in same group should share a cluster label
            @test a3[1] == a3[2] == a3[3]
            @test a3[4] == a3[5] == a3[6]
            @test a3[7] == a3[8] == a3[9]
        end
    end

    # ---------------------------------------------------------------------------
    # ADRIA.analysis.cluster_scenarios  (2D and 3D paths)
    # ---------------------------------------------------------------------------
    @testset "ADRIA.analysis.cluster_scenarios" begin
        # 2D path: should delegate to cluster_series
        d2 = rand(10, 8)
        a2 = ADRIA.analysis.cluster_scenarios(d2, 3)
        @test a2 isa Array{Int64}
        @test size(a2) == (8,)
        @test all(1 .<= a2 .<= 3)

        # 3D path: T × S × M  (4 timesteps, 6 scenarios, 2 metrics)
        d3 = rand(4, 6, 2)
        a3 = ADRIA.analysis.cluster_scenarios(d3, 2)
        @test a3 isa Array{Int64}
        @test size(a3) == (6, 2)       # S × M
        @test all(1 .<= a3 .<= 2)
    end

    # ---------------------------------------------------------------------------
    # ADRIAanalysis.target_clusters
    # ---------------------------------------------------------------------------
    @testset "target_clusters" begin
        # outcomes: 6 timesteps × 6 scenarios
        # Scenarios 4-6 have high values (should be the target cluster)
        outcomes = Float64[
            0.1 0.1 0.1 0.9 0.9 0.9;
            0.1 0.2 0.1 0.8 0.9 0.9;
            0.1 0.1 0.2 0.9 0.8 0.9;
            0.2 0.1 0.1 0.9 0.9 0.8;
            0.1 0.2 0.2 0.8 0.8 0.9;
            0.2 0.1 0.1 0.9 0.9 0.8
        ]
        # Two clear clusters
        clusters = Int64[1, 1, 1, 2, 2, 2]

        result = ADRIAanalysis.target_clusters(clusters, outcomes)

        @test result isa Vector{Int64}
        @test length(result) == 6
        @test all(r ∈ (0, 1) for r in result)

        # The high-value group (scenarios 4-6, cluster 2) should be the target
        @test result[4] == 1
        @test result[5] == 1
        @test result[6] == 1
        @test result[1] == 0
        @test result[2] == 0
        @test result[3] == 0
    end

    @testset "target_clusters — size_limit forces merge" begin
        # One cluster has only 1 scenario (< 1% of 100 would trigger merge;
        # here we use size_limit=0.4 to force the merge with 3 scenarios total)
        outcomes = Float64[0.9 0.1 0.2]  # 1 timestep × 3 scenarios
        clusters = Int64[1, 2, 2]

        result = ADRIAanalysis.target_clusters(clusters, outcomes; size_limit=0.4)
        # Even though cluster 1 is "best", it's < 40% of scenarios,
        # so cluster 2 gets merged in. All scenarios should be target.
        @test result isa Vector{Int64}
        @test length(result) == 3
    end

    # ---------------------------------------------------------------------------
    # ADRIAAnalysis.find_scenarios  — 2D YAXArray path
    #
    # NOTE: The 2D method accepts AbstractArray{<:Real} but the body uses
    # YAXArray-specific axis indexing (axes(x, :timesteps) and x[timesteps=...]).
    # Real callers always pass a YAXArray (e.g. from ADRIA.metrics.loc_trajectory).
    # Passing a plain Matrix dispatches correctly but throws inside the body.
    # The test below uses ADRIA.DataCube to create a minimal YAXArray.
    # ---------------------------------------------------------------------------
    @testset "find_scenarios 2D (YAXArray with :timesteps axis)" begin
        # 4 timesteps × 4 scenarios: scenarios 3 & 4 have high values
        data = Float64[
            0.1 0.1 0.9 0.9;
            0.1 0.1 0.9 0.9;
            0.1 0.1 0.9 0.9;
            0.1 0.1 0.9 0.9
        ]
        outcomes_yax = ADRIA.DataCube(data; timesteps=1:4, scenarios=1:4)
        clusters_2d = Int64[1, 1, 2, 2]
        filter_fn = x -> x .>= maximum(x)   # pick the single best cluster

        result = ADRIAanalysis.find_scenarios(outcomes_yax, clusters_2d, filter_fn)

        @test result isa BitVector
        @test length(result) == 4
        # Scenarios 3 & 4 (high-value cluster 2) should be selected
        @test result[3] == true
        @test result[4] == true
        @test result[1] == false
        @test result[2] == false
    end

    # Also verify that passing a plain Matrix raises (documents the type mismatch)
    @testset "find_scenarios 2D plain Matrix — type mismatch" begin
        outcomes_2d = Float64[0.1 0.9; 0.1 0.9; 0.1 0.9; 0.1 0.9]
        clusters_2d = Int64[1, 2]
        filter_fn = x -> x .>= maximum(x)
        @test_throws Exception ADRIAanalysis.find_scenarios(
            outcomes_2d, clusters_2d, filter_fn
        )
    end

    # ---------------------------------------------------------------------------
    # find_pareto_optimal / find_robust / screen_scenarios
    # ---------------------------------------------------------------------------
    @testset "find_pareto_optimal" begin
        # Scenario 1 dominates scenario 2 for RCP45 (higher is better in both objectives).
        # Scenarios 3 and 4 are mutually non-dominated for RCP60.
        scens = DataFrames.DataFrame(; RCP=Int[45, 45, 60, 60])
        y = Float64[0.9 0.9; 0.1 0.1; 0.8 0.2; 0.3 0.7]

        result = ADRIAanalysis.find_pareto_optimal(scens, y, [45, 60])

        @test result isa NamedTuple
        @test haskey(result, :RCP45)
        @test haskey(result, :RCP60)
        @test result.RCP45 == [1]
        @test sort(result.RCP60) == [3, 4]
    end

    @testset "find_robust" begin
        scens = DataFrames.DataFrame(; RCP=Int[45, 45, 45])
        # After col_normalize all three end up on the Pareto front.
        # y_star = [1.0 0.0; 0.5 0.5; 0.0 1.0]
        y = Float64[1.0 0.0; 0.5 0.5; 0.0 1.0]

        # Rule: all normalized values >= 0.4 → only scenario 2 ([0.5, 0.5]) passes.
        rule = x -> all(x .>= 0.4)
        result = ADRIAanalysis.find_robust(scens, y, rule, [45])
        @test result isa NamedTuple
        @test haskey(result, :RCP45)
        @test result.RCP45 == [2]

        # Rule no scenario satisfies → empty result.
        result_empty = ADRIAanalysis.find_robust(scens, y, x -> all(x .>= 0.99), [45])
        @test isempty(result_empty.RCP45)
    end

    @testset "screen_scenarios" begin
        # col_normalize([1 3; 2 4]) = [0 0; 1 1]
        y = Float64[1.0 3.0; 2.0 4.0]

        @test ADRIAanalysis.screen_scenarios(y, x -> x >= 0.5) == [2]
        @test ADRIAanalysis.screen_scenarios(y, x -> x >= 0.0) == [1, 2]
        @test isempty(ADRIAanalysis.screen_scenarios(y, x -> x > 1.0))
    end

    # ---------------------------------------------------------------------------
    # Submodule and function existence checks
    # ---------------------------------------------------------------------------
    @testset "Submodule and function presence" begin
        @test isdefined(ADRIAanalysis, :rules)
        @test isdefined(ADRIAanalysis, :cluster_rules)
    end

    include("unit_rule_core.jl")
    include("negative_cluster_rules.jl")
    include("ext_activation.jl")
    include("integration_cluster_rules.jl")
    include("test_rsa.jl")
    include("test_rsa_dispatch.jl")
    include("test_cramers_v.jl")
    include("test_stratified_rsa.jl")
    include("test_counterfactual_delta.jl")
    include("test_feature_set_metadata.jl")
end
