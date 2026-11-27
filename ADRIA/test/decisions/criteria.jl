using Test
using LinearAlgebra
using Random
using SparseArrays

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

@testset "Decision criteria" begin
    @testset "Criterion" begin
        c = ADRIA.decision.Criterion(:heat_stress)
        @test c.name == :heat_stress
        @test c.spec.val == 1.0
        @test c.spec.direction === minimum

        c = ADRIA.decision.Criterion(:heat_stress; val=0.9, dist_params=(0.8, 1.0))
        @test c.spec.val == 0.9
        @test c.spec.dist_params == (0.8, 1.0)

        @test_throws ArgumentError ADRIA.decision.Criterion(:not_a_criterion)
        @test_throws ArgumentError ADRIA.decision.Criterion(:heat_stress; name="Other")
        @test_throws ArgumentError ADRIA.decision.Criterion(:heat_stress; description="Other")
    end

    @testset "Weights are generated from the criteria lists" begin
        for (W, prefix) in (
            (ADRIA.CAqCriteriaWeights, "iv_CAq"),
            (ADRIA.FogCriteriaWeights, "iv_Fog"),
            (ADRIA.LvMCriteriaWeights, "iv_LvM")
        )
            names = ADRIA.decision.criterion_names(W)
            fields = ADRIA.decision.criterion_fields(W)

            @test collect(fieldnames(W)) == fields
            @test fields == Symbol.(prefix, "_", names)
            @test :out_spread in names
            @test names[findfirst(==(:out_spread), names) - 1] == :out_connectivity
        end

        listed = reduce(
            vcat,
            ADRIA.decision.criterion_fields.([
                ADRIA.CAqCriteriaWeights, ADRIA.FogCriteriaWeights, ADRIA.LvMCriteriaWeights
            ])
        )
        @test ADRIA._GROUP_MEMBERS[:criteria_weights] == listed
        @test all(startswith.(string.(ADRIA._CRITERIA_WEIGHT_GROUPS[:CAq_weights]), "iv_CAq_"))
    end

    @testset "out_spread weights" begin
        for (W, val, range) in (
            (ADRIA.CAqCriteriaWeights, 0.8, [0.2, 1.0]),
            (ADRIA.FogCriteriaWeights, 1.0, [0.0, 1.0]),
            (ADRIA.LvMCriteriaWeights, 0.8, [0.5, 1.0])
        )
            weights = W()
            p = getfield(weights, only(filter(f -> endswith(string(f), "_out_spread"), fieldnames(W))))

            @test p.val == val
            @test collect(p.dist_params) == range
            @test p.direction === maximum
            @test startswith(p.name, "Outgoing Spread (")
            @test endswith(p.description, "Prefers higher values.")
        end
    end

    @testset "Labels state the direction in use" begin
        fog_depth = ADRIA.FogCriteriaWeights().iv_Fog_depth
        caq_depth = ADRIA.CAqCriteriaWeights().iv_CAq_depth

        @test endswith(fog_depth.description, "Prefers lower values.")
        @test endswith(caq_depth.description, "Prefers higher values.")
    end

    @testset "out_spread" begin
        # Independent reference: Katz centrality of the matrix without its diagonal
        function reference(C; fraction=0.5)
            C0 = Matrix{Float64}(C)
            C0[diagind(C0)] .= 0.0
            rho = maximum(abs, eigvals(C0))
            katz = (I - (fraction / rho) * C0) \ ones(size(C0, 1)) .- 1.0
            return katz ./ maximum(katz)
        end

        C = [
            0.3 0.4 0.0 0.1
            0.0 0.2 0.5 0.0
            0.2 0.0 0.1 0.3
            0.0 0.0 0.0 0.6
        ]
        s = ADRIA.out_spread(C)

        @test s ≈ reference(C)
        @test maximum(s) == 1.0
        @test all(0.0 .<= s .<= 1.0)

        @testset "Self links are ignored" begin
            C_diag = copy(C)
            C_diag[diagind(C_diag)] .= [0.9, 0.0, 0.5, 0.2]
            @test ADRIA.out_spread(C_diag) ≈ s
        end

        @testset "Scale of the matrix does not matter" begin
            @test ADRIA.out_spread(10.0 .* C) ≈ s
            @test ADRIA.out_spread(Float32.(C)) ≈ s rtol = 1e-5
        end

        @testset "Locations that only send larvae to themselves score zero" begin
            @test isapprox(s[4], 0.0; atol=1e-12)
            @test all(s[1:3] .> 0.0)
        end

        @testset "No links to other locations" begin
            @test ADRIA.out_spread(zeros(3, 3)) == zeros(3)
            @test ADRIA.out_spread(Matrix(0.5I, 3, 3)) == zeros(3)
            @test ADRIA.out_spread(fill(0.1, 1, 1)) == [0.0]
        end

        @testset "Chain without cycles" begin
            chain = [0.0 0.4 0.0; 0.0 0.0 0.4; 0.0 0.0 0.0]
            s_chain = ADRIA.out_spread(chain)

            @test all(isfinite, s_chain)
            @test s_chain[1] > s_chain[2] > s_chain[3]
            @test isapprox(s_chain[3], 0.0; atol=1e-12)
        end

        @testset "Sending larvae to well connected locations scores higher" begin
            # 1 -> 2 -> 5 and 3 -> 4, where 4 is a dead end
            links = zeros(5, 5)
            links[1, 2] = 0.2
            links[2, 5] = 0.5
            links[3, 4] = 0.2
            s_links = ADRIA.out_spread(links)

            @test s_links[1] > s_links[3] > 0.0
        end

        @testset "Large matrices use the sparse eigen-solver" begin
            rng = Xoshiro(7)
            n = 300
            big = sprand(rng, n, n, 0.05)
            @test ADRIA.out_spread(big) ≈ reference(big) rtol = 1e-5
            @test ADRIA.out_spread(Matrix(big)) ≈ ADRIA.out_spread(big)
        end

        @testset "Attenuation fraction" begin
            @test ADRIA.out_spread(C; fraction=0.5) ≈ s
            @test ADRIA.out_spread(C; fraction=0.9) ≈ reference(C; fraction=0.9)
            @test ADRIA.out_spread(C; fraction=0.9) != s
        end
    end

    @testset "Domain" begin
        dom = ADRIA_DOM_45
        scores = dom.loc_data.out_spread

        @test length(scores) == length(dom.loc_ids)
        @test all(0.0 .<= scores .<= 1.0)
        @test scores ≈ ADRIA.out_spread(dom.conn.data)
        @test ADRIA.location_out_spread(dom) === scores

        # Domains without the stored column compute it
        no_column = deepcopy(dom)
        DataFrames.select!(no_column.loc_data, DataFrames.Not(:out_spread))
        @test !hasproperty(no_column.loc_data, :out_spread)
        @test ADRIA.location_out_spread(no_column) ≈ scores
    end

    @testset "Decision preferences" begin
        dom = ADRIA_DOM_45
        scens = ADRIA.sample_guided(dom, 2^2)
        row = scens[1, :]
        cube = ADRIA.DataCube(
            Matrix(scens); scenarios=1:DataFrames.nrow(scens), factors=names(scens)
        )

        for W in (ADRIA.CAqCriteriaWeights, ADRIA.FogCriteriaWeights, ADRIA.LvMCriteriaWeights)
            expected_names = ADRIA.decision.criterion_names(W)
            expected_weights = Float64[row[f] for f in ADRIA.decision.criterion_fields(W)]

            from_row = ADRIA.decision.DecisionPreferences(W, dom, row)
            from_cube = ADRIA.decision.DecisionPreferences(W, dom, cube[1, :])
            defaults = ADRIA.decision.DecisionPreferences(W, dom)

            @test from_row.names == expected_names
            @test from_row.weights ≈ expected_weights
            @test from_cube.weights ≈ expected_weights
            @test defaults.names == expected_names
            @test length(defaults.weights) == length(expected_names)
            @test all(d -> d === maximum || d === minimum, from_row.directions)
        end

        prefs = ADRIA.decision.DecisionPreferences(ADRIA.CAqCriteriaWeights, dom, row)
        @test prefs.directions[findfirst(==(:out_spread), prefs.names)] === maximum
    end

    @testset "Guided setup fills out_spread in the decision matrix" begin
        dom = ADRIA_DOM_45
        scens = ADRIA.sample_selection(dom, 2^2)
        cube = ADRIA.DataCube(
            Matrix(scens); scenarios=1:DataFrames.nrow(scens), factors=names(scens)
        )
        params = cube[1, :]

        # Keep every other location, to check the rows of out_spread follow the valid locations
        depth_ok = BitVector(isodd.(1:length(dom.loc_ids)))

        for W in (ADRIA.CAqCriteriaWeights, ADRIA.FogCriteriaWeights, ADRIA.LvMCriteriaWeights)
            pref, dm, _ = ADRIA.setup_guided_intervention(
                dom, params, depth_ok, W, String.(dom.loc_ids), false, identity
            )

            valid = (dom.loc_data.k .> 0.0) .& depth_ok
            @test size(dm.data, 1) == count(valid)
            @test :out_spread in pref.names
            @test dm.data[:, findfirst(==(:out_spread), pref.names)] ==
                  ADRIA.location_out_spread(dom)[valid]
        end
    end
end
