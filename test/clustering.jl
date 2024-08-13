
@testset "Temporal clustering" begin
    @testset "Variable series" begin
        d1 = [1.0; 2.0; 3.0]
        d2 = [10.0; 20.0; 30.0]
        d3 = [1.0; 5.0; 8.0]

        test_data::Matrix = [d1 d2 d3]

        @testset "Compute CE (Complexity)" begin
            # Compute CD for test_data
            ce = ADRIA.analysis._complexity(test_data)

            # CE is a N Vector, where N is the number of rows in test_data
            @test ce isa Vector
            @test length(ce) == size(test_data, 2)

            # Expected results
            @test ce[1] == sqrt(2.0) + 1
            @test ce[2] == sqrt(200.0) + 1
            @test ce[3] == sqrt(25.0) + 1
        end

        @testset "Compute CF (Correction Factor)" begin
            # mock ce vector
            ce = [2.5, 207.0, 25.0, 25.0]

            # Expected Results
            @test ADRIA.analysis.correction_factor(ce[1], ce[2]) == 207.0 / 2.5
            @test ADRIA.analysis.correction_factor(ce[2], ce[3]) == 207.0 / 25.0
            @test ADRIA.analysis.correction_factor(ce[1], ce[3]) == 25.0 / 2.5
            @test ADRIA.analysis.correction_factor(ce[3], ce[4]) == 1
        end

        @testset "Comput CID Matrix (Complexity Invariance Matrix)" begin
            complexity = ADRIA.analysis._complexity(test_data)
            cid = ADRIA.analysis.complexity_invariance_distance(test_data)

            # CID is a Matrix (N,N)
            @test size(cid, 1) == size(cid, 2) == size(test_data, 2)

            # All CID are positive
            @testset "CID positivity" for i in cid
                @test i >= 0
            end

            # CID ij and ji entries are the same
            @test cid[1, 2] == cid[2, 1] >= 0
            @test cid[1, 3] == cid[3, 1] >= 0
            @test cid[2, 3] == cid[3, 2] >= 0

            # CID (i,i) is null
            @test cid[1, 1] == cid[2, 2] == cid[3, 3] == 0.0
        end

        @testset "Call cluster_series function" begin
            # Since test_data is a 3x3 matrix, Clustering.kmeioids requires num_clusters ≤ 3
            num_clusters = 3
            clusters = ADRIA.analysis.cluster_series(test_data, num_clusters)

            @test length(clusters) == size(test_data, 2)
            @test -1 ∉ clusters
            @test 0 ∉ clusters
            @test 1 ∈ clusters
        end
    end

    @testset "Data with some constant series" begin
        d1 = [1.0; 2.0; 3.0]
        d2 = [2.0; 3.0; 4.0]
        d3 = [6.0; 3.0; 1.0]
        d4 = [0.0; 0.0; 0.0]
        d5 = [1.0; 1.0; 1.0]

        const_test_data::Matrix = [d1 d2 d3 d4 d5]

        @testset "Compute CE (Complexity)" begin
            # Compute CD for test_data
            ce = ADRIA.analysis._complexity(const_test_data)

            # CE is a N Vector, where N is the number of rows in test_data
            @test ce isa Vector
            @test length(ce) == size(const_test_data, 2)

            # Expected results
            @test ce[1] == sqrt(2.0) + 1
            @test ce[2] == sqrt(2.0) + 1
            @test ce[3] == sqrt(13.0) + 1
            @test ce[4] == 1
            @test ce[5] == 1
        end

        @testset "Call cluster_series function" begin
            num_clusters = 5
            clusters = ADRIA.analysis.cluster_series(const_test_data, num_clusters)

            @test length(clusters) == size(const_test_data, 2)
        end
    end
end
