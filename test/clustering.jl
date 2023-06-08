
@testset "Temporal clustering" begin
    d1 = [1.; 2.; 3.]
    d2 = [10.; 20.; 30.]
    d3 = [1.; 5.; 8.]

    test_data::Matrix = [d1 d2 d3]
    
    @testset "Compute CE (Complexity)" begin
        # Compute CD for test_data
        ce = ADRIA.analysis.complexity(test_data)

        # CE is a N Vector, where N is the number of rows in test_data
        @test ce isa Vector
        @test length(ce) == size(test_data, 2)

        # Expected results
        @test ce[1] == sqrt(2.)
        @test ce[2] == sqrt(200.)
        @test ce[3] == sqrt(25.)
    end

    @testset "Compute CF (Correlation Factor)" begin
        # mock ce vector
        ce = [2.5, 207.0, 25.0, 25.0]

        # Expected Results
        @test ADRIA.analysis.correlation_factor(ce[1], ce[2]) == 207.0 / 2.5
        @test ADRIA.analysis.correlation_factor(ce[2], ce[3]) == 207.0 / 25.0
        @test ADRIA.analysis.correlation_factor(ce[1], ce[3]) == 25.0 / 2.5
        @test ADRIA.analysis.correlation_factor(ce[3], ce[4]) == 1
    end

    @testset "Comput CID Matrix (Complexity Invariance Matrix)" begin
        cid = ADRIA.analysis.complexity_invariance_distance(test_data)

        # CID is a Matrix (N,N)
        @test size(cid, 1) == size(cid, 2) == size(test_data, 2)

        # All CID are positive
        @testset "CID positivity" for i in cid
            @test i >= 0
        end

        # CID ij and ji entries are the same
        @test cid[1,2] == cid[2,1] >= 0
        @test cid[1,3] == cid[3,1] >= 0
        @test cid[2,3] == cid[3,2] >= 0

        # CID (i,i) is null
        @test cid[1,1] == cid[2,2] == cid[3,3] == 0.0
    end
end
