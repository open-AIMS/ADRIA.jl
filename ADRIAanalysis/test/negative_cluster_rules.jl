using Test, ADRIAanalysis

@testset "cluster_rules/rules without SIRUS: stub methods only" begin
    # Assert method count = 0 for typed methods (only the varargs stub exists).
    # Using methods() count rather than call-and-catch avoids false positives
    # from argument-type mismatches.
    @test length(methods(ADRIAanalysis.cluster_rules)) == 1
    @test length(methods(ADRIAanalysis.rules)) == 1
    # Confirm stub raises informative error, not MethodError.
    @test_throws ErrorException ADRIAanalysis.cluster_rules()
    @test_throws ErrorException ADRIAanalysis.rules()
end
