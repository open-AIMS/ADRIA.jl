using ADRIA: DataCube
using ADRIAviz
using Random

@testset "sensitivity" begin
    @testset "pawn(Si) returns Figure" begin
        # Create synthetic sensitivity indices (factors × stats)
        factors = ["Factor1", "Factor2", "Factor3"]
        stats = [:median, :mean, :std]
        Si = DataCube(
            rand(length(factors), length(stats));
            factors=factors, stats=stats
        )

        @test ADRIA.viz.pawn(Si) isa Figure
    end
end
