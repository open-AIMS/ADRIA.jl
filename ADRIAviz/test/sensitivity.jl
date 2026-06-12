using ADRIA: DataCube
using ADRIAviz
using Random

@testset "sensitivity" begin
    @testset "pawn(Si) returns Figure" begin
        # Fixture matches the structure returned by ADRIA.sensitivity.pawn():
        # dimensions are :factors (Symbol) × :Si (Symbol col names)
        factor_names = Symbol.(["Factor1", "Factor2", "Factor3"])
        Si_cols = [:min, :lb, :mean, :median, :ub, :max, :std, :cv]
        Si = DataCube(
            rand(length(factor_names), length(Si_cols));
            factors=factor_names, Si=Si_cols
        )

        @test ADRIA.viz.pawn(Si) isa Figure
    end
end
