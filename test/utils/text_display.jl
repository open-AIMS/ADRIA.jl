@testset "text_display" begin
    @testset "get_scientific_factors" begin
        get_scientific_factors = ADRIA.get_scientific_factors

        @testset "Edgecases" begin
            @test get_scientific_factors(0.0) == (0.0, 0)
            @test get_scientific_factors(1.0) == (1.0, 0)
            @test get_scientific_factors(-1.0) == (-1.0, 0)
            @test get_scientific_factors(1.01) == (1.01, 0)
            @test get_scientific_factors(-1.01) == (-1.01, 0)
            @test get_scientific_factors(10.0) == (1.00, 1)
            @test get_scientific_factors(100.0) == (1.00, 2)
            @test get_scientific_factors(10000.0) == (1.00, 4)
            @test get_scientific_factors(-10.0) == (-1.00, 1)
            @test get_scientific_factors(-1000.0) == (-1.00, 3)
        end

        @testset "Integers" begin
            @test get_scientific_factors(123) == get_scientific_factors(123.0)
        end

        @testset "x > 10" begin
            @test get_scientific_factors(123.0) == (1.23, 2)
            @test get_scientific_factors(12.0) == (1.2, 1)
            @test get_scientific_factors(10001.0; digits=3) == (1.000, 4)
        end

        @testset "0 > x > 1" begin
            @test get_scientific_factors(0.01) == (1.00, -2)
            @test get_scientific_factors(0.01000001) == (1.00, -2)
            @test get_scientific_factors(0.000123) == (1.23, -4)
            @test get_scientific_factors(0.000123456) == (1.23, -4)
            @test get_scientific_factors(0.000123456; digits=4) == (1.2345, -4)
            @test get_scientific_factors(0.0093879; digits=3) == (9.387, -3)
        end

        @testset "x < -10" begin
            @test get_scientific_factors(-123.0) == (-1.23, 2)
            @test get_scientific_factors(-12.0) == (-1.2, 1)
        end
    end
end
