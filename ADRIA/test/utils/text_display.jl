using ADRIA: human_readable_name, get_scientific_factors

@testset "text_display" begin
    @testset "human_readable_name" begin
        @testset "Symbol inputs" begin
            input_vector_symbols = [:time_steps, :locations]
            input_single_symbol = :time_steps

            @test human_readable_name(input_vector_symbols) ==
                ["time steps", "locations"]
            @test human_readable_name(input_vector_symbols; title_case=true) ==
                ["Time Steps", "Locations"]
            @test human_readable_name(input_single_symbol) == "time steps"
            @test human_readable_name(input_single_symbol; title_case=true) ==
                "Time Steps"
        end

        @testset "Strings inputs" begin
            input_vector_strings = ["coral_cover", "1234_test"]
            input_string = "coral_cover"

            @test human_readable_name(
                input_vector_strings
            ) == ["coral cover", "1234 test"]
            @test human_readable_name(input_vector_strings; title_case=true) ==
                ["Coral Cover", "1234 Test"]
            @test human_readable_name(input_string) == "coral cover"
            @test human_readable_name(input_string; title_case=true) ==
                "Coral Cover"
        end
    end

    @testset "get_scientific_factors" begin
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

        @testset "0 < x < 1" begin
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

        @testset "0 > x > -1" begin
            @test get_scientific_factors(-0.01) == (-1.00, -2)
            @test get_scientific_factors(-0.01000001) == (-1.00, -2)
            @test get_scientific_factors(-0.000123) == (-1.23, -4)
            @test get_scientific_factors(-0.000123456) == (-1.23, -4)
            @test get_scientific_factors(-0.000123456; digits=4) == (-1.2345, -4)
            @test get_scientific_factors(-0.0093879; digits=3) == (-9.387, -3)
        end
    end
end
