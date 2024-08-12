@testset "scale" begin
    @testset "linear_scale" begin
        linear_scale = ADRIA.linear_scale

        @test linear_scale(100, :cm, :m) == round(10^0)   # 100cm = 1m
        @test linear_scale(0.1, :m, :cm) == round(10^1)  # 0.1m = 10cm

        @test linear_scale(:m, :cm) == round(10^2)
        @test linear_scale(:cm, :m) == round(10^-2; digits=2)
    end

    @testset "quadratic_scale" begin
        quadratic_scale = ADRIA.quadratic_scale

        @test quadratic_scale(100, :cm, :m) == round(10^-2; digits=2)   # 100cm^2 = 1m
        @test quadratic_scale(0.1, :m, :cm) == round(10^3)  # 0.1m = 10cm

        @test quadratic_scale(:m, :cm) == round(10^4)
        @test quadratic_scale(:cm, :m) == round(10^-4; digits=4)
    end
end
