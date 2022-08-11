@testset "metric modifications" begin
    tac_metric = @extend_metric(TAC, total_absolute_cover, [site_area(domain)])

    @test tac_metric.func == total_absolute_cover.func
    @test tac_metric.dims == total_absolute_cover.dims
end