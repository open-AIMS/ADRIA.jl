if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

@testset "Extracting component parameters" begin
    # Run full example to make sure nothing errors
    ADRIA.setup()  # Load and apply configuration options

    # Use a temporary directory for result location
    ENV["ADRIA_OUTPUT_DIR"] = mktempdir()

    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, "45")

    x = ADRIA.component_params(dom.model, Intervention)
    @test size(x, 1) > 0

    x = ADRIA.component_params(dom.model, ADRIA.Coral)
    @test size(x, 1) > 0

    x = ADRIA.component_params(dom.model, ADRIA.CriteriaWeights)
    @test size(x, 1) > 0

end
