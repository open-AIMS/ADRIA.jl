if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

@testset "Extracting component parameters" begin
    # Run full example to make sure nothing errors
    ADRIA.setup()  # Load and apply configuration options

    # Use a temporary directory for result location
    ENV["ADRIA_OUTPUT_DIR"] = mktempdir()

    x = ADRIA.component_params(ADRIA_DOM_45.model, Intervention)
    @test size(x, 1) > 0

    x = ADRIA.component_params(ADRIA_DOM_45.model, ADRIA.Coral)
    @test size(x, 1) > 0

    x = ADRIA.component_params(ADRIA_DOM_45.model, ADRIA.SeedCriteriaWeights)
    @test size(x, 1) > 0

    x = ADRIA.component_params(ADRIA_DOM_45.model, ADRIA.FogCriteriaWeights)
    @test size(x, 1) > 0
end
