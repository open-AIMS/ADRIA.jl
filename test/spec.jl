@testset "Extracting component parameters" begin
    orig_loc = pwd()
    this_loc = @__DIR__
    cd(this_loc)

    # Run full example to make sure nothing errors
    ADRIA.setup()  # Load and apply configuration options

    # Use a temporary directory for result location
    ENV["ADRIA_OUTPUT_DIR"] = mktempdir()

    dom = ADRIA.load_domain("../examples/Example_domain", "45")

    x = ADRIA.component_params(dom.model, Intervention)
    @test size(x, 1) > 0

    x = ADRIA.component_params(dom.model, ADRIA.Coral)
    @test size(x, 1) > 0

    x = ADRIA.component_params(dom.model, Criteria)
    @test size(x, 1) > 0

end