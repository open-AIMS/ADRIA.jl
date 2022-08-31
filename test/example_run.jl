using Test
using ADRIA


@testset "Full example run" begin
    orig_loc = pwd()
    this_loc = @__DIR__
    cd(this_loc)

    # Run full example to make sure nothing errors
    ADRIA.setup()  # Load and apply configuration options

    # Use a temporary directory for result location
    ENV["ADRIA_OUTPUT_DIR"] = mktempdir()

    ex_domain = ADRIA.load_domain("../examples/Example_domain", "45")
    p_df = ADRIA.load_scenarios(ex_domain, "../examples/example_scenarios.csv")
    ex_domain = ADRIA.run_scenarios(p_df, ex_domain)
    rs = ADRIA.load_results(ex_domain)

    cd(orig_loc)
end


