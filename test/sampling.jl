using ADRIA


if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "examples", "Example_domain")
end

@testset "sample" begin
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH)
    num_samples = 32
    scens = ADRIA.sample(dom, num_samples)

    ms = ADRIA.model_spec(dom)
    @test all(values(scens[1, ms.is_constant.==true]) == values(scens[end, ms.is_constant.==true])) || "Constant params are not constant!"

    min_x = values(ms[:, :lower_bound])
    max_x = values(ms[:, :upper_bound])
    eco = findall((ms.component .== "Coral") .& (ms.is_constant .== false))
    for i in 1:num_samples
        msg = "Sampled values were not in expected bounds! Expected "
        x = values(scens[i, :])

        if scens[i, :guided] > 0
            @test all(min_x .<= x .<= max_x) || "Sampled values were not in expected bounds! $(min_x .<= x .<= max_x)"
        else
            # When guided == 0, intervention parameters are set to 0 so only check ecological values
            @test all(min_x[eco] .<= x[eco] .<= max_x[eco]) || "Sampled coral values were not in expected bounds! $(min_x[eco] .<= x[eco] .<= max_x[eco])"
        end
    end
end