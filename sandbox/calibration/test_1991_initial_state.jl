using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ADRIA, COTSMod, CSV, DataFrames, Test

@testset "Opt-in site-keyed COTS initial states" begin
    mktempdir() do directory
        path = joinpath(directory, "initial.csv")
        table = DataFrame(site_id=["b", "a"], juvenile_ha=[0.0, 1.0],
            subadult_ha=[0.0, 2.0], adult_ha=[0.5, 3.0])
        CSV.write(path, table)
        state = ADRIA.load_cots_initial_state_csv(path, ["a", "b"])
        @test state == [1.0 2.0 3.0; 0.0 0.0 0.5]
        models = COTSMod.initialize_cots(2, COTSMod.COTSParams();
            spatial_initial_density=[0.0, 0.0])
        ADRIA.apply_cots_initial_state!(models, state)
        @test collect(models[1].N) == [1.0, 2.0, 3.0]
        @test collect(models[2].N) == [0.0, 0.0, 0.5]
        @test_throws ErrorException ADRIA.load_cots_initial_state_csv(path, ["a", "c"])
        table.site_id[2] = "b"
        CSV.write(path, table)
        @test_throws ErrorException ADRIA.load_cots_initial_state_csv(path, ["a", "b"])
        table.site_id[2] = "a"
        table.adult_ha[2] = -1.0
        CSV.write(path, table)
        @test_throws ErrorException ADRIA.load_cots_initial_state_csv(path, ["a", "b"])
    end
end

@testset "Opt-in diagnostic grazing bypass" begin
    params = COTSMod.COTSParams()
    coral = fill(0.08, 5, 1, 2)
    prey = ADRIA.CotsPreyMap([1, 2, 3], [4, 5])
    models = COTSMod.initialize_cots(2, params;
        spatial_initial_density=[1.0, 1.0])
    ADRIA.apply_cots_initial_state!(models, [1.0 1.0 10.0; 1.0 1.0 10.0])
    direct_models = deepcopy(models)
    default_models = deepcopy(models)
    bypass_models = deepcopy(models)
    direct_coral = copy(coral)
    default_coral = copy(coral)
    bypass_coral = copy(coral)
    ADRIA.apply_predation!(direct_coral, direct_models, prey)
    ADRIA.apply_cots_predation_diagnostic!(default_coral, default_models, prey,
        2005, nothing)
    ADRIA.apply_cots_predation_diagnostic!(bypass_coral, bypass_models, prey,
        2005, 2004:2012)
    @test default_coral == direct_coral
    @test bypass_coral == coral
    @test any(direct_coral .< coral)
    @test all(default_models[i].N == direct_models[i].N for i in eachindex(models))
    @test all(bypass_models[i].N == direct_models[i].N for i in eachindex(models))
end
