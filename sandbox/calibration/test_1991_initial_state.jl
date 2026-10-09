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

@testset "Opt-in adult size initialization and diagnostics" begin
    params = COTSMod.COTSParams(size_weighted_fecundity=true,
        size_initial_large_fraction=0.25, separated_recruitment=true,
        allee_threshold=3.0, IMM=0.0)
    models = COTSMod.initialize_cots(2, params;
        spatial_initial_density=[0.0, 0.0])
    ADRIA.apply_cots_initial_state!(models, [1.0 2.0 4.0; 0.0 0.0 2.0])
    @test [m.large_adults for m in models] == [1.0, 0.5]
    COTSMod.cots_timestep!(models[1], 0.4, 0.1)
    flows = COTSMod.cots_flow_diagnostics(models)
    @test flows.small_adults[1] + flows.large_adults[1] ≈ models[1].N[3]
    @test flows.effective_breeders[1] < 4.0
    @test flows.small_adult_growth[1] > 0.0

    legacy = COTSMod.initialize_cots(1, COTSMod.COTSParams();
        spatial_initial_density=[0.0])
    ADRIA.apply_cots_initial_state!(legacy, [0.0 0.0 4.0])
    @test legacy[1].large_adults == 0.0
    COTSMod.cots_timestep!(legacy[1], 0.4, 0.1)
    @test legacy[1].last_effective_breeders == 4.0
end

@testset "Opt-in lagged adult hazard adapter and state reset" begin
    params = COTSMod.COTSParams(IMM=0.0, p_tilde=0.0,
        lagged_adult_hazard=true, hazard_food_coupled=true,
        hazard_max=1.5, hazard_threshold_ha=1.5, hazard_width_ha=0.5)
    models = COTSMod.initialize_cots(2, params;
        spatial_initial_density=[0.0, 0.0])
    ADRIA.apply_cots_initial_state!(models, [0.0 1.0 3.0; 0.0 1.0 0.5])
    @test all(model.burden_ha == 0.0 for model in models)
    for model in models
        COTSMod.cots_timestep!(model, 0.1, 0.0)
    end
    diagnostics = ADRIA.cots_hazard_diagnostics(models)
    @test length(diagnostics.burden_ha) == 2
    @test all(iszero, diagnostics.hazard)
    @test all(iszero, diagnostics.adult_deaths_ha)
    @test models[1].burden_ha > models[2].burden_ha
    ADRIA.apply_cots_initial_state!(models, [0.0 1.0 3.0; 0.0 1.0 0.5])
    @test all(model.burden_ha == 0.0 && model.last_hazard == 0.0 &&
        model.last_hazard_deaths_ha == 0.0 for model in models)
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
