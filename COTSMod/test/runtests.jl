using Test
using Random
using SparseArrays
using StaticArrays
using COTSMod

const MVector = StaticArrays.MVector

test_cots_params() = COTSParams(
    a = 1.5,
    b = 0.5,
    IMM = 0.002,
    p_tilde = 0.97,
    C_max = 0.8,
    m1 = 0.4,
    m2 = 0.2,
    m3 = 0.1,
    a_F = 0.6,
    a_S = 0.15,
    h = 0.0,
    eta_F = 1.0,
    eta_S = 1.0,
    eta_starve = 2.0,
    eta_imm = 2.0,
    imm_threshold = 0.35,
    fecundity_gate = false,
    a_ricker = 6.0,
    b_ricker = 0.1,
    tau_condition = 5.0,
    allee_threshold = 1.0
)

@testset "COTSMod" begin
    @testset "COTSHuman construction" begin
        p = test_cots_params()
        model = COTSHuman(MVector{3,Float64}(0.0, 0.0, 0.5), 0.8, p)
        @test model.N[1] == 0.0
        @test model.N[2] == 0.0
        @test model.N[3] == 0.5
        @test model.body_condition == 0.8
        @test model.params.p_tilde == 0.97
    end

    @testset "cots_timestep! basic dynamics" begin
        p = test_cots_params()
        model = COTSHuman(MVector{3,Float64}(0.1, 0.2, 1.0), 0.8, p)
        Cons_F, Cons_S = cots_timestep!(model, 0.3, 0.1)
        @test all(model.N .>= 0.0)
        @test all(isfinite.(model.N))
        @test 0.0 <= Cons_F <= 0.3
        @test 0.0 <= Cons_S <= 0.1
    end

    @testset "cots_timestep! starvation with zero coral" begin
        p = test_cots_params()
        model = COTSHuman(MVector{3,Float64}(0.0, 0.5, 2.0), 0.8, p)
        cots_timestep!(model, 0.0, 0.0)
        @test model.N[3] < 2.0
        @test all(model.N .>= 0.0)
    end

    @testset "Allee effect suppresses recruitment at low density" begin
        p = test_cots_params()
        model_low = COTSHuman(MVector{3,Float64}(0.0, 0.0, 0.1), 0.8, p)
        model_high = COTSHuman(MVector{3,Float64}(0.0, 0.0, 3.0), 0.8, p)
        cots_timestep!(model_low, 0.5, 0.2)
        cots_timestep!(model_high, 0.5, 0.2)
        @test model_high.N[1] > model_low.N[1] * 5.0
    end

    @testset "cots_mortality! reduces cover" begin
        p = test_cots_params()
        C_cover_t = fill(0.05, 5, 7, 3)
        prey_map = COTSPreyMap([1, 2, 3], [4, 5])
        models = [COTSHuman(MVector{3,Float64}(0.0, 0.0, 2.0), 0.8, p) for _ in 1:3]
        cover_before = sum(C_cover_t)
        cots_mortality!(C_cover_t, models, prey_map)
        @test sum(C_cover_t) < cover_before
        @test all(0.0 .<= C_cover_t .<= 1.0)
    end

    @testset "disperse_cots_larvae!" begin
        p = test_cots_params()
        models = [
            COTSHuman(MVector{3,Float64}(10.0, 0.0, 0.0), 0.8, p),
            COTSHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p),
            COTSHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p),
            COTSHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p)
        ]
        conn = sparse([1, 1, 2, 3], [2, 3, 2, 3], [0.5, 0.3, 0.1, 0.1], 4, 4)
        disperse_cots_larvae!(models, conn)
        @test models[2].N[1] > 0.0
        @test models[3].N[1] > 0.0
        @test all(m -> all(isfinite.(m.N)) && all(m.N .>= 0.0), models)
    end

    @testset "init_cots_populations" begin
        p = test_cots_params()
        seed = Set([1, 5, 10])
        models = init_cots_populations(20, p; seed_locs=seed, init_density=2.0)
        @test length(models) == 20
        @test models[1].N[3] == 2.0
        @test models[5].N[3] == 2.0
        @test models[10].N[3] == 2.0
        @test models[2].N[3] == 0.0
        models_rand = init_cots_populations(20, p; rng=Random.MersenneTwister(42))
        @test count(m -> m.N[3] > 0.0, models_rand) >= 1
    end

    @testset "inject_upstream_pulse!" begin
        p = test_cots_params()
        models = [COTSHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p) for _ in 1:10]
        inject_upstream_pulse!(models, Set([2, 5, 8]); pulse_val=2.0)
        @test models[2].N[1] == 2.0
        @test models[5].N[1] == 2.0
        @test models[8].N[1] == 2.0
        @test models[1].N[1] == 0.0
    end

    @testset "cots_mortality! preserves non-prey groups and applies proportional mortality" begin
        p = test_cots_params()
        prey_map = COTSPreyMap([1], [2])
        C_cover_t = zeros(3, 2, 1)
        C_cover_t[1, :, 1] .= [0.2, 0.4]
        C_cover_t[2, :, 1] .= [0.1, 0.3]
        C_cover_t[3, :, 1] .= [0.7, 0.9]
        non_prey_before = copy(C_cover_t[3, :, 1])
        models = [COTSHuman(MVector{3,Float64}(0.0, 0.0, 1.5), 0.8, p)]
        cots_mortality!(C_cover_t, models, prey_map)
        @test isapprox(C_cover_t[1, 1, 1] / 0.2, C_cover_t[1, 2, 1] / 0.4)
        @test isapprox(C_cover_t[2, 1, 1] / 0.1, C_cover_t[2, 2, 1] / 0.3)
        @test C_cover_t[3, :, 1] == non_prey_before
        @test sum(C_cover_t[1:2, :, 1]) < 1.0
        @test all(0.0 .<= C_cover_t .<= 1.0)
    end

    @testset "disperse_cots_larvae! skips self-connectivity and applies scalar" begin
        p = test_cots_params()
        models = [
            COTSHuman(MVector{3,Float64}(10.0, 0.0, 0.0), 0.8, p),
            COTSHuman(MVector{3,Float64}(1.0, 0.0, 0.0), 0.8, p),
            COTSHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p),
            COTSHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p)
        ]
        conn = sparse([1, 1, 1, 2, 2], [1, 2, 3, 2, 4], [0.9, 0.5, 0.3, 0.1, 0.25], 4, 4)
        disperse_cots_larvae!(models, conn; immigration_scalar=2.0)
        @test models[1].N[1] == 10.0
        @test models[2].N[1] == 11.0
        @test models[3].N[1] == 6.0
        @test models[4].N[1] == 0.5
    end

    @testset "inject_upstream_pulse! vector pulse values" begin
        p = test_cots_params()
        models = [COTSHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p) for _ in 1:10]
        pulse_vals = collect(0.1:0.1:1.0)
        inject_upstream_pulse!(models, Set([2, 5, 8, 99]), pulse_vals)
        @test models[2].N[1] == pulse_vals[2]
        @test models[5].N[1] == pulse_vals[5]
        @test models[8].N[1] == pulse_vals[8]
        @test models[1].N[1] == 0.0
        @test models[3].N[1] == 0.0
    end

    @testset "package-shaped API wrappers" begin
        p = test_cots_params()
        cots_state = initialize_cots(4, p; seed_locs=Set([1, 3]), init_density=1.2)
        @test cots_state isa Vector{COTSHuman}
        @test cots_state[1].N[3] == 1.2
        @test cots_state[2].N[3] == 0.0
        @test cots_state[3].N[3] == 1.2
        disabled_state = initialize_cots(3, p; enabled=false)
        @test all(m -> all(m.N .== 0.0), disabled_state)
        C_cover_t = fill(0.05, 5, 2, 4)
        cover_before = sum(C_cover_t)
        apply_predation!(C_cover_t, cots_state, COTSPreyMap([1, 2, 3], [4, 5]))
        @test sum(C_cover_t) < cover_before
        conn = sparse([1], [2], [0.5], 4, 4)
        before_recruits = cots_state[2].N[1]
        disperse_larvae!(cots_state, conn; scalar=2.0)
        @test cots_state[2].N[1] >= before_recruits
        pulse_vals = [0.0, 0.2, 0.4, 0.6]
        apply_external_supply!(cots_state, Set([2, 4]), pulse_vals)
        @test cots_state[4].N[1] >= pulse_vals[4]
    end
end