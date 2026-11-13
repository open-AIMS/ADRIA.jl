using Test
using ADRIA
using ADRIA: StaticArrays, SparseArrays, Random

# Access internal types via ADRIA module
const CotsHuman = ADRIA.CotsHuman
const CotsPreyMap = ADRIA.CotsPreyMap
const MVector = StaticArrays.MVector

# Default COTS parameter set for testing
function test_cots_params()
    return ADRIA.COTSMod.COTSParams(;
        a=1.5,
        b=0.5,
        IMM=0.002,
        p_tilde=0.97,
        C_max=0.8,
        m1=0.4,
        m2=0.2,
        m3=0.1,
        a_F=0.6,
        a_S=0.15,
        h=0.0,
        eta_F=1.0,
        eta_S=1.0,
        eta_starve=2.0,
        eta_imm=2.0,
        imm_threshold=0.35,
        fecundity_gate=false,
        a_ricker=6.0,
        b_ricker=0.1,
        tau_condition=5.0,
        allee_threshold=1.0
    )
end

@testset "COTS Submodel" begin
    @testset "CotsParams factors" begin
        defaults = ADRIA.CotsParams()
        @test defaults.cots_on.val == 0
        for mech in (:adult_senescence, :allee_stretched, :juvenile_predation,
            :per_capita_consumption, :preferred_condition, :preferred_prey_starvation)
            @test getfield(defaults, mech).val == 0
        end
        @test ADRIA.COTS_FACTOR_NAMES == Set(string.(fieldnames(ADRIA.CotsParams)))

        inst = ADRIA.create_cots_params_instance(;
            overrides=Dict("cots_on" => 1.0, "allee_threshold" => 3.0)
        )
        @test inst.cots_on.val == 1
        @test inst.allee_threshold.val == 3.0
        # Constant at the value by default; the full bounds are kept for sampling studies
        @test inst.allee_threshold.dist_params == [3.0, 3.0, 3.0]
        @test inst.allee_threshold.default_dist_params ==
            defaults.allee_threshold.default_dist_params
        @test defaults.cots_on.dist_params == [0.0, 0.0]
    end

    @testset "calibrated COTS values load by label" begin
        YA = ADRIA.YAXArrays
        ds = YA.Dataset(;
            cots=YA.YAXArray(
                (
                    ADRIA.DimensionalData.Dim{:cots_param}([
                        "juvenile_predation_capacity", "cots_on"
                    ]),
                ),
                [25.0, 1.0]
            )
        )
        overrides = ADRIA._cots_calib_overrides(ds)
        @test overrides["juvenile_predation_capacity"] == 25.0
        @test overrides["cots_on"] == 1.0
        @test isempty(
            ADRIA._cots_calib_overrides(
                YA.Dataset(;
                    other=YA.YAXArray(
                        (ADRIA.DimensionalData.Dim{:x}(1:2),), [1.0, 2.0])
                )
            )
        )
    end

    @testset "factors map onto COTSMod parameters" begin
        names = collect(string.(fieldnames(ADRIA.CotsParams)))
        vals = [ADRIA.COTS_FACTOR_SPECS[Symbol(n)].val for n in names]
        ps = ADRIA.DataCube(Float64.(vals); factors=names)
        p0 = ADRIA.cots_runtime_params(ps)
        @test !p0.per_capita_consumption && !p0.juvenile_predation
        @test p0.senescence_age isa Int && p0.senescence_age == 6

        vals[findfirst(==("per_capita_consumption"), names)] = 1.0
        vals[findfirst(==("juvenile_predation_capacity"), names)] = 25.0
        p1 = ADRIA.cots_runtime_params(ADRIA.DataCube(Float64.(vals); factors=names))
        @test p1.per_capita_consumption
        @test p1.juvenile_predation_capacity == 25.0
    end

    @testset "derived RNG seed ignores COTS factors" begin
        YA = ADRIA.YAXArrays
        base_names = ["guided", "eff_dhw_base", "RCP"]
        base_vals = [1.0, 0.63, 45.0]
        cots_names = sort(collect(ADRIA.COTS_FACTOR_NAMES))
        base = ADRIA.DataCube(base_vals; factors=base_names)
        with_cots = ADRIA.DataCube(vcat(base_vals, fill(7.5, length(cots_names)));
            factors=vcat(base_names, cots_names))
        withenv("ADRIA_RNG_SEED" => nothing) do
            @test rand(ADRIA.set_random_seed(base), 4) ==
                rand(ADRIA.set_random_seed(with_cots), 4)
        end
    end

    @testset "CotsHuman construction" begin
        p = test_cots_params()
        N = MVector{3,Float64}(0.0, 0.0, 0.5)
        model = CotsHuman(N, 0.8, p)

        @test model.N[1] == 0.0
        @test model.N[2] == 0.0
        @test model.N[3] == 0.5
        @test model.body_condition == 0.8
        @test model.params.p_tilde == 0.97
    end

    @testset "cots_timestep! basic dynamics" begin
        p = test_cots_params()
        N = MVector{3,Float64}(0.1, 0.2, 1.0)
        model = CotsHuman(N, 0.8, p)

        # Run with moderate coral cover
        Cons_F, Cons_S = ADRIA.cots_timestep!(model, 0.3, 0.1)

        # Population should remain non-negative and finite
        @test all(model.N .>= 0.0)
        @test all(isfinite.(model.N))

        # Consumption should be non-negative and bounded by available cover
        @test Cons_F >= 0.0
        @test Cons_S >= 0.0
        @test Cons_F <= 0.3
        @test Cons_S <= 0.1
    end

    @testset "cots_timestep! starvation with zero coral" begin
        p = test_cots_params()
        N = MVector{3,Float64}(0.0, 0.5, 2.0)
        model = CotsHuman(N, 0.8, p)

        # Run with zero coral — starvation should kick in
        ADRIA.cots_timestep!(model, 0.0, 0.0)

        # Adults should survive but be severely reduced
        @test model.N[3] < 2.0
        @test all(model.N .>= 0.0)
    end

    @testset "Allee effect suppresses recruitment at low density" begin
        p = test_cots_params()

        # Low adult density — below Allee threshold of 1.0
        N_low = MVector{3,Float64}(0.0, 0.0, 0.1)
        model_low = CotsHuman(N_low, 0.8, p)
        ADRIA.cots_timestep!(model_low, 0.5, 0.2)
        recruits_low = model_low.N[1]

        # High adult density — above Allee threshold
        N_high = MVector{3,Float64}(0.0, 0.0, 3.0)
        model_high = CotsHuman(N_high, 0.8, p)
        ADRIA.cots_timestep!(model_high, 0.5, 0.2)
        recruits_high = model_high.N[1]

        # Recruitment should be much higher when adults are above threshold
        @test recruits_high > recruits_low * 5.0
    end

    @testset "cots_mortality! reduces cover" begin
        n_groups = 5
        n_sizes = 7
        n_locs = 3
        p = test_cots_params()

        C_cover_t = fill(0.05, n_groups, n_sizes, n_locs)
        prey_map = CotsPreyMap([1, 2, 3], [4, 5])

        # Create models with active adult COTS
        models = [
            CotsHuman(MVector{3,Float64}(0.0, 0.0, 2.0), 0.8, p) for _ in 1:n_locs
        ]

        cover_before = sum(C_cover_t)
        ADRIA.cots_mortality!(C_cover_t, models, prey_map)
        cover_after = sum(C_cover_t)

        # Cover should decrease due to predation
        @test cover_after < cover_before

        # All cover values should remain in [0, 1]
        @test all(0.0 .<= C_cover_t .<= 1.0)
    end

    @testset "disperse_cots_larvae!" begin
        p = test_cots_params()
        n_locs = 4

        models = [
            CotsHuman(MVector{3,Float64}(10.0, 0.0, 0.0), 0.8, p),  # Source with recruits
            CotsHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p),
            CotsHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p),
            CotsHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p)
        ]

        # Simple connectivity: loc 1 sends to locs 2 and 3
        conn = SparseArrays.sparse(
            [1, 1, 2, 3],    # rows (sources)
            [2, 3, 2, 3],    # cols (sinks)
            [0.5, 0.3, 0.1, 0.1], # weights
            n_locs, n_locs   # dims
        )

        ADRIA.disperse_cots_larvae!(models, conn)

        # Loc 2 and 3 should have received recruits
        @test models[2].N[1] > 0.0
        @test models[3].N[1] > 0.0
        # All values should be finite and non-negative
        @test all(m -> all(isfinite.(m.N)) && all(m.N .>= 0.0), models)
    end

    @testset "init_cots_populations" begin
        p = test_cots_params()
        n_locs = 20

        # Test with explicit seed locations
        seed = Set([1, 5, 10])
        models = ADRIA.init_cots_populations(
            n_locs, p; seed_locs=seed, init_density=2.0
        )

        @test length(models) == n_locs

        # Seeded locations should have adults
        @test models[1].N[3] == 2.0
        @test models[5].N[3] == 2.0
        @test models[10].N[3] == 2.0

        # Non-seeded locations should be empty
        @test models[2].N[3] == 0.0
        @test models[3].N[3] == 0.0

        # Test with default random seeding (no explicit seed_locs)
        models_rand = ADRIA.init_cots_populations(
            n_locs, p; rng=Random.MersenneTwister(42)
        )
        @test length(models_rand) == n_locs
        # At least some should be seeded (25% default = 5 locations)
        n_seeded = count(m -> m.N[3] > 0.0, models_rand)
        @test n_seeded >= 1
    end

    @testset "cots_mortality! preserves non-prey groups and applies proportional mortality" begin
        p = test_cots_params()
        prey_map = CotsPreyMap([1], [2])
        C_cover_t = zeros(3, 2, 1)
        C_cover_t[1, :, 1] .= [0.2, 0.4]
        C_cover_t[2, :, 1] .= [0.1, 0.3]
        C_cover_t[3, :, 1] .= [0.7, 0.9]
        non_prey_before = copy(C_cover_t[3, :, 1])

        models = [CotsHuman(MVector{3,Float64}(0.0, 0.0, 1.5), 0.8, p)]
        ADRIA.cots_mortality!(C_cover_t, models, prey_map)

        @test C_cover_t[1, 1, 1] / 0.2 ≈ C_cover_t[1, 2, 1] / 0.4
        @test C_cover_t[2, 1, 1] / 0.1 ≈ C_cover_t[2, 2, 1] / 0.3
        @test C_cover_t[3, :, 1] == non_prey_before
        @test sum(C_cover_t[1:2, :, 1]) < 1.0
        @test all(0.0 .<= C_cover_t .<= 1.0)
    end

    @testset "disperse_cots_larvae! skips self-connectivity and applies scalar" begin
        p = test_cots_params()
        n_locs = 4
        models = [
            CotsHuman(MVector{3,Float64}(10.0, 0.0, 0.0), 0.8, p),
            CotsHuman(MVector{3,Float64}(1.0, 0.0, 0.0), 0.8, p),
            CotsHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p),
            CotsHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, p)
        ]

        conn = SparseArrays.sparse(
            [1, 1, 1, 2, 2],
            [1, 2, 3, 2, 4],
            [0.9, 0.5, 0.3, 0.1, 0.25],
            n_locs, n_locs
        )

        ADRIA.disperse_cots_larvae!(models, conn; immigration_scalar=2.0)

        @test models[1].N[1] == 10.0
        @test models[2].N[1] == 11.0
        @test models[3].N[1] == 6.0
        @test models[4].N[1] == 0.5
    end

    @testset "package-shaped API wrappers" begin
        p = test_cots_params()
        seed = Set([1, 3])
        cots_state = ADRIA.initialize_cots(4, p; seed_locs=seed, init_density=1.2)

        @test cots_state isa Vector{CotsHuman}
        @test cots_state[1].N[3] == 1.2
        @test cots_state[2].N[3] == 0.0
        @test cots_state[3].N[3] == 1.2

        disabled_state = ADRIA.initialize_cots(3, p; enabled=false)
        @test all(m -> all(m.N .== 0.0), disabled_state)

        C_cover_t = fill(0.05, 5, 2, 4)
        cover_before = sum(C_cover_t)
        ADRIA.apply_predation!(C_cover_t, cots_state, CotsPreyMap([1, 2, 3], [4, 5]))
        @test sum(C_cover_t) < cover_before

        conn = SparseArrays.sparse([1], [2], [0.5], 4, 4)
        before_recruits = cots_state[2].N[1]
        ADRIA.disperse_larvae!(cots_state, conn; scalar=2.0)
        @test cots_state[2].N[1] >= before_recruits
    end
end
