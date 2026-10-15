@testset "scenario_groups" begin
    df = DataFrames.DataFrame(;
        iv_CAq_N_TA=Int[0, 100, 100],
        iv_Fog=Float64[0.0, 0.0, 0.0],
        iv_Shd=Float64[0.0, 0.0, 0.0],
        iv_LvM_N_settlers=Int[0, 0, 0],
        iv_MCB_duration=Float64[0.0, 0.0, 0.0],
        guided=Int[0, 0, 1],
        RCP=Int[45, 45, 60]
    )

    @testset "scenario_types" begin
        result = ADRIA.analysis.scenario_types(df)

        @test result isa Dict{Symbol,BitVector}
        @test haskey(result, :counterfactual)
        @test haskey(result, :unguided)
        @test haskey(result, :guided)

        @test result[:counterfactual] == BitVector([true, false, false])
        @test result[:unguided] == BitVector([false, true, false])
        @test result[:guided] == BitVector([false, false, true])
    end

    @testset "scenario_rcps" begin
        result = ADRIA.analysis.scenario_rcps(df)

        @test result isa Dict{Symbol,BitVector}
        @test haskey(result, :RCP45)
        @test haskey(result, :RCP60)

        @test result[:RCP45] == BitVector([true, true, false])
        @test result[:RCP60] == BitVector([false, false, true])
    end

    @testset "absent type is excluded" begin
        df2 = DataFrames.DataFrame(;
            iv_CAq_N_TA=Int[0, 100],
            iv_Fog=Float64[0.0, 0.0],
            iv_Shd=Float64[0.0, 0.0],
            iv_LvM_N_settlers=Int[0, 0],
            iv_MCB_duration=Float64[0.0, 0.0],
            guided=Int[0, 0],
            RCP=Int[45, 45]
        )
        result = ADRIA.analysis.scenario_types(df2)
        @test !haskey(result, :guided)
    end
end
