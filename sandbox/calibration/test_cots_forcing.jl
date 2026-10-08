using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ADRIA
using Test

@testset "COTS annual forcing selection and external boundary" begin
    outbreak_potential = zeros(2, 2, 2)
    outbreak_pelagic = zeros(2, 2, 2)
    outbreak_potential[:, 1, 1] .= [0.5, 1.0]
    outbreak_pelagic[:, 1, 1] .= [0.1, 0.25]
    forcing = ADRIA.COTSConnectivityForcing(
        [2010, 2011],
        [[0.2 0.3; 0.4 0.5], [0.1 0.6; 0.2 0.7]],
        [1, 1, 2],
        [2.0, 1.0],
        [0.5 0.6; 0.7 0.8],
        [0.1 0.2; 0.3 0.4],
        [1991, 1992],
        outbreak_potential,
        outbreak_pelagic
    )
    @test ADRIA.cots_forcing_index(forcing, 1, "cycle", 42) == 1
    @test ADRIA.cots_forcing_index(forcing, 2, "cycle", 42) == 2
    @test ADRIA.cots_forcing_index(forcing, 3, "cycle", 42) == 1
    @test ADRIA.cots_forcing_index(forcing, 4, "sample", 42) ==
        ADRIA.cots_forcing_index(forcing, 4, "sample", 42)
    @test ADRIA.cots_external_supply_by_site(forcing, 1, 2.0) == [0.1, 0.1, 0.6]
    potential, pelagic = ADRIA.cots_external_outbreak_supply_by_site(
        forcing, 1, 1991, 2.0
    )
    @test potential == [1.0, 1.0, 2.0]
    @test pelagic == [0.2, 0.2, 0.5]
    mean_potential, mean_pelagic = ADRIA.cots_external_outbreak_supply_by_site(
        forcing, 0, 1991, 2.0
    )
    @test mean_potential == [0.5, 0.5, 1.0]
    @test mean_pelagic == [0.1, 0.1, 0.25]
    missing_potential, missing_pelagic = ADRIA.cots_external_outbreak_supply_by_site(
        forcing, 1, 1985, 2.0
    )
    @test all(iszero, missing_potential)
    @test all(iszero, missing_pelagic)
    @test_throws ErrorException ADRIA.cots_forcing_index(forcing, 1, "invalid", 42)
end

@testset "Opt-in low-cover mortality adapter" begin
    base = (a=1.5, b=0.5, IMM=0.0, p_tilde=0.95, C_max=0.5,
        m1=0.4, m2=0.2, m3=0.1, a_F=0.6, a_S=0.15)
    legacy = ADRIA._cots_runtime_params(base)
    @test !legacy.low_cover_mortality
    @test legacy.low_cover_threshold == 0.10
    @test legacy.low_cover_strength == 0.0
    treatment = ADRIA._cots_runtime_params(merge(base, (
        low_cover_mortality=true, low_cover_threshold=0.10,
        low_cover_strength=0.5)))
    @test treatment.low_cover_mortality
    @test treatment.low_cover_threshold == 0.10
    @test treatment.low_cover_strength == 0.5
end

@testset "Opt-in internal settlement blackout" begin
    @test ADRIA.cots_internal_settlement_scalar(1.0, 2005, nothing) == 1.0
    @test ADRIA.cots_internal_settlement_scalar(0.4, 1999, 2000:2010) == 0.4
    @test ADRIA.cots_internal_settlement_scalar(0.4, 2000, 2000:2010) == 0.0
    @test ADRIA.cots_internal_settlement_scalar(0.4, 2010, 2000:2010) == 0.0
    @test ADRIA.cots_internal_settlement_scalar(0.4, 2011, 2000:2010) == 0.4
end
