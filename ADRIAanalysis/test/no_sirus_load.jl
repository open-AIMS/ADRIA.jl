# IMPORTANT: This file must run in a fresh Julia process that has NOT loaded SIRUS or MLJ.
# Do NOT include this from a process that has already triggered the extension.
# CI job: unit-no-sirus (standalone `julia --project=ADRIAanalysis test/no_sirus_load.jl`)

using Test
using ADRIAanalysis

@testset "Package loads without SIRUS/MLJ" begin
    @test !any(string(k.name) == "SIRUS" for k in keys(Base.loaded_modules))
    @test !any(string(k.name) == "MLJ" for k in keys(Base.loaded_modules))
    @test isdefined(ADRIAanalysis, :Rule)
    @test isdefined(ADRIAanalysis, :cluster_rules)
    @test isdefined(ADRIAanalysis, :rules)
    @test isnothing(Base.get_extension(ADRIAanalysis, :ADRIAanalysisRulesExt))
    r = ADRIAanalysis.Rule(Vector{Vector}([["f", :L, 0.1f0]]), [0.9, 0.1])
    @test r isa ADRIAanalysis.Rule
    @test_nowarn ADRIAanalysis.print_rules([r])
    @test ADRIAanalysis.maximum_probability([r]) >= 0.0
end
