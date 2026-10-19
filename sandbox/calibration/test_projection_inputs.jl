"""Integration test for the 1991-2100 projection input set v0.5 (build_projection_inputs.py;
crash_mechanism_protocol.md, "Preregistered: projection design to 2100", 2026-10-19).

Every input is parsed the way the model reads it and must equal the v0.4 / W-wave-lf04 input in
1985/1991-2026 exactly; coverage must extend to 2100 (the adapter silently returns zero external
supply and unit chl multipliers for missing years); the closed supply scenario must be zero after 2026.

    julia --project=sandbox sandbox/calibration/test_projection_inputs.jl
"""

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "sandbox"))

using ADRIA, CSV, DataFrames, Test

const V4 = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.4")
const PROJ = joinpath(ROOT, "sandbox", "data", "derived", "projection_inputs_v0.5")
const DATASET = "owen_global_2018_2023_corr_w3"
const CHL = joinpath(ROOT, "sandbox", "data", "derived", "larval_survival_index_lizard_full_1990s_corrected",
    "larval_survival_ap_rma_k.csv")
const MELLIN = joinpath(ROOT, "sandbox", "data", "derived", "mellin_background_mortality_1991_2026", "background_mortality.csv")
const HIST = 1985:2026
const ALL = 1985:2100

ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = DATASET
ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = "cots"
domain = ADRIA.load_domain(ADRIA.LizardDomain, V4, "historical")
loc_ids = domain.loc_ids
groups = ADRIA.functional_group_names()
n_hist = length(HIST)

@testset "projection inputs v0.5" begin
    @testset "DHW" begin
        scenarios = CSV.read(joinpath(PROJ, "dhw", "dhw_scenarios.csv"), DataFrame)
        base = Array(domain.dhw_scens)  # time x location x scenario, 1985-2026
        for ssp in unique(scenarios.ssp)
            cube = ADRIA.YAXArrays.Cube(joinpath(PROJ, "dhw", "dhw_$(ssp).nc"))
            values = Float64.(Array(cube))
            @test collect(cube.time) == collect(ALL)
            @test size(values, 3) == count(==(ssp), scenarios.ssp)
            for s in axes(values, 3)
                @test values[1:n_hist, :, s] == base[:, :, 1]
            end
            @test all(isfinite, values) && all(>=(0.0), values)
        end
    end
    @testset "cyclones and chl (environment replicates)" begin
        base = Array(domain.cyclone_mortality_scens)
        chl_base = CSV.read(CHL, DataFrame; types=Dict(:site_id => String))
        year_map = CSV.read(joinpath(PROJ, "env", "env_year_map.csv"), DataFrame)
        for r in sort(unique(year_map.replicate))
            tag = lpad(r, 2, '0')
            cyc = Array(ADRIA.load_lizard_cyclone_mortality(joinpath(PROJ, "env", "cyclone_mortality_r$(tag).csv"),
                loc_ids, ALL, groups))
            @test cyc[1:n_hist, :, :, :] == base
            m = year_map[year_map.replicate .== r, :]
            @test m.year == collect(2027:2100)
            @test all(1993 .<= m.source_year .<= 2024)
            # A projected year carries its source year's mortality at every site.
            y = m.year[1]; s = m.source_year[1]
            @test cyc[y - 1984, :, 1, 1] == base[s - 1984, :, 1, 1]
            proj = ADRIA.load_cots_larval_survival(joinpath(PROJ, "env", "larval_survival_r$(tag).csv"), loc_ids, 1991:2100)
            hist = ADRIA.load_cots_larval_survival(CHL, loc_ids, 1991:2026)
            @test proj[1:36, :] == hist
            @test proj[y - 1990, :] == hist[s - 1990, :]
        end
    end
    @testset "Mellin background" begin
        old = CSV.read(MELLIN, DataFrame; types=Dict(:site_id => String))
        new = CSV.read(joinpath(PROJ, "mellin", "background_mortality.csv"), DataFrame; types=Dict(:site_id => String))
        @test new.site_id == old.site_id
        @test Matrix{Float64}(new[:, string.(1991:2026)]) == Matrix{Float64}(old[:, string.(1991:2026)])
        @test all(Matrix{Float64}(new[:, string.(2027:2100)]) .== Matrix{Float64}(new[:, ["2027"]])[1])
    end
    @testset "external supply" begin
        base = domain.cots_forcing
        for name in readdir(joinpath(PROJ, "supply"))
            isdir(joinpath(PROJ, "supply", name)) || continue
            f = ADRIA.load_cots_connectivity_forcing(joinpath(PROJ, "supply", name, "annual_forcing"), loc_ids)
            @test f.outbreak_years == collect(1991:2100)
            @test f.external_outbreak_potential[:, 1:36, :] == base.external_outbreak_potential
            @test f.external_outbreak_pelagic[:, 1:36, :] == base.external_outbreak_pelagic
            @test f.reef_connectivity == base.reef_connectivity && f.site_to_reef == base.site_to_reef
            if name == "closed"
                @test all(iszero, f.external_outbreak_potential[:, 37:end, :])
            elseif name == "mean"
                @test f.external_outbreak_potential[:, 37, :] ≈ dropdims(sum(base.external_outbreak_potential; dims=2); dims=2) ./ 36
            end
        end
    end
end
