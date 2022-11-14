using ADRIA, DataFrames, CSV
import GeoDataFrames as GDF
using CSV


@testset "Connectivity loading" begin
    # dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)
    # "C:/development/ADRIA_data/data_packages/Moore_2022-11-02/site_data/Moore_2022-11-02.gpkg"

    site_data = GDF.read(joinpath(@__DIR__, "..", "examples", "Example_domain", "site_data", "Example_domain.gpkg"))
    sort!(site_data, :reef_siteid)

    unique_site_ids = site_data.reef_siteid

    conn_files = joinpath(@__DIR__, "..", "examples", "Example_domain", "connectivity")
    conn_data = CSV.read(joinpath(conn_files, "2000", "example_conn.csv"), DataFrame, comment="#", drop=[1], types=Float64)

    conn_details = ADRIA.site_connectivity(conn_files, unique_site_ids)

    TP_data = conn_details.TP_base
    @test all(names(TP_data, 2) .== site_data.reef_siteid) || "Sites do not match expected order!"
    @test all(unique_site_ids .== conn_details.site_ids) || "Included site ids do not match length/order in geospatial file."
end

@testset "Environmental data" begin
    site_data = GDF.read(joinpath(@__DIR__, "..", "examples", "Example_domain", "site_data", "Example_domain.gpkg"))

    # site_data = GDF.read("C:/development/ADRIA_data/data_packages/Moore_2022-11-02/site_data/Moore_2022-11-02.gpkg")

    sort!(site_data, :reef_siteid)

    wave_fn = joinpath(@__DIR__, "..", "examples", "Example_domain", "waves", "wave_RCP45.nc")
    waves = ADRIA.load_env_data(wave_fn, "Ub", site_data)
    @test all(names(waves, 2) .== site_data.reef_siteid) || "Wave data not aligned with order specified in geospatial data"

    dhw_fn = joinpath(@__DIR__, "..", "examples", "Example_domain", "DHWs", "dhwRCP45.nc")
    dhw = ADRIA.load_env_data(dhw_fn, "dhw", site_data)
    @test all(names(dhw, 2) .== site_data.reef_siteid) || "Wave data not aligned with order specified in geospatial data"
end

@testset "Initial covers" begin
    site_data = GDF.read(joinpath(@__DIR__, "..", "examples", "Example_domain", "site_data", "Example_domain.gpkg"))
    sort!(site_data, :reef_siteid)

    coral_cover_fn = joinpath(@__DIR__, "..", "examples", "Example_domain", "site_data", "coral_cover.nc")
    coral_covers = ADRIA.load_covers(coral_cover_fn, "covers", site_data)

    @test all(names(coral_covers, 2) .== site_data.reef_siteid) || "Coral cover data not aligned with order specified in geospatial data"
end
