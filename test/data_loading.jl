using ADRIA, DataFrames, CSV, NamedDims, AxisKeys
import GeoDataFrames as GDF
using CSV


if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "examples", "Example_domain")
end


@testset "Connectivity loading" begin
    location_data = GDF.read(joinpath(EXAMPLE_DOMAIN_PATH, "location_data", "Example_domain.gpkg"))
    sort!(location_data, :reef_siteid)

    unique_location_ids = location_data.reef_siteid

    conn_files = joinpath(EXAMPLE_DOMAIN_PATH, "connectivity")
    conn_data = CSV.read(joinpath(conn_files, "2000", "example_conn.csv"), DataFrame, comment="#", drop=[1], types=Float64)

    conn_details = ADRIA.location_connectivity(conn_files, unique_location_ids)

    TP_data = conn_details.TP_base
    @test all(axiskeys(TP_data, 1) .== axiskeys(TP_data, 2)) || "Site order does not match between rows/columns."
    @test all(axiskeys(TP_data, 2) .== location_data.reef_siteid) || "Sites do not match expected order."
    @test all(unique_location_ids .== conn_details.location_ids) || "Included location ids do not match length/order in geospatial file."
end

@testset "Environmental data" begin
    location_data = GDF.read(joinpath(EXAMPLE_DOMAIN_PATH, "location_data", "Example_domain.gpkg"))

    sort!(location_data, :reef_siteid)

    wave_fn = joinpath(EXAMPLE_DOMAIN_PATH, "waves", "wave_RCP45.nc")
    waves = ADRIA.load_env_data(wave_fn, "Ub", location_data)
    @test all(axiskeys(waves, 2) .== location_data.reef_siteid) || "Wave data not aligned with order specified in geospatial data"

    dhw_fn = joinpath(EXAMPLE_DOMAIN_PATH, "DHWs", "dhwRCP45.nc")
    dhw = ADRIA.load_env_data(dhw_fn, "dhw", location_data)
    @test all(axiskeys(dhw, 2) .== location_data.reef_siteid) || "Wave data not aligned with order specified in geospatial data"
end

@testset "Initial covers" begin
    location_data = GDF.read(joinpath(EXAMPLE_DOMAIN_PATH, "location_data", "Example_domain.gpkg"))
    sort!(location_data, :reef_siteid)

    coral_cover_fn = joinpath(EXAMPLE_DOMAIN_PATH, "location_data", "coral_cover.nc")
    coral_covers = ADRIA.load_covers(coral_cover_fn, "covers", location_data)

    @test all(axiskeys(coral_covers, 2) .== location_data.reef_siteid) || "Coral cover data not aligned with order specified in geospatial data"
end
