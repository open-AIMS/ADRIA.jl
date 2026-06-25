using ADRIA, ADRIA.DataFrames, ADRIA.CSV
import ADRIA.GeoDataFrames as GDF

@testset "Connectivity loading" begin
    loc_data = GDF.read(
        joinpath(@__DIR__, "..", "examples", "Test_domain", "spatial", "Test_domain.gpkg")
    )
    sort!(loc_data, [:reef_siteid])

    unique_loc_ids = loc_data.reef_siteid
    conn_files = joinpath(@__DIR__, "..", "examples", "Test_domain", "connectivity")
    conn_data = CSV.read(
        joinpath(conn_files, "2000", "test_conn_data.csv"),
        DataFrame;
        comment="#",
        drop=[1],
        types=Float64
    )

    conn_details = ADRIA.location_connectivity(conn_files, unique_loc_ids)

    conn = conn_details.conn
    @test all(names(conn, 2) .== loc_data.reef_siteid) ||
        "Sites do not match expected order!"
end
