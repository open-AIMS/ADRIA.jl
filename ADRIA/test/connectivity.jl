using ADRIA, ADRIA.DataFrames, ADRIA.CSV
import ADRIA.GDF as GDF

if !@isdefined(TEST_DOMAIN_PATH)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DATA_DIR = joinpath(ADRIA_DIR, "test", "data")
    const TEST_DOMAIN_PATH = joinpath(TEST_DATA_DIR, "Test_domain")
end

@testset "Connectivity loading" begin
    loc_data = GDF.read(
        joinpath(TEST_DOMAIN_PATH, "spatial", "Test_domain.gpkg")
    )
    sort!(loc_data, [:reef_siteid])

    unique_loc_ids = loc_data.reef_siteid
    conn_files = joinpath(TEST_DOMAIN_PATH, "connectivity")
    conn_data = CSV.read(
        joinpath(conn_files, "example_conn.csv"),
        DataFrame;
        comment="#",
        drop=[1],
        types=Float64
    )

    conn_details = ADRIA.location_connectivity(conn_files, unique_loc_ids)

    conn = conn_details.conn
    d1, d2 = axes(conn)
    @test all(d1.dim.val.data .== d2.dim.val.data) ||
        "Site order does not match between rows/columns."
    @test all(d2.dim.val.data .== loc_data.reef_siteid) ||
        "Sites do not match expected order!"
end
