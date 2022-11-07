using ADRIA, DataFrames, CSV
import GeoDataFrames as GDF
using CSV


@testset "Connectivity loading" begin
    # dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)
    # "C:/development/ADRIA_data/data_packages/Moore_2022-11-02/site_data/Moore_2022-11-02.gpkg"

    site_data = GDF.read(joinpath(@__DIR__, "..", "examples", "Example_domain", "site_data", "Example_domain.gpkg"))
    sort!(site_data, [:reef_siteid])

    conn_ids = site_data.site_id
    unique_site_ids = site_data.reef_siteid

    conn_files = joinpath(@__DIR__, "..", "examples", "Example_domain", "connectivity")
    conn_data = CSV.read(joinpath(conn_files, "2000", "test_conn_data.csv"), DataFrame, comment="#", drop=[1], types=Float64)
    # con_site_ids = names(conn_data)
    conn_ids = site_data.site_id

    conn_details = ADRIA.site_connectivity(conn_files, conn_ids, unique_site_ids)

    TP_data = conn_details.TP_base
    @test all(names(TP_data, 2) .== site_data.reef_siteid) || "Sites do not match expected order!"
    @test all(conn_ids .== conn_details.site_ids) || "Included site ids do not match length/order in geospatial file."
end
