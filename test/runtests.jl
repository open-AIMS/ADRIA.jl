using ADRIA
using Test


@testset "ADRIA.jl" begin
    # Write your tests here.
end


@testset "site_selection" begin
    site_path = joinpath("./data", "test_site_data.gpkg")
    conn_path = joinpath("./data", "test_conn_data.csv")

    test_domain = Domain(
        "Test",
        site_path,
        "siteref",
        "reef_siteid",
        "",            # empty coral cover
        conn_path,     # test connectivity data
        "",            # empty DHW
        ""             # empty wave
    );

    p_tbl = ADRIA.param_table(test_domain)
    p_tbl.depth_offset .= 7.0
    ranks = ADRIA.site_selection(test_domain, p_tbl, 1, 10, 1)

end
