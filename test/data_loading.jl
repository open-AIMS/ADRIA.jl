using ADRIA, ADRIA.DataFrames, ADRIA.CSV
using ADRIA.YAXArrays
import ADRIA.GDF as GDF

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DATA_DIR = joinpath(ADRIA_DIR, "test", "data")
    const TEST_DOMAIN_PATH = joinpath(TEST_DATA_DIR, "Test_domain")
end

@testset "Domain loading" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
    @test dom isa Domain
    @test all(dom.dhw_scens .== 0.0)

    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
    @test dom isa Domain
    @test maximum(dom.dhw_scens) > 0.0
end

@testset "Connectivity loading" begin
    site_data = GDF.read(joinpath(TEST_DOMAIN_PATH, "spatial", "Test_domain.gpkg"))
    sort!(site_data, :reef_siteid)

    unique_site_ids = site_data.reef_siteid

    conn_files = joinpath(TEST_DOMAIN_PATH, "connectivity")
    conn_data = CSV.read(
        joinpath(conn_files, "example_conn.csv"),
        DataFrame;
        comment="#",
        drop=[1],
        types=Float64
    )

    conn_details = ADRIA.site_connectivity(conn_files, unique_site_ids)

    conn = conn_details.conn
    d1, d2 = axes(conn)
    @test all(d1.dim .== d2.dim) || "Site order does not match between rows/columns."
    @test all(d2.dim .== site_data.reef_siteid) || "Sites do not match expected order."
    @test all(unique_site_ids .== conn_details.site_ids) ||
        "Included site ids do not match length/order in geospatial file."
end

@testset "Environmental data" begin
    site_data = GDF.read(joinpath(TEST_DOMAIN_PATH, "spatial", "Test_domain.gpkg"))

    sort!(site_data, :reef_siteid)

    wave_fn = joinpath(TEST_DOMAIN_PATH, "waves", "wave_RCP45.nc")
    waves = ADRIA.load_env_data(wave_fn, "Ub")
    @test all(axes(waves, 2).dim .== site_data.reef_siteid) ||
        "Wave data not aligned with order specified in geospatial data"

    dhw_fn = joinpath(TEST_DOMAIN_PATH, "DHWs", "dhwRCP45.nc")
    dhw = ADRIA.load_env_data(dhw_fn, "dhw")
    @test all(axes(dhw, 2).dim .== site_data.reef_siteid) ||
        "Wave data not aligned with order specified in geospatial data"
end

@testset "Initial covers" begin
    site_data = GDF.read(joinpath(TEST_DOMAIN_PATH, "spatial", "Test_domain.gpkg"))
    sort!(site_data, :reef_siteid)

    coral_cover_fn = joinpath(TEST_DOMAIN_PATH, "spatial", "coral_cover.nc")
    coral_covers = ADRIA.load_initial_cover(coral_cover_fn)

    @test all(axes(coral_covers, 2).dim .== site_data.reef_siteid) ||
        "Coral cover data not aligned with order specified in geospatial data"
end

@testset "Cyclone mortality data" begin
    site_data = GDF.read(joinpath(TEST_DOMAIN_PATH, "spatial", "Test_domain.gpkg"))
    sort!(site_data, :reef_siteid)

    cyclone_mortality_fn = joinpath(TEST_DOMAIN_PATH, "cyclones", "cyclone_mortality.nc")
    cyclone_mortality = ADRIA.load_cyclone_mortality(cyclone_mortality_fn)

    expected_species_order = [
        "abhorescent_acropora",
        "tabular_acropora",
        "corymbose_acropora",
        "corymbose_non_acropora",
        "small_massives",
        "large_massives"
    ]

    @test all(axes(cyclone_mortality, 2).dim .== site_data.reef_siteid) ||
        "Cyclone mortality locations do not align with location order specified in geospatial data"
    @test all(axes(cyclone_mortality, 3).dim .== expected_species_order) ||
        "Cyclone mortality data does not list species in expected order"
end
