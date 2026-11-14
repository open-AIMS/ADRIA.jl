using ADRIA.CSV, ADRIA.DataFrames

if !@isdefined(TEST_DOMAIN_PATH)
    const TEST_DOMAIN_PATH = joinpath(pkgdir(ADRIA), "test", "data", "Test_domain")
end

"""
    make_cots_test_package(dest::String)::String

Copy Test_domain to `dest` (whose folder name must be `Test_domain`, which names the
package's geopackage) and add every optional COTS input (see `ADRIA.CotsInputs`):
locations 1-5 form reef A (50 ha) and 6-10 reef B (80 ha); two annual reef matrices
(hydrodynamic years 2001 and 2002); outside outbreak supply in 2026 and 2027.
"""
function make_cots_test_package(dest::String)::String
    cp(TEST_DOMAIN_PATH, dest)
    conn_fn = only(readdir(joinpath(dest, "connectivity"); join=true))
    ids = names(CSV.read(conn_fn, DataFrame; comment="#"))[2:end]
    n = length(ids)

    conn_dir = mkpath(joinpath(dest, "cots_connectivity"))
    conn = DataFrame(; source=ids)
    for (j, id) in enumerate(ids)
        conn[!, id] = [i == j ? 0.3 : (j == mod1(i + 1, n) ? 0.05 : 0.0) for i in 1:n]
    end
    CSV.write(joinpath(conn_dir, "cots_conn.csv"), conn)

    forcing_dir = mkpath(joinpath(conn_dir, "annual_forcing"))
    reef = [i <= 5 ? 1 : 2 for i in 1:n]
    CSV.write(
        joinpath(forcing_dir, "site_reef_map.csv"),
        DataFrame(; site_id=ids, reef_index=reef,
            reef_habitable_area_ha=[r == 1 ? 50.0 : 80.0 for r in reef])
    )
    for (year, m) in ((2001, [0.2 0.1; 0.05 0.3]), (2002, [0.1 0.2; 0.2 0.1]))
        CSV.write(
            joinpath(forcing_dir, "reef_connectivity_$(year).csv"),
            DataFrame(; reef_id=["A", "B"], A=m[:, 1], B=m[:, 2])
        )
    end
    boundary = DataFrame(;
        reef_id=String[], calendar_year=Int[], hydrodynamic_year=Int[],
        probability_area_connectivity=Float64[]
    )
    for (r, id) in enumerate(["A", "B"]), cal in (2026, 2027), hydro in (2001, 2002)
        push!(boundary, (id, cal, hydro, 1e-4 * r * (hydro - 2000)))
    end
    CSV.write(joinpath(forcing_dir, "external_outbreak_boundary.csv"), boundary)
    write(
        joinpath(forcing_dir, "provenance.toml"),
        """
        years = [2001, 2002]
        site_map = "site_reef_map.csv"
        matrix_files = ["reef_connectivity_2001.csv", "reef_connectivity_2002.csv"]
        external_outbreak_boundary = "external_outbreak_boundary.csv"
        """
    )

    CSV.write(
        joinpath(dest, "cots_initial_state.csv"),
        DataFrame(; site_id=ids, juvenile_ha=collect(1.0:n), subadult_ha=fill(0.5, n),
            adult_ha=fill(2.0, n))
    )
    CSV.write(
        joinpath(dest, "cots_larval_survival.csv"),
        DataFrame(; site_id=[ids[1], ids[2], ids[3]], year=[2026, 2030, 1990],
            multiplier=[0.5, 0.0, 9.0])
    )
    CSV.write(
        joinpath(dest, "cots_settlement_modifier.csv"),
        DataFrame(; site_id=ids, multiplier=[i == n ? 0.5 : 0.8 for i in 1:n])
    )
    return dest
end
