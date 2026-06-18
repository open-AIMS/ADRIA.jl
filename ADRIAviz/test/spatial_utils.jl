using Test
using DataFrames
using ADRIAviz.viz

@testset "Spatial Utilities" begin
    @testset "Haversine Distance" begin
        # Test basic functionality with known distances
        # San Francisco to Los Angeles (approximately 559 km)
        sf_lon, sf_lat = -122.419, 37.775
        la_lon, la_lat = -118.243, 34.052
        dist_sf_la = _haversine_km(sf_lon, sf_lat, la_lon, la_lat)
        @test isapprox(dist_sf_la, 559; atol=15)

        # Same point should be zero distance
        @test isapprox(_haversine_km(0.0, 0.0, 0.0, 0.0), 0.0; atol=1e-9)

        # Equator: 1 degree longitude ≈ 111.32 km at equator
        @test isapprox(_haversine_km(0.0, 0.0, 1.0, 0.0), 111.32; atol=1.0)

        # Prime meridian: 1 degree latitude ≈ 111.32 km anywhere
        @test isapprox(_haversine_km(0.0, 0.0, 0.0, 1.0), 111.32; atol=0.5)

        # Test symmetry: distance should be the same regardless of direction
        d1 = _haversine_km(10.0, 20.0, 30.0, 40.0)
        d2 = _haversine_km(30.0, 40.0, 10.0, 20.0)
        @test isapprox(d1, d2; atol=1e-9)

        # Test with negative coordinates (western hemisphere, southern hemisphere)
        d_neg = _haversine_km(-120.0, -30.0, -100.0, -20.0)
        @test d_neg > 0  # Distance should be positive

        # Antipodal points (opposite sides of Earth, ~20015 km)
        # Use points not exactly at poles to avoid singularities
        d_antipodal = _haversine_km(0.0, 45.0, 180.0, -45.0)
        @test isapprox(d_antipodal, 14142; atol=100)  # Great circle distance
    end

    @testset "Nice Scale Bar Length" begin
        # Test rounding to nice values
        @test _nice_length(50.0) == 50
        @test _nice_length(37.0) == 25  # Round down to nearest nice value
        @test _nice_length(0.5) == 1    # Minimum value
        @test _nice_length(150.0) == 100  # Round down from 150
        @test _nice_length(250.0) == 200  # Round down from 250
        @test _nice_length(1500.0) == 1000  # Largest value

        # Test edge cases
        @test _nice_length(0.0) == 1  # Very small values
        @test _nice_length(0.1) == 1
        @test _nice_length(2000.0) == 1000  # Beyond max, should return max

        # Test that function is monotonically non-decreasing
        prev = _nice_length(1.0)
        for val in [2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0]
            curr = _nice_length(val)
            @test curr >= prev
            prev = curr
        end
    end

    @testset "Location ID Resolution" begin
        # Test _loc_id_col is defined and accessible
        @test isdefined(ADRIAviz.viz, :_loc_id_col)

        # Test _get_site_ids fallback
        df_simple = DataFrames.DataFrame(site_id=[1, 2, 3])
        ids = _get_site_ids(df_simple)
        @test ids == ["1", "2", "3"]

        # Test with multiple candidate columns (should warn but use first present)
        df_ambiguous = DataFrames.DataFrame(
            reef_siteid=[10, 20, 30],
            site_id=[1, 2, 3]
        )
        @test_logs (:warn,) ids_ambig = _get_site_ids(df_ambiguous)
        @test length(ids_ambig) == 3

        # Test fallback to row indices when no site_id columns present
        df_noind = DataFrames.DataFrame(data=[1, 2, 3])
        ids_fallback = _get_site_ids(df_noind)
        @test ids_fallback == ["1", "2", "3"]
    end

    @testset "MapDecorationData Structure" begin
        # Test that MapDecorationData can be created and accessed
        places = [(name="Test", lon=150.0, lat=-20.0)]
        deco = MapDecorationData(
            places, 50, 0.45, 145.0, -25.0,
            [145.0, 155.0], [-25.0, -15.0]
        )

        @test deco.scale_bar_km == 50
        @test deco.scale_bar_deg ≈ 0.45
        @test deco.scale_bar_x0 ≈ 145.0
        @test deco.scale_bar_y ≈ -25.0
        @test length(deco.lon_range) == 2
        @test length(deco.lat_range) == 2
        @test length(deco.places) == 1
    end

    @testset "Coastal Places Constant" begin
        # Test that GBR_COASTAL_PLACES is accessible from viz module
        @test :GBR_COASTAL_PLACES in names(ADRIAviz.viz; all=true)
        places = ADRIAviz.viz.GBR_COASTAL_PLACES
        @test length(places) > 0

        # Each entry should be a 3-tuple of (name, lon, lat)
        for place in places
            @test length(place) == 3
            @test isa(place[1], String)  # name
            @test isa(place[2], Number)  # lon
            @test isa(place[3], Number)  # lat

            # Check reasonable coordinate ranges for Australia
            @test place[2] > 140 && place[2] < 155  # Longitude
            @test place[3] < -10 && place[3] > -30  # Latitude
        end
    end

    @testset "Coordinate System Assumptions" begin
        # Test that haversine produces reasonable results for GBR region
        # Brisbane to Cairns is roughly 1400 km
        brisbane = (153.025, -27.470)
        cairns = (145.770, -16.920)
        dist_brisbane_cairns = _haversine_km(brisbane..., cairns...)

        @test 1350 < dist_brisbane_cairns < 1450

        # Verify that larger distances scale appropriately
        # (not testing exact values, just that haversine is working correctly)
        d1 = _haversine_km(145.0, -20.0, 145.0, -21.0)  # 1 degree latitude
        d2 = _haversine_km(145.0, -20.0, 145.0, -22.0)  # 2 degrees latitude
        @test isapprox(d2 / d1, 2.0; atol=0.01)  # Should be 2x distance
    end
end
