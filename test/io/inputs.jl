@testset "NetCDF" begin
    @testset "Handling Vector of Strings and Matrix{ASCIIChar}" begin
        char_fp = joinpath(TEST_DATA_DIR, "test_ncdf_char.nc")
        str_fp = joinpath(EXAMPLE_DOMAIN_PATH, "DHWs", "dhwRCP45.nc")

        ADRIA.NetCDF.open(char_fp; mode=ADRIA.NetCDF.NC_NOWRITE) do char_nc
            @test ADRIA._site_labels(char_nc) isa Vector{String}
        end

        ADRIA.NetCDF.open(char_fp; mode=ADRIA.NetCDF.NC_NOWRITE) do str_nc
            @test ADRIA._site_labels(str_nc) isa Vector{String}
        end
    end
end
