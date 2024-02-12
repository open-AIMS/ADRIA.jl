using DimensionalData
using Random
using YAXArrays

@testset "NetCDF" begin
    @testset "Handling Vector of Strings and Matrix{ASCIIChar}" begin
        char_fp = joinpath(TEST_DATA_DIR, "test_ncdf_char.nc")
        str_fp = joinpath(TEST_DOMAIN_PATH, "DHWs", "dhwRCP45.nc")

        ADRIA.NetCDF.open(char_fp; mode=ADRIA.NetCDF.NC_NOWRITE) do char_nc
            char2str_vector = ADRIA._site_labels(char_nc)
            @test char2str_vector isa Vector{String}
            @test char2str_vector[1] == "abcd"
            @test char2str_vector[2] == "something"
            @test char2str_vector[3] == "something else"
        end

        ADRIA.NetCDF.open(char_fp; mode=ADRIA.NetCDF.NC_NOWRITE) do str_nc
            @test ADRIA._site_labels(str_nc) isa Vector{String}
        end
    end
end

@testset "axes_names" begin
    axes_names = ADRIA.axes_names

    @testset "No name dims" begin
        cube = YAXArray(rand(5, 5))

        @test axes_names(cube) == (:Dim_1, :Dim_2)
    end

    @testset "Named dims" begin
        axlist = (Dim{:time}(1:10), Dim{:locations}(1:20))
        data = rand(10, 20)
        cube = YAXArray(axlist, data)

        @test axes_names(cube) == (:time, :locations)
    end
end

@testset "sort_axis" begin
    sort_axis = ADRIA.sort_axis

    locations = shuffle(1:50)
    axlist = (Dim{:timesteps}(1:10), Dim{:locations}(locations))
    cube = YAXArray(axlist, rand(10, 50))
    sorted_cube = sort_axis(cube, :locations)

    @test collect(sorted_cube.locations) == collect(1:50)

    match = true
    for location in locations
        if cube[locations=At(location)] != sorted_cube[locations=At(location)]
            match = false
        end
    end

    @test match == true
end
