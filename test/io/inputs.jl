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

@testset "ZeroDataCube" begin

    ZeroDataCube = ADRIA.ZeroDataCube
    
    dim_1_name = :timesteps
    dim_1_vals = rand(10)
    dim_1 = Dim{dim_1_name}(dim_1_vals)
    
    dim_2_name = :location
    dim_2_vals = ["loc 1", "loc 2", "loc 3", "loc 4", "loc 5"]
    dim_2 = Dim{dim_2_name}(dim_2_vals)

    yax_res = ZeroDataCube(Int, timesteps=dim_1_vals)

    @test typeof(yax_res) <: YAXArray{Int, 1} || 
        "Incorrect return type. Expected a subtype of YAXArray{Int, 1} \
         but received $(typeof(yax_res))"

    @test all(0 .== yax_res) || 
        "Incorrect data contained in YAXArray. Expected all 0."

    @test size(yax_res) == (10,) ||
        "Incorrect YAXArray shape. Expected (10,) but received $(size(yax_res))"

    @test name(dim_1) == name(yax_res.axes[1]) ||
        "Incorrect dimension name. Expected $(name(dim_1)) \
         but received $(name(yax_res.axes[1]))"

    @test dim_1_vals == collect(yax_res.axes[1]) || 
        "Incorrect axis indices. Expected $(dim_1_vals) but received $(collect(yax_res.axes[1]))"
      
    yax_res = ZeroDataCube(Float64, timesteps=dim_1_vals, location=dim_2_vals)

    @test typeof(yax_res) <: YAXArray{Float64, 2} || 
        "Incorrect return type. Expected a subtype of YAXArray{Int, 2} \
         but received $(typeof(yax_res))"

    @test all(0 .== yax_res) || 
        "Incorrect data contained in YAXArray. Expected all 0."

    @test size(yax_res) == (10, 5) ||
        "Incorrect YAXArray shape. Expected (10, 5) but received $(size(yax_res))"

    @test name(dim_1) == name(yax_res.axes[1]) ||
        "Incorrect dimension name. Expected $(name(dim_1)) \
         but received $(name(yax_res.axes[1]))"

    @test name(dim_2) == name(yax_res.axes[2]) ||
        "Incorrect dimension name. Expected $(name(dim_2)) \
         but received $(name(yax_res.axes[2]))"

    @test dim_1_vals == collect(yax_res.axes[1]) || 
        "Incorrect axis indices. Expected $(dim_1_vals) but received $(collect(yax_res.axes[1]))"

    @test dim_2_vals == collect(yax_res.axes[2]) || 
        "Incorrect axis indices. Expected $(dim_2_vals) but received $(collect(yax_res.axes[2]))"
end
