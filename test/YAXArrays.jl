using DimensionalData
using Random
using YAXArrays

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

@testset "sort" begin
    sort = ADRIA.sort

    locations = shuffle(1:50)
    axlist = (Dim{:timesteps}(1:10), Dim{:locations}(locations))
    cube = YAXArray(axlist, rand(10, 50))
    sorted_cube = sort(cube, :locations)

    @test collect(sorted_cube.locations) == collect(1:50)

    match = true
    for location in locations
        if cube[locations=At(location)] != sorted_cube[locations=At(location)]
            match = false
        end
    end

    @test match == true
end
