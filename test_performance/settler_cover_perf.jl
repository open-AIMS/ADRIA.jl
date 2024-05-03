using ADRIA
using BenchmarkTools
using Test
using Random

include("generate_data.jl")
include("growth_cuda.jl")

# Note: run with Moore_2024-02-14_v060_rc1 has 334 locations
# may have 6,000 max
location_nums = [3_000]

@testset "check settler_cover results" begin
    for num_locs in location_nums
        println("$(num_locs) locations")
        local args = generate_data(num_locs)

        result = ADRIA.settler_cover(args...)
        @test size(result) == (num_corals, num_locs)
        @test all(v -> v > 0, result)

        # TODO run settler_cover_cuda and verify same output
    end
end

# Benchmark different location sizes
for num_locs in location_nums
    @info "CPU - $(num_locs) locations benchmark"
    local args = generate_data(num_locs)

    cpu_bm = @benchmark ADRIA.settler_cover($args...)
    display(cpu_bm)

    @info "GPU - $(num_locs) locations benchmark"
    local (fec_scope, conn, leftover_space, alpha, beta, basal_area_per_settler, potential_settlers) = generate_data(
        num_locs
    )
    c_fec_scope = CuArray(fec_scope)
    c_conn = CuArray(conn)
    c_potential_settlers = CuArray(potential_settlers)

    valid_sources = CuArray(vec(sum(conn; dims=2) .> 0.0))
    valid_sinks = CuArray(vec(sum(conn; dims=1) .> 0.0))

    gpu_bm = @benchmark begin
        settler_cover_cuda(
            $c_fec_scope,
            $c_conn,
            $leftover_space,
            $alpha,
            $beta,
            $basal_area_per_settler,
            $c_potential_settlers,
            $valid_sources,
            $valid_sinks
        )
    end

    display(gpu_bm)

    @info "ratio of median CPU/GPU"
    display(ratio(median(cpu_bm), median(gpu_bm)))

    @info "GPU settler_cover_cuda2"
    gpu2_bm = @benchmark begin
        settler_cover_cuda2(
            $c_fec_scope,
            $c_conn,
            $leftover_space,
            $alpha,
            $beta,
            $basal_area_per_settler,
            $potential_settlers
        )
    end

    display(gpu2_bm)

    #= 	recruitment_rate ~500 microseconds, negligible

    @info "recruitment_rate no multiply"
    rr_bm = @benchmark begin
        ADRIA.recruitment_rate($potential_settlers, $leftover_space; α=$alpha, β=$beta)
    end
    display(rr_bm)

	@info "recruitment_rate with multiply"
	rr_bm = @benchmark begin
        ADRIA.recruitment_rate($potential_settlers, $leftover_space; α=$alpha, β=$beta) .*
        $basal_area_per_settler
    end
    display(rr_bm)
    =#
end
