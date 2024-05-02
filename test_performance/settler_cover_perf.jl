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

	perf = @benchmark begin 
		# Mock 75 time steps
		for t in 1:75
			ADRIA.settler_cover($args...)
		end
	end

	display(perf)
end

# Benchmark different location sizes
for num_locs in location_nums
	@info "GPU - $(num_locs) locations benchmark"
	local args = generate_data(num_locs)
	fec_scope = CuArray(args[1])
	conn = CuArray(args[2])

	perf = @benchmark begin
		# Mock 75 time steps
		for t in 1:75
			settler_cover_cuda($fec_scope, $conn, $args[3:end]...)
		end
	end

	display(perf)
end
