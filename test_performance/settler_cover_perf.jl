
using ADRIA
using BenchmarkTools
using Test
using Random

const num_corals = 6

#=
Generate random values from 0.0eX to 1.0eY (X = e_min, Y=e_max)

This should work, but doesn't have the beahavior I'd expect:
rand(1.0e6:9.9e11, (15, 16))
it mostly generates values around e11??
=#
function generate_matrix(nrows::Int, ncols::Int, e_min::Int, e_max::Int)
	m = rand(nrows, ncols)
	# increment since rand is generating values 0.0:1.0
	e_min += 1
	e_max += 1
	# discrete e values
	e_scales = rand([10^x for x in e_min:e_max], nrows, ncols)
	return m .* e_scales
end

#=
Generate arguments for the setter_cover function.
data 

returns (fec_scope, conn, leftover_space, alpha, beta, basal_area_per_settler, potential_settlers)
=#
function generate_data(num_locs::Int)
	# fec_scope
	# see cache.fec_scope and fecundity_scope!()
	# num_corals, num_locs = size(fec_scope)
	# generate values of a similar magnitude 1e6:1e10
	fec_scope = generate_matrix(num_corals, num_locs, 9, 11)

	# conn - square of locations 0.00091:0.067
	# see domain.conn
	# ~98% of values are 0.0
	# TODO should set values to 0.0, not sure if it matters
	conn = rand(num_locs, num_locs) ./ 100

	# leftover_space 1,1 to 1,10 values[2709:85_760]
	#leftover_space = [[2709.3686645488638, 8606.852496743664, 3129.1949880571137, 4540.825085215324, 85760.03954534238, 81485.9754509628, 2855.635238583546, 40542.730053980136, 18016.42807901491, 6802.382677733639]]
	leftover_space = rand(2500.0:86_000.0, 1, num_locs)

	# Coral settings from test small spec
	alpha = [0.75, 3.75, 3.75, 2.25, 2.25, 2.25]

	beta = [5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0]

	basal_area_per_settler = [7.823175026020862e-5, 8.208148940404582e-5, 7.375108157700818e-5, 8.265213751220709e-5, 8.057901929322589e-5, 8.176323743537886e-5]

	# potential_settlers: either 0.0, or large number
	# in test small spec, some loactions 0.0
	# example row: [0.0 0.0 0.0 4.4097883711344564e8 4.937007418080003e8 1.8201398503492028e8 1.5191456268188483e8 0.0 0.0 0.0],
	potential_settlers = generate_matrix(num_corals, num_locs, 8, 10)

	return (fec_scope, conn, leftover_space, alpha, beta, basal_area_per_settler, potential_settlers)
end

@testset "check settler_cover results" begin
	for i in 1:3
		num_locs = 10^i
		println("$(num_locs) locations")
		args = generate_data(num_locs)

		result = ADRIA.settler_cover(args...)
		@test size(result) == (num_corals, num_locs)
		@test all(v -> v > 0, result)
	end
end

# Note: run with Moore_2024-02-14_v060_rc1 has 334 locations
# Benchmark different location sizes
for i in 1:3
	num_locs = 10^i
	println("$(num_locs) locations benchmark")
	args = generate_data(num_locs)

	display(@benchmark ADRIA.settler_cover($args...))
end
