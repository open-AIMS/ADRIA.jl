using Revise

includet("generate_data.jl")
includet("growth_cuda.jl")

args = generate_data(1_000)

settler_cover_cuda(args...)
