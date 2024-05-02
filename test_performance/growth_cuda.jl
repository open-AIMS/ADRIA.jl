using CUDA
using ADRIA

#= 
Temporary location for CUDA variations of ADRIA functions.
TBD if we put these in a Julia extension or separate package.
We should make this decision and move this before merging to main.
=#

#=
Wraps the settler_cover function, converting all arguments to CuArray
=#
function settler_cover_cuda(
    fec_scope::T,
    conn::AbstractMatrix{Float64},
    leftover_space::T,
    α::V,
    β::V,
    basal_area_per_settler::V,
    potential_settlers::T
)::T where {T<:Matrix{Float64},V<:Vector{Float64}}
    fec_scope = CuArray(fec_scope)
    conn = CuArray(conn)
    leftover_space = CuArray(leftover_space)
    α = CuArray(α)
    β = CuArray(β)
    basal_area_per_settler = CuArray(basal_area_per_settler)
    potential_settlers = CuArray(potential_settlers)

    # need to construct these as CuArray and provide them to _settler_cover()
    # settler_cover's corresponding BitVector lines causes LoadError: Scalar indexing is disallowed.
    valid_sources = CuArray(vec(sum(conn, dims=2) .> 0.0))
    valid_sinks = CuArray(vec(sum(conn, dims=1) .> 0.0))

    return ADRIA._settler_cover(fec_scope, conn, leftover_space, α, β, basal_area_per_settler, potential_settlers, valid_sources, valid_sinks)
end
