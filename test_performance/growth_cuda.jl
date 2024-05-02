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

    # Potential improvement: make `conn` and `fec_scope` arguments CuArrays
    # which means creation of indices will run on GPU and subsetting will just work.
    # cu_conn = CUDA.rand(1000, 1000)
    # valid_sources = vec(sum(cu_conn, dims=2) .> 0.0)
    # valid_sinks = vec(sum(cu_conn, dims=1) .> 0.0)
    valid_sources = vec(sum(conn, dims=2) .> 0.0)
    valid_sinks = vec(sum(conn, dims=1) .> 0.0)

    c_fec_scope = CuArray(fec_scope[:, valid_sources])
    c_conn = CuArray(conn[valid_sources, valid_sinks])

    # Calculate settler cover and copy result back to host device (i.e., RAM)
    # This matrix multiplication is the most time-consuming part
    # (`recruitment_rate()` takes < 1ms)
    copyto!(c_fec_scope * c_conn, potential_settlers[:, valid_sinks]) 

    return ADRIA.recruitment_rate(potential_settlers, leftover_space; α=α, β=β) .* basal_area_per_settler
end
