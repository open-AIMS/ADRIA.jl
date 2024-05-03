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
    fec_scope::CA,
    conn::CA,
    leftover_space::T,
    α::V,
    β::V,
    basal_area_per_settler::V,
    potential_settlers::CA,
    valid_sources::CuArray{Bool},
    valid_sinks::CuArray{Bool}
)::T where {T<:Matrix{Float64},V<:Vector{Float64},CA<:CuArray{Float64}}
    potential_settlers[:, valid_sinks] .= (
        fec_scope[:, valid_sources] * conn[valid_sources, valid_sinks]
    )

    # move from GPU to CPU
    potential_settlers = Array(potential_settlers)

    return ADRIA.recruitment_rate(potential_settlers, leftover_space; α=α, β=β) .*
           basal_area_per_settler
end

function settler_cover_cuda2(
    fec_scope::CuArray{Float64},
    conn::CuArray{Float64},
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

    c_fec_scope = fec_scope[:, valid_sources]
    c_conn = conn[valid_sources, valid_sinks]

    # Copy index of valid sinks from GPU to host
    # `sink_idx` could be preallocated or some other optimization...
    sink_idx = zeros(Bool, length(valid_sinks))
    copyto!(valid_sinks, sink_idx)

    # Calculate settler cover and copy result back to host
    # This matrix multiplication is the most time-consuming part
    # (`recruitment_rate()` takes < 1ms)
    # FIXME reversed? dest param is first.
    copyto!(c_fec_scope * c_conn, view(potential_settlers, :, sink_idx))

    return ADRIA.recruitment_rate(potential_settlers, leftover_space; α=α, β=β) .* basal_area_per_settler
end