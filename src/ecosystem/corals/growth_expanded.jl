"""
Not intended for use in production.

The ODE defined here is identical to `growth.jl::growthODE()` in function.
The terms are expanded out to cover each coral taxa/size class.

# Arguments
du : matrix holding derivatives
X : Current coral cover, relative to `k`
p : additional parameters
t : time, unused, so marking with `_`
"""
function growthODE_expanded(du::Array{Float64,2}, X::Array{Float64,2}, p::NamedTuple, _::Real)::Nothing
    # `s` refers to sigma holding leftover space for each site in form of: 1 x n_sites
    s = p.sigma[:, :]

    s .= max.(p.k' .- sum(X, dims=1), 0.0)  # Make relative to k (max. carrying capacity)
    s = vec(s)

    rec = @view p.rec[:, :]  # recruitment values
    r = @view p.r[:]  # growth rates
    X_mb = X .* p.mb   # current cover * background mortality

    # 0.3 is the competition rate
    # Tabular Acropora Enhanced
    du[1, :] .= rec[1, :] .- s .* X[1, :] .* r[1] .- X_mb[1, :]
    du[2, :] .= (s .* X[1, :] .* r[1]) .- (s .* X[2, :] .* r[2]) .- X_mb[2, :]
    du[3, :] .= (s .* X[2, :] .* r[2]) .- (s .* X[3, :] .* r[3]) .- X_mb[3, :]
    du[4, :] .= (s .* X[3, :] .* r[3]) .- (s .* X[4, :] .* r[4]) .- X_mb[4, :]
    du[5, :] .= (s .* X[4, :] .* r[4]) .- (s .* X[5, :]) .* (r[5] .+ 0.3 .* sum(X[[26, 27, 28], :], dims=1)') .- X_mb[5, :]
    du[6, :] .= (s .* X[5, :]) .* (r[5] .+ 0.3 .* sum(X[[26, 27, 28], :], dims=1)') .+ (s .* X[5, :] .* r[5]) .- X_mb[5, :]

    # Tabular Acropora Unenhanced
    du[7, :] .= rec[2, :] .- (s .* X[7, :] .* r[7]) .- X_mb[2, :]
    du[8, :] .= (s .* X[7, :] .* r[7]) .- (s .* X[8, :] .* r[8]) .- X_mb[8, :]
    du[9, :] .= (s .* X[8, :] .* r[8]) .- (s .* X[9, :] .* r[9]) .- X_mb[9, :]
    du[10, :] .= (s .* X[9, :] .* r[9]) .- (s .* X[10, :] .* r[10]) .- X_mb[10, :]
    du[11, :] .= (s .* X[10, :] .* r[10]) .- (s .* X[11, :]) .* (r[11] .+ 0.3 .* sum(X[[26, 27, 28], :], dims=1)') .- X_mb[11, :]
    du[12, :] .= (s .* X[11, :]) .* (r[11] .+ 0.3 .* sum(X[[26, 27, 28], :], dims=1)') .+ (s .* X[11, :] * r[11]) .- X_mb[11, :]

    # Corymbose Acropora Enhanced
    du[13, :] .= rec[3, :] .- s .* X[13, :] .* r[13] .- X_mb[13, :]
    du[14, :] .= (s .* X[13, :] .* r[13]) .- (s .* X[14, :] .* r[14]) .- X_mb[14, :]
    du[15, :] .= (s .* X[14, :] .* r[14]) .- (s .* X[15, :] .* r[15]) .- X_mb[15, :]
    du[16, :] .= (s .* X[15, :] .* r[15]) .- (s .* X[16, :] .* r[16]) .- X_mb[16, :]
    du[17, :] .= (s .* X[16, :] .* r[16]) .- (s .* X[17, :] .* r[17]) .- X_mb[17, :]
    du[18, :] .= (s .* X[17, :] .* r[17]) .+ (s .* X[18, :] .* r[18]) .- X_mb[18, :]

    # Corymbose Acropora Unenhanced
    du[19, :] .= rec[4, :] .- s .* X[19, :] .* r[19] .- X_mb[4, :]
    du[20, :] .= (s .* X[19, :] .* r[19]) .- (s .* X[20, :] .* r[20]) .- X_mb[20, :]
    du[21, :] .= (s .* X[20, :] .* r[20]) .- (s .* X[21, :] .* r[21]) .- X_mb[21, :]
    du[22, :] .= (s .* X[21, :] .* r[21]) .- (s .* X[22, :] .* r[22]) .- X_mb[22, :]
    du[23, :] .= (s .* X[22, :] .* r[22]) .- (s .* X[23, :] .* r[23]) .- X_mb[23, :]
    du[24, :] .= (s .* X[23, :] .* r[23]) .+ (s .* X[24, :] .* r[24]) .- X_mb[24, :]

    # Small massives and encrusting Unenhanced
    du[25, :] .= rec[5, :] .- s .* X[25, :] .* r[25] .- X_mb[25, :]
    du[26, :] .= (s .* X[25, :] .* r[25, :]) .- (s .* X[26, :] .* r[26, :]) .- (X[26, :] .* p.mb[26] .+ 0.3 .* (X[6, :] .+ X[12, :]))
    du[27, :] .= (s .* X[26, :] .* r[26, :]) .- (s .* X[27, :] .* r[27, :]) .- (X[27, :] .* p.mb[27] .+ 0.3 .* (X[6, :] .+ X[12, :]))
    du[28, :] .= (s .* X[27, :] .* r[27, :]) .- (s .* X[28, :] .* r[28, :]) .- (X[28, :] .* p.mb[28] .+ 0.3 .* (X[6, :] .+ X[12, :]))
    du[29, :] .= (s .* X[28, :] .* r[28]) .- (s .* X[29, :] .* r[29]) .- X_mb[29, :]
    du[30, :] .= (s .* X[29, :] .* r[29]) .+ (s .* X[30, :] .* r[30]) .- X_mb[30, :]

    # Large massives Unenhanced
    du[31, :] .= rec[6, :] .- s .* X[31, :] .* r[31] .- X_mb[6, :]
    du[32, :] .= (s .* X[31, :] .* r[31]) .- (s .* X[32, :] .* r[32]) .- X_mb[32, :]
    du[33, :] .= (s .* X[32, :] .* r[32]) .- (s .* X[33, :] .* r[33]) .- X_mb[33, :]
    du[34, :] .= (s .* X[33, :] .* r[33]) .- (s .* X[34, :] .* r[34]) .- X_mb[34, :]
    du[35, :] .= (s .* X[34, :] .* r[34]) .- (s .* X[35, :] .* r[35]) .- X_mb[35, :]
    du[36, :] .= (s .* X[35, :] .* r[35]) .+ (s .* X[36, :] .* r[36]) .- X_mb[36, :]

    # Ensure no non-negative values
    du .= max.(du, 0.0)

    return
end
