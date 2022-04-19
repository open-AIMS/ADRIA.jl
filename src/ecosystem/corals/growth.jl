"""Growth ODE functions"""


"""
    growthODE(du, X, p, _)

Base coral growth function.
"""
function growthODE(du::Array{Float64, 2}, X::Array{Float64, 2}, p::NamedTuple, _::Real)::Nothing
    k = @view p.k[:, :]
    k .= max.(p.P .- sum(X, dims=1), 0.0)

    # Use temporary caches
    k_X_r = @view p.kXr[:, :]
    k_rec = @view p.k_rec[:, :]
    X_mb = @view p.X_mb[:, :]
    kX_sel_en = @view p.kX_sel_en[:, :]
    X_tab = @view p.X_tab[:, :]
    @. k_X_r = k * X * p.r
    @. k_rec = k * p.rec
    @. X_mb = X * p.mb

    @views @. kX_sel_en = k * X[p.sel_en, :]
    @views @. X_tab = (p.mb[26] .+ p.comp .* (X[6, :] .+ X[12, :])')

    @views r_comp = (p.r[p.sel_en] .+ (p.comp .* sum(X[p.small_massives, :])))

    @views @. du[p.sel_en, :] = k_X_r[p.sel_en - 1, :] - kX_sel_en * r_comp - X_mb[p.sel_en, :]
    @views @. du[p.sel_unen, :] = kX_sel_en * r_comp + k_X_r[p.sel_unen, :] - X_mb[p.sel_unen, :]

    @views @. du[p.encrusting, :] = k_rec[p.enc, :] - k_X_r[p.encrusting, :] - X_mb[p.encrusting, :]

    @views @. du[p.small_massives, :] = k_X_r[p.small_massives - 1, :] - k_X_r[p.small_massives, :] - X[p.small_massives, :] * X_tab

    @views @. du[p.small, :] = k_rec[p.small_r, :] - k_X_r[p.small, :] - X_mb[p.small, :]
    @views @. du[p.mid, :] = k_X_r[p.mid - 1, :] - k_X_r[p.mid, :] - X_mb[p.mid, :]
    @views @. du[p.large, :] = k_X_r[p.large - 1 , :] + k_X_r[p.large, :] - X_mb[p.large, :]

    # Ensure no non-negative values
    du .= max.(du, 0.0)

    return
end