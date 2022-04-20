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
    @views @. X_tab = (p.mb[26] + p.comp * (X[6, :] + X[12, :])')

    r_comp = @views p.r[p.sel_en] .+ (p.comp .* sum(X[p.small_massives, :]))

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


"""
    slow_ODE(Y, X, p, t)

Slow version of the growth model, ported directly from MATLAB.

Proportions of corals within a size class transitioning to the next size
class up (r) is based on the assumption that colony sizes within each size
bin are evenly distributed within bins. Transitions are then a simple
ratio of the change in colony size to the width of the bin. See
coralParms for further explanation of these coral metrics.

Note that recruitment pertains to coral groups (n = 6) and represents
the contribution to the cover of the smallest size class within each
group.  While growth and mortality metrics pertain to groups (6) as well
as size classes (6) across all sites (total of 36 by nsites), recruitment is
a 6 by nsites array.

Reshape flattened input from ODE back to expected matrix shape
Dims: (coral species, sites)
"""
function slow_ODE(Y, X, p, t)

    r, P, mb, rec, comp = p

    ## Density dependent growth and recruitment
    # P - sum over coral covers within each site
    # This sets the carrying capacity k := 0.0 <= k <= P
    # resulting in a matrix of (species * sites)
    # ensuring that density dependence is applied per site
    k = max(P - sum(X, 1), 0.0);

    # Total cover of small massives and encrusting
    X_sm = sum(X[26:28, :]);

    # Total cover of largest tabular Acropora
    X_tabular = (X[6, :] + X[12, :]); # this is for enhanced and unenhanced

    k_X_r = k .* X .* r;
    k_rec = k .* rec;
    X_mb = X .* mb;

    # Tabular Acropora Enhanced
    kX5 = k .* X[5, :]
    Y[5, :] = k_X_r[4, :] - kX5 .* (r[5] + comp .* X_sm) - X_mb[5, :];
    Y[6, :] = kX5 .* (r[5] + comp .* X_sm) + k_X_r[6, :] - X_mb[6, :];

    # Tabular Acropora Unenhanced
    kX11 = k .* X[11, :]
    Y[11, :] = k_X_r[10, :] - kX11 .* (r[11] + comp .* X_sm) - X_mb[11, :];
    Y[12, :] = kX11 .* (r[11] + comp .* X_sm) + k_X_r[12, :] - X_mb[12, :];

    # Corymbose Acropora Enhanced
    Y[13, :] = k_rec[3, :] - k .* X[13, :] .* r[13] - X_mb[13,:];

    # Small massives and encrusting Unenhanced
    Y[25, :] = k_rec[5, :] - k .* X[25, :] .* r[25] - X_mb[25,:];
    Y[26:28, :] = k_X_r[25:27, :] - k_X_r[26:28, :] - X[26:28, :] .* (mb[26] + comp .* X_tabular);

    # Small size classes
    Y[[1,7,19,31], :] = k_rec[[1,2,4,6], :] - k_X_r[[1,7,19,31], :] - X_mb[[1,7,19,31], :];

    # Mid size classes
    Y[[2:4,8:10,14:17,20:23,29,32:35], :] = k_X_r[[1:3,7:9,13:16,19:22,28,31:34], :] - k_X_r[[2:4,8:10,14:17,20:23,29,32:35], :] - X_mb[[2:4,8:10,14:17,20:23,29,32:35], :];

    # Larger size classes
    Y[[18,24,30,36], :] = k_X_r[[17,23,29,35], :] + k_X_r[[18,24,30,36], :] - X_mb[[18,24,30,36], :];

    # Ensure no non-negative values
    Y = max(Y, 0);
end
