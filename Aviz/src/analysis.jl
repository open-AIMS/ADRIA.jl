import ADRIA.sensitivity: pawn, relative_importance

function relative_sensitivities(X, y; S=10, stat=:median)::Vector{Float64}
    med_pawn_idx = pawn(X, Array(y); S=S)[:, stat]
    rel_pawn_idx = relative_importance(med_pawn_idx)

    return vec(rel_pawn_idx)
end

