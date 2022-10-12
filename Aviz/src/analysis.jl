function relative_sensitivities(X, y; S=10, stat=:median)
    med_pawn_idx = ADRIA.sensitivity.pawn(X, Array(y); S=S)[:, stat]
    rel_pawn_idx = ADRIA.sensitivity.relative_importance(med_pawn_idx)
    # rel_pawn_idx = rel_pawn_idx[interv_idx]
    return vec(rel_pawn_idx)
end

