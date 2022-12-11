import ADRIA.sensitivity: pawn, relative_importance


function relative_sensitivities(X, y; S=10, stat=:median)::Vector{Float64}
    return relative_importance(pawn(X, Array(y); S=S)[:, stat])
end

