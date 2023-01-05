import ADRIA.analysis: normalize
import ADRIA.sensitivity: pawn


function relative_sensitivities(X, y; S=10, stat=:median)::Vector{Float64}
    return normalize(pawn(X, Array(y); S=S)[:, stat])
end

