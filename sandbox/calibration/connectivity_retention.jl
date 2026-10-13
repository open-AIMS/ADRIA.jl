"""
Source-conserving local-retention transform for reef connectivity matrices
(crash_mechanism_protocol.md, 2026-10-13: Lizard Island local-retention screen).

Orientation: rows are SOURCE reefs, columns are SINK reefs; entries are the fraction of
source production arriving at the sink (pelagic mortality included), so a row sum is
the source's total larval survival.

For each source `s` in a retention group `G`, the share of its surviving larvae that
stays in `G` is set to `rho`. Within `G` the original pattern is kept; exports outside
`G` are scaled proportionally. Row sums (source survival) are preserved, so larvae are
redistributed, never created.
"""
function retain_within_groups(matrix::AbstractMatrix{<:Real}, groups::Vector{Vector{Int}},
    rho::Float64)::Matrix{Float64}
    0.0 <= rho <= 1.0 || throw(ArgumentError("rho must lie in [0, 1]"))
    out = Matrix{Float64}(matrix)
    n = size(out, 1)
    size(out, 2) == n || throw(DimensionMismatch("connectivity must be square"))
    all(>=(0.0), out) || throw(ArgumentError("connectivity has negative entries"))
    seen = falses(n)
    for group in groups
        isempty(group) && throw(ArgumentError("empty retention group"))
        all(1 .<= group .<= n) || throw(ArgumentError("retention group index out of range"))
        any(seen[group]) && throw(ArgumentError("retention groups overlap"))
        seen[group] .= true
        outside = setdiff(1:n, group)
        for s in group
            total = sum(@view out[s, :])
            total > 0.0 || continue                      # zero rows (no source) stay zero
            kept = sum(@view out[s, group])
            exported = total - kept
            if kept > 0.0
                out[s, group] .*= rho * total / kept
            else                                          # no Owen path within G: keep on self
                out[s, s] = rho * total
            end
            if exported > 0.0
                out[s, outside] .*= (1.0 - rho) * total / exported
            elseif rho < 1.0
                throw(ArgumentError("source $s has no exports to reduce retention below 1"))
            end
        end
    end
    return out
end

"""Share of each source's surviving larvae that stays within its group (diagnostic)."""
function retained_share(matrix::AbstractMatrix{<:Real}, group::Vector{Int})::Vector{Float64}
    return [let total = sum(matrix[s, :]); total > 0 ? sum(matrix[s, group]) / total : NaN end
            for s in group]
end
