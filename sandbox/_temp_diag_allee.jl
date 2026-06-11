using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using ADRIA, SparseArrays, Statistics

dom = ADRIA.load_domain(ADRIA.RMEDomain, joinpath(@__DIR__, "data", "rme_ml_2025_06_05"), "45")
conn = sparse(dom.conn.data)
n = size(conn, 1)
diag_vals = [conn[i,i] for i in 1:n]
println("GBR Self-connectivity (diagonal):")
println("  min=$(minimum(diag_vals)), max=$(maximum(diag_vals))")
println("  mean=$(mean(diag_vals)), median=$(median(diag_vals))")
println("  zeros=$(sum(diag_vals .== 0.0)) out of $n")
println("  >0.01: $(sum(diag_vals .> 0.01))")
println("  >0.1: $(sum(diag_vals .> 0.1))")

row_sums = vec(sum(conn; dims=2))
retention = diag_vals ./ max.(row_sums, 1e-12)
println("\nSelf-retention fraction (diag / row_sum):")
println("  mean=$(mean(retention)), median=$(median(retention))")

println("\nAllee effect N^2/(A^2 + N^2) with A=1.0:")
for d in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0]
    allee = d^2 / (1.0^2 + d^2)
    println("  N=$d: allee=$(round(allee, digits=4))")
end
println("\nAllee effect with A=0.1:")
for d in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0]
    allee = d^2 / (0.1^2 + d^2)
    println("  N=$d: allee=$(round(allee, digits=4))")
end

# What the dispersal ACTUALLY does to self-recruitment:
# cots_timestep! sets N[1] = R + IMM  (say ~0.97)
# disperse_cots_larvae! REPLACES N[1] with sum(source * conn)
# For self: N[1] becomes R * conn[i,i]
# So most of the Ricker recruitment is DESTROYED
println("\nWith Ricker R=0.97, what remains at source after dispersal?")
for d in [0.0, 0.001, 0.01, 0.1, 0.28, 0.32]
    remaining = 0.97 * d
    println("  conn[i,i]=$d: remaining=$(round(remaining, digits=4))")
end
