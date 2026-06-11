using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using ADRIA
using Statistics
using SparseArrays
using Printf
using StaticArrays

println("=" ^ 80)
println("COTS Dispersal Connectivity Diagnostic")
println("=" ^ 80)

# --- Load Moore Domain ---
println("\n--- Moore Domain ---")
MOORE_PATH = joinpath(@__DIR__, "data", "Moore_2025-08-22_v080")
dom_moore = ADRIA.load_domain(MOORE_PATH, "45")
conn_moore = sparse(dom_moore.conn.data)
n_moore = size(conn_moore, 1)

println("  Locations: $n_moore")
println("  Conn shape: $(size(conn_moore))")
println("  Conn nnz: $(nnz(conn_moore))")
println("  Conn density: $(nnz(conn_moore) / (n_moore^2))")
println("  Conn value range: [$(minimum(nonzeros(conn_moore))), $(maximum(nonzeros(conn_moore)))]")
println("  Conn mean nonzero: $(mean(nonzeros(conn_moore)))")
println("  Conn median nonzero: $(median(nonzeros(conn_moore)))")

# Row sums (total outgoing from each source)
row_sums_moore = vec(sum(conn_moore; dims=2))
println("  Row sums (source -> all sinks):")
println("    min=$(minimum(row_sums_moore)), max=$(maximum(row_sums_moore)), mean=$(mean(row_sums_moore))")

# Col sums (total incoming to each sink)
col_sums_moore = vec(sum(conn_moore; dims=1))
println("  Col sums (all sources -> sink):")
println("    min=$(minimum(col_sums_moore)), max=$(maximum(col_sums_moore)), mean=$(mean(col_sums_moore))")

# How many non-zero entries per row (how many sinks does each source connect to?)
nnz_per_row_moore = [length(findnz(conn_moore[i, :])[1]) for i in 1:n_moore]
println("  Non-zero entries per row (outgoing connections):")
println("    min=$(minimum(nnz_per_row_moore)), max=$(maximum(nnz_per_row_moore)), mean=$(mean(nnz_per_row_moore))")

# How many non-zero entries per col (how many sources feed into each sink?)
nnz_per_col_moore = diff(conn_moore.colptr)
println("  Non-zero entries per col (incoming connections):")
println("    min=$(minimum(nnz_per_col_moore)), max=$(maximum(nnz_per_col_moore)), mean=$(mean(nnz_per_col_moore))")

# --- Load GBR Domain ---
println("\n--- GBR Domain (RMEDomain) ---")
GBR_PATH = joinpath(@__DIR__, "data", "rme_ml_2025_06_05")
dom_gbr = ADRIA.load_domain(ADRIA.RMEDomain, GBR_PATH, "45")
conn_gbr = sparse(dom_gbr.conn.data)
n_gbr = size(conn_gbr, 1)

println("  Locations: $n_gbr")
println("  Conn shape: $(size(conn_gbr))")
println("  Conn nnz: $(nnz(conn_gbr))")
println("  Conn density: $(nnz(conn_gbr) / (n_gbr^2))")
println("  Conn value range: [$(minimum(nonzeros(conn_gbr))), $(maximum(nonzeros(conn_gbr)))]")
println("  Conn mean nonzero: $(mean(nonzeros(conn_gbr)))")
println("  Conn median nonzero: $(median(nonzeros(conn_gbr)))")

# Row sums
row_sums_gbr = vec(sum(conn_gbr; dims=2))
println("  Row sums (source -> all sinks):")
println("    min=$(minimum(row_sums_gbr)), max=$(maximum(row_sums_gbr)), mean=$(mean(row_sums_gbr))")

# Col sums
col_sums_gbr = vec(sum(conn_gbr; dims=1))
println("  Col sums (all sources -> sink):")
println("    min=$(minimum(col_sums_gbr)), max=$(maximum(col_sums_gbr)), mean=$(mean(col_sums_gbr))")

# NNZ per row/col
nnz_per_col_gbr = diff(conn_gbr.colptr)
println("  Non-zero entries per col (incoming connections):")
println("    min=$(minimum(nnz_per_col_gbr)), max=$(maximum(nnz_per_col_gbr)), mean=$(mean(nnz_per_col_gbr))")

# --- COTS Larval Production Diagnostic ---
println("\n--- COTS Larval Production Test ---")
println("Simulating one timestep of COTS recruitment at an outbreak reef...")

# Using the default COTS params from scenario.jl
cots_params = (
    a_F = 0.3, a_S = 0.05, IMM = 0.001, p_tilde = 0.9, C_max = 0.6,
    m1 = 0.6, m2 = 0.2, m3 = 0.1,
    a_ricker = 4.0, b_ricker = 0.8,
    allee_threshold = 1.0,
    imm_threshold = 0.35, eta_imm = 2.0,
    tau_condition = 3.0,
    h = 0.0, eta_F = 1.0, eta_S = 1.0
)

# Simulate a single reef with adult COTS = 2.0 and healthy coral
N_init = MVector{3, Float64}(0.02, 0.02, 2.0)
model = ADRIA.CotsHuman(N_init, 0.8, cots_params)
F = 0.3  # healthy fast coral cover
S = 0.15 # healthy slow coral cover

println("  Before step: N = $(model.N)")
Cons_F, Cons_S = ADRIA.cots_timestep!(model, F, S)
println("  After step:  N = $(model.N)")
println("  Recruits (N[1], age-0) = $(model.N[1])")
println("  Consumption: F=$(Cons_F), S=$(Cons_S)")

# Now check: if this reef produces N[1] recruits, how much reaches neighbor sinks?
recruit_output = model.N[1]

# For Moore domain: pick a well-connected source reef
println("\n--- Dispersal Propagation Test ---")
# Find the 14- prefix reefs in GBR 
seed_locs_gbr = findall(x -> startswith(String(x), "14-"), dom_gbr.loc_data.GBRMPA_ID)
println("  GBR 14- prefix locations: $(length(seed_locs_gbr))")

# Check row sums for these source reefs (total outgoing connectivity)
if !isempty(seed_locs_gbr)
    src_row_sums = row_sums_gbr[seed_locs_gbr]
    println("  Row sums for 14- reefs (total outgoing probability):")
    println("    min=$(minimum(src_row_sums)), max=$(maximum(src_row_sums)), mean=$(mean(src_row_sums))")
    
    # How many non-14- sinks do these sources connect to?
    non_seed = setdiff(1:n_gbr, seed_locs_gbr)
    local total_outgoing_to_others = 0.0
    for src in seed_locs_gbr
        for sink in non_seed
            val = conn_gbr[src, sink]
            if val > 0
                total_outgoing_to_others += val
            end
        end
    end
    avg_outgoing_to_others = total_outgoing_to_others / length(seed_locs_gbr)
    println("  Average total connectivity from each 14- reef to non-14- sinks: $avg_outgoing_to_others")
    
    # Expected incoming larvae at each non-seed sink after one dispersal step
    # with all 14- reefs producing `recruit_output` recruits
    println("\n  If each 14- reef produces $(recruit_output) recruits:")
    println("  Expected incoming at each non-14- sink = recruit_output * conn[src,sink]")
    
    # Find the top 10 non-seed sinks that would receive the most
    incoming_at_non_seed = zeros(n_gbr)
    for sink in 1:n_gbr
        if sink in seed_locs_gbr
            continue
        end
        for src in seed_locs_gbr
            incoming_at_non_seed[sink] += recruit_output * conn_gbr[src, sink]
        end
    end
    
    non_zero_sinks = findall(incoming_at_non_seed .> 0)
    println("  Non-seed sinks receiving ANY larvae: $(length(non_zero_sinks)) out of $(length(non_seed))")
    
    if !isempty(non_zero_sinks)
        sorted_idx = sortperm(incoming_at_non_seed[non_zero_sinks], rev=true)
        top_n = min(10, length(sorted_idx))
        println("  Top $top_n receiving sinks:")
        for i in 1:top_n
            sink = non_zero_sinks[sorted_idx[i]]
            id = dom_gbr.loc_data.GBRMPA_ID[sink]
            println("    Sink $sink ($id): incoming = $(incoming_at_non_seed[sink])")
        end
        println("  Max incoming: $(maximum(incoming_at_non_seed))")
        println("  Mean incoming (non-zero only): $(mean(incoming_at_non_seed[non_zero_sinks]))")
    end
end

# For comparison, check Moore domain dispersal
println("\n--- Moore Domain Dispersal Comparison ---")
# Pick source reef 1 in Moore
src_moore = 1
row_sum_moore_src = row_sums_moore[src_moore]
println("  Moore reef 1 row sum: $row_sum_moore_src")
println("  If Moore reef 1 produces $(recruit_output) recruits:")
outgoing_moore = [recruit_output * conn_moore[src_moore, sink] for sink in 1:n_moore]
non_zero_moore = findall(outgoing_moore .> 0)
println("  Non-zero recipient sinks: $(length(non_zero_moore))")
if !isempty(non_zero_moore)
    println("  Max incoming: $(maximum(outgoing_moore))")
    println("  Mean incoming (non-zero): $(mean(outgoing_moore[non_zero_moore]))")
end

println("\n" * "=" ^ 80)
println("DIAGNOSTIC COMPLETE")
println("=" ^ 80)
