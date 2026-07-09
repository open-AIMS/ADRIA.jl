using Pkg; Pkg.activate(".")
using ADRIA
using CSV, DataFrames, Statistics

println("Loading Lizard Island historical domain...")
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "historical")
site_to_reef = CSV.read("sandbox/data/Lizard_Historical_v0.1/site_to_reef.csv", DataFrame)

clean_reef_names = [split(r, " (")[1] for r in site_to_reef.reef_name]
site_to_reef.reef_name_clean = clean_reef_names
unique_reefs_clean = unique(clean_reef_names)

# Best parameter candidate from Phase 9b (run_id = 2, loss = 23.0897):
best_a_F = 1.34378
best_a_S = 0.0957831
best_IMM = 0.0947791
best_seed_mult = 2.23695

# Set up environment variables for Cairns Initiation Box pulse & initial density
ENV["COTS_INITIAL_MULTIPLIER"] = string(best_seed_mult)
ENV["COTS_EXTERNAL_PULSE"] = "true"
ENV["COTS_PULSE_PERIOD"] = "15"
ENV["COTS_PULSE_OFFSET"] = "1"
ENV["COTS_PULSE_VAL"] = "1.5"
# The pulse code injects into the explicit COTS seed locations. Use the first
# sites as a simple upstream initiation box for this diagnostic workflow.
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = "10"

# Run 16 stochastic simulation runs of our single best parameter candidate (power of 2 required by Sobol)
N_scens = 16
println("Sampling $N_scens stochastic scenarios...")
scens = ADRIA.sample(dom, N_scens)

p_df = ADRIA.param_table(dom)
for col in names(scens)
    if !(col in ["dhw_scenario", "wave_scenario", "cyclone_mortality_scenario"]) && !startswith(col, "surv_") && !in(col, ["fecundity", "coral_recruitment"])
        scens[!, col] .= p_df[1, col]
    end
end

# Lock COTS parameters across all scenarios to our single best candidate
scens.a_F .= best_a_F
scens.a_S .= best_a_S
scens.IMM .= best_IMM

println("Running $N_scens stochastic simulation runs sequentially...")
sim_df = DataFrame(
    sim_id = Int[],
    year = Int[],
    reef_name = String[],
    sim_cots_adult = Float64[],
    sim_cots_norm = Float64[],
    sim_coral = Float64[]
)
site_df = DataFrame(
    sim_id = Int[],
    year = Int[],
    reef_name = String[],
    site_index = Int[],
    sim_cots_adult = Float64[],
    sim_coral = Float64[]
)

for s in 1:N_scens
    println("Simulating stochastic run $s / $N_scens ...")
    p = scens[s, :]
    rs = ADRIA.run_scenario(dom, p)
    adult_cots_site = rs.cots_log[:, 3, :]
    total_cover_site = dropdims(sum(rs.raw, dims=(2, 3)), dims=(2, 3))
    
    for reef in unique_reefs_clean
        site_indices = findall(site_to_reef.reef_name_clean .== reef)
        if isempty(site_indices) continue end
        
        reef_sim_cots = [mean(adult_cots_site[t, site_indices]) for t in 1:40]
        reef_sim_coral = [mean(total_cover_site[t, site_indices]) for t in 1:40]
        
        for t in 1:40
            push!(sim_df, (s, 1984 + t, reef, reef_sim_cots[t], 0.0, reef_sim_coral[t]))
            for site_idx in site_indices
                push!(
                    site_df,
                    (
                        s,
                        1984 + t,
                        reef,
                        site_idx,
                        adult_cots_site[t, site_idx],
                        total_cover_site[t, site_idx]
                    )
                )
            end
        end
    end
end

# Normalize COTS with one shared maximum per reef across all stochastic runs.
for reef in unique(sim_df.reef_name)
    reef_rows = sim_df.reef_name .== reef
    max_sim = maximum(sim_df[reef_rows, :sim_cots_adult])
    if max_sim > 0
        sim_df[reef_rows, :sim_cots_norm] .= sim_df[reef_rows, :sim_cots_adult] ./ max_sim
    end
end

CSV.write("sandbox/data/best_stochastic_trajectories.csv", sim_df)
CSV.write("sandbox/data/best_stochastic_site_trajectories.csv", site_df)
println("Saved stochastic trajectories to sandbox/data/best_stochastic_trajectories.csv")
println("Saved site-level stochastic trajectories to sandbox/data/best_stochastic_site_trajectories.csv")
