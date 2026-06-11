using Pkg; Pkg.activate(".")
using ADRIA
using CSV, DataFrames, Statistics

# Best candidate from optimization: [a_F, a_S, IMM]
best_params = [1.7950474114611144, 0.09580041820105212, 0.08860289578720529]

dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "historical")
site_to_reef = CSV.read("sandbox/data/Lizard_Historical_v0.1/site_to_reef.csv", DataFrame)

# Clean mapping reef names by removing the " (XX-XXX)" part
clean_reef_names = [split(r, " (")[1] for r in site_to_reef.reef_name]
site_to_reef.reef_name_clean = clean_reef_names
unique_reefs_clean = unique(clean_reef_names)

df_cots_raw = CSV.read("sandbox/data/reef_cots.csv", DataFrame)
df_cots = filter(row -> row.reef_name in unique_reefs_clean && row.year >= 1985, df_cots_raw)

# Normalize empirical
df_cots[!, :cotsptow_norm] .= 0.0
for reef in unique(df_cots.reef_name)
    idx = df_cots.reef_name .== reef
    max_val = maximum(df_cots[idx, :cotsptow])
    if max_val > 0
        df_cots[idx, :cotsptow_norm] = df_cots[idx, :cotsptow] ./ max_val
    end
end

p_df = ADRIA.param_table(dom)
p = p_df[1, :]
p.a_F = best_params[1]
p.a_S = best_params[2]
p.IMM = best_params[3]

rs = ADRIA.run_scenario(dom, p)
adult_cots_site = rs.cots_log[:, 3, :] 
# Extract total absolute coral cover: sum over size and species dimensions (dims 2 and 3)
total_cover_site = dropdims(sum(rs.raw, dims=(2, 3)), dims=(2, 3))

# Save simulated data
sim_df = DataFrame(year = Int[], reef_name = String[], sim_cots_norm = Float64[], sim_coral = Float64[])

for reef in unique_reefs_clean
    site_indices = findall(site_to_reef.reef_name_clean .== reef)
    if isempty(site_indices) continue end
    
    reef_sim_cots = [mean(adult_cots_site[t, site_indices]) for t in 1:40]
    max_sim = maximum(reef_sim_cots)
    if max_sim > 0
        reef_sim_cots_norm = reef_sim_cots ./ max_sim
    else
        reef_sim_cots_norm = zeros(40)
    end
    
    # Calculate mean coral cover across sites for this reef
    reef_sim_coral = [mean(total_cover_site[t, site_indices]) for t in 1:40]
    
    for t in 1:40
        push!(sim_df, (1984 + t, reef, reef_sim_cots_norm[t], reef_sim_coral[t]))
    end
end

CSV.write("sandbox/data/calibration_results_sim.csv", sim_df)
CSV.write("sandbox/data/calibration_results_emp.csv", df_cots)

# Empirical Coral Cover
df_coral_raw = CSV.read("sandbox/data/reef_manta.csv", DataFrame)
df_coral = filter(row -> row.reef_name in unique_reefs_clean && row.report_year >= 1985, df_coral_raw)
CSV.write("sandbox/data/calibration_results_emp_coral.csv", df_coral)

println("Results saved to CSV!")
