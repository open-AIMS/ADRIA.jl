using Pkg; Pkg.activate(".")
using ADRIA
using CSV, DataFrames, Statistics
using LatinHypercubeSampling

# 1. Load domain
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "45")

# 2. Setup mapping
site_to_reef = CSV.read("sandbox/data/Lizard_Historical_v0.1/site_to_reef.csv", DataFrame)
clean_reef_names = [split(r, " (")[1] for r in site_to_reef.reef_name]
site_to_reef.reef_name_clean = clean_reef_names
unique_reefs_clean = unique(clean_reef_names)

df_cots_raw = CSV.read("sandbox/data/reef_cots.csv", DataFrame)
df_cots = filter(row -> row.reef_name in unique_reefs_clean && row.year >= 1985, df_cots_raw)

df_cots[!, :cotsptow_norm] .= 0.0
for reef in unique(df_cots.reef_name)
    idx = df_cots.reef_name .== reef
    max_val = maximum(df_cots[idx, :cotsptow])
    if max_val > 0
        df_cots[idx, :cotsptow_norm] = df_cots[idx, :cotsptow] ./ max_val
    end
end

# 3. Create LHS Ensemble
N_samples = 250
# Bounds: a_F (0.1, 2.0), a_S (0.01, 0.9), IMM (0.0, 0.1), initial_seed_multiplier (0.5, 3.0)
lhs_plan = scaleLHC(randomLHC(N_samples, 4), [(0.1, 2.0), (0.01, 0.9), (0.0, 0.1), (0.5, 3.0)])

println("Starting sequential LHS sweep of $N_samples scenarios...")

p_df = ADRIA.param_table(dom)
results_df = DataFrame(run_id=Int[], a_F=Float64[], a_S=Float64[], IMM=Float64[], seed_mult=Float64[], loss=Float64[])
losses = zeros(Float64, N_samples)

for i in 1:N_samples
    println("Running scenario $i / $N_samples ...")
    p = p_df[1, :]
    p.a_F = lhs_plan[i, 1]
    p.a_S = lhs_plan[i, 2]
    p.IMM = lhs_plan[i, 3]
    
    # Inject initial_seed_multiplier via ENV variable
    ENV["COTS_INITIAL_MULTIPLIER"] = string(lhs_plan[i, 4])
    
    rs = ADRIA.run_scenario(dom, p)
    adult_cots_site = rs.cots_log[:, 3, :] 
    
    loss = 0.0
    valid_reefs = 0
    
    for reef in unique_reefs_clean
        site_indices = findall(site_to_reef.reef_name_clean .== reef)
        if isempty(site_indices) continue end
        
        reef_sim_cots = [mean(adult_cots_site[t, site_indices]) for t in 1:40]
        max_sim = maximum(reef_sim_cots)
        reef_sim_cots_norm = max_sim > 0 ? reef_sim_cots ./ max_sim : zeros(40)
        
        reef_emp = filter(row -> row.reef_name == reef, df_cots)
        if nrow(reef_emp) > 0
            valid_reefs += 1
            
            sim_peak_idx = argmax(reef_sim_cots_norm)
            sim_peak_year = 1984 + sim_peak_idx
            
            emp_max_idx = argmax(reef_emp.cotsptow_norm)
            emp_peak_year = reef_emp.year[emp_max_idx]
            
            phase_penalty = abs(sim_peak_year - emp_peak_year) * 2.0
            
            sim_peak2_idx = 20 + argmax(reef_sim_cots_norm[21:40])
            sim_peak2_year = 1984 + sim_peak2_idx
            
            emp_after_2005 = filter(row -> row.year > 2005, reef_emp)
            if nrow(emp_after_2005) > 0
                emp_max2_idx = argmax(emp_after_2005.cotsptow_norm)
                emp_peak2_year = emp_after_2005.year[emp_max2_idx]
                phase_penalty += abs(sim_peak2_year - emp_peak2_year) * 2.0
            end
            
            loss += phase_penalty
            
            for row in eachrow(reef_emp)
                t_idx = row.year - 1984
                if 1 <= t_idx <= 40
                    loss += (reef_sim_cots_norm[t_idx] - row.cotsptow_norm)^2
                end
            end
        end
    end
    
    losses[i] = valid_reefs > 0 ? loss / valid_reefs : 9999.0
end

for i in 1:N_samples
    push!(results_df, (i, lhs_plan[i, 1], lhs_plan[i, 2], lhs_plan[i, 3], lhs_plan[i, 4], losses[i]))
end

CSV.write("sandbox/data/formal_calibration_ensemble.csv", results_df)
println("Sweep complete! Results saved to sandbox/data/formal_calibration_ensemble.csv")
