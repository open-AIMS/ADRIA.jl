using Pkg; Pkg.activate(".")
using ADRIA
using CSV, DataFrames, Statistics
using BlackBoxOptim

# 1. Load historical domain
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "historical")

# 2. Load Mapping and Empirical Data
site_to_reef = CSV.read("sandbox/data/Lizard_Historical_v0.1/site_to_reef.csv", DataFrame)

# Clean mapping reef names by removing the " (XX-XXX)" part
clean_reef_names = [split(r, " (")[1] for r in site_to_reef.reef_name]
site_to_reef.reef_name_clean = clean_reef_names

unique_reefs_clean = unique(clean_reef_names)

df_cots_raw = CSV.read("sandbox/data/reef_cots.csv", DataFrame)
df_cots = filter(row -> row.reef_name in unique_reefs_clean && row.year >= 1985, df_cots_raw)

# We normalize empirical COTS per reef to [0,1]
df_cots[!, :cotsptow_norm] .= 0.0
for reef in unique(df_cots.reef_name)
    idx = df_cots.reef_name .== reef
    max_val = maximum(df_cots[idx, :cotsptow])
    if max_val > 0
        df_cots[idx, :cotsptow_norm] = df_cots[idx, :cotsptow] ./ max_val
    end
end

# 3. Create Loss Function
function run_and_evaluate(params::Vector{Float64})
    # params: [a_F, a_S, IMM]
    p_df = ADRIA.param_table(dom)
    p = p_df[1, :]
    
    # Update COTS parameters
    p.a_F = params[1]
    p.a_S = params[2]
    p.IMM = params[3]
    # Minor tuning for coral
    # p.coral_growth_rate_modifier = params[4] ... etc if needed
    
    rs = ADRIA.run_scenario(dom, p)
    
    # Extract adult COTS over time (shape: time, life_stages, sites) -> adults are index 3
    # rs.cots_log is (31, 3, 2914)
    adult_cots_site = rs.cots_log[:, 3, :] 
    
    loss = 0.0
    valid_reefs = 0
    
    for reef in unique_reefs_clean
        # find sites belonging to this reef
        site_indices = findall(site_to_reef.reef_name_clean .== reef)
        if isempty(site_indices) continue end
        
        # Calculate mean adult COTS density for this reef over time (1:40)
        # We assume site_indices directly map to 1:2914 locations in ADRIA
        reef_sim_cots = [mean(adult_cots_site[t, site_indices]) for t in 1:40]
        
        # Normalize simulated COTS
        max_sim = maximum(reef_sim_cots)
        if max_sim > 0
            reef_sim_cots_norm = reef_sim_cots ./ max_sim
        else
            reef_sim_cots_norm = zeros(40)
        end
        
        # Compare with empirical
        reef_emp = filter(row -> row.reef_name == reef, df_cots)
        if nrow(reef_emp) > 0
            valid_reefs += 1
            
            # Phase Shift Penalty: force peak timing alignment
            sim_peak_idx = argmax(reef_sim_cots_norm)
            sim_peak_year = 1984 + sim_peak_idx
            
            emp_max_idx = argmax(reef_emp.cotsptow_norm)
            emp_peak_year = reef_emp.year[emp_max_idx]
            
            phase_penalty = abs(sim_peak_year - emp_peak_year) * 2.0
            loss += phase_penalty
            
            for row in eachrow(reef_emp)
                t_idx = row.year - 1984 # 1985 -> 1
                if 1 <= t_idx <= 40
                    sim_val = reef_sim_cots_norm[t_idx]
                    emp_val = row.cotsptow_norm
                    loss += (sim_val - emp_val)^2
                end
            end
        end
    end
    
    return valid_reefs > 0 ? loss / valid_reefs : 9999.0
end

println("Loss function ready! Running optimization...")

# Run optimization
# Parameters: [a_F, a_S, IMM]
# Bounds:
# a_F: 0.1 to 2.0 (starvation)
# a_S: 0.01 to 0.9 (survival)
# IMM: 0.0 to 0.1 (immigration)
res = bboptimize(run_and_evaluate; SearchSpace=[(0.1, 2.0), (0.01, 0.9), (0.0, 0.1)], NumDimensions=3, MaxSteps=50)
println("Best candidate: ", best_candidate(res))
println("Best fitness: ", best_fitness(res))
