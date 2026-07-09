using Pkg; Pkg.activate(".")
using ADRIA
using CSV, DataFrames, Statistics

# 1. Load domain and mapping
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "historical")

site_to_reef = CSV.read("sandbox/data/Lizard_Historical_v0.1/site_to_reef.csv", DataFrame)
site_to_reef.reef_name_clean = [split(r, " (")[1] for r in site_to_reef.reef_name]
unique_reefs_clean = unique(site_to_reef.reef_name_clean)

# 2. Load sweep results
res_df = CSV.read("sandbox/data/formal_calibration_ensemble.csv", DataFrame)
sort!(res_df, :loss)
top5 = res_df[1:5, :]

println("Top 5 scenarios:")
println(top5)

# 3. Re-run top 5 and extract trajectories
p_df = ADRIA.param_table(dom)
sim_df = DataFrame(run_id=Int[], year=Int[], reef_name=String[], sim_cots_norm=Float64[], sim_coral=Float64[])

for i in 1:5
    p = p_df[1, :]
    p.a_F = top5.a_F[i]
    p.a_S = top5.a_S[i]
    p.IMM = top5.IMM[i]
    ENV["COTS_INITIAL_MULTIPLIER"] = string(top5.seed_mult[i])
    
    # Enable external larval pulse (Cairns Initiation Box)
    ENV["COTS_EXTERNAL_PULSE"] = "true"
    ENV["COTS_PULSE_PERIOD"] = "15"
    ENV["COTS_PULSE_OFFSET"] = "1"
    ENV["COTS_PULSE_VAL"] = "1.5"
    
    rs = ADRIA.run_scenario(dom, p)
    adult_cots_site = rs.cots_log[:, 3, :]
    total_cover_site = dropdims(sum(rs.raw, dims=(2, 3)), dims=(2, 3))
    
    for reef in unique_reefs_clean
        site_indices = findall(site_to_reef.reef_name_clean .== reef)
        if isempty(site_indices) continue end
        
        reef_sim_cots = [mean(adult_cots_site[t, site_indices]) for t in 1:40]
        max_sim = maximum(reef_sim_cots)
        reef_sim_cots_norm = max_sim > 0 ? reef_sim_cots ./ max_sim : zeros(40)
        
        reef_sim_coral = [mean(total_cover_site[t, site_indices]) for t in 1:40]
        
        for t in 1:40
            push!(sim_df, (top5.run_id[i], 1984 + t, reef, reef_sim_cots_norm[t], reef_sim_coral[t]))
        end
    end
end

CSV.write("sandbox/data/top5_trajectories.csv", sim_df)
println("Top 5 trajectories saved to sandbox/data/top5_trajectories.csv")
