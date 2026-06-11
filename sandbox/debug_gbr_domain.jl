using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using ADRIA
using Plots
using Statistics
using DataFrames

println("Loading whole of GBR domain (RMEDomain)...")
DOMAIN_PATH = joinpath(@__DIR__, "data", "rme_ml_2025_06_05")
# Explicitly load as RMEDomain, which is how ADRIA handles ReefMod Engine data
dom = ADRIA.load_domain(ADRIA.RMEDomain, DOMAIN_PATH, "45")
println("Loaded domain successfully!")

# Define our 4 key reefs by GBRMPA ID
target_ids = ["14-116", "18-075", "16-029", "21-566"]
target_names = ["Lizard Island", "John Brewer Reef", "Batt Reef", "Gannet Cay"]

# Find location indices for the target reefs
key_reef_locs = Dict{String, Vector{Int}}()
for (name, id_prefix) in zip(target_names, target_ids)
    # Some reefs like Lizard Island have multiple sub-polygons (e.g. 14-116a, 14-116b)
    locs = findall(x -> startswith(String(x), id_prefix), dom.loc_data.GBRMPA_ID)
    if !isempty(locs)
        key_reef_locs[name] = locs
    else
        println("Warning: Target reef $name ($id_prefix) not found in domain.")
    end
end

# The user requested to seed COTS at ALL reefs with the "14-" prefix
# to observe if they spread via the connectivity network.
seed_locs = findall(x -> startswith(String(x), "14-"), dom.loc_data.GBRMPA_ID)

println("Seeding COTS at key reefs...")
ENV["ADRIA_COTS_ENABLED"] = "true"
ENV["COTS_EXTERNAL_PULSE"] = "true"
ENV["ADRIA_DEBUG_INIT_DENSITY"] = "2.0"
ENV["ADRIA_DEBUG_SEED_LOCATIONS"] = join(seed_locs, ",")
ENV["COTS_IMMIGRATION_SCALAR"] = "100.0"  # Amplify cross-reef larval immigration

tf = length(ADRIA.timesteps(dom))
scens = ADRIA.sample(dom, 2)
scens.guided .= 0
scens.seed_strategy .= 0
scens.fog_strategy .= 0
scens.mc_strategy .= 0
scens.N_seed_TA .= 0.0
scens.N_seed_CA .= 0.0
scens.N_seed_CNA .= 0.0
scens.N_seed_SM .= 0.0
scens.N_seed_LM .= 0.0
scens.N_mc_settlers .= 0.0
scens.fogging .= 0.0

println("Running baseline scenario (No COTS)...")
ENV["ADRIA_COTS_ENABLED"] = "false"
res_no_cots = ADRIA.run_scenario(dom, scens[1, :])

println("Running scenario with COTS enabled...")
ENV["ADRIA_COTS_ENABLED"] = "true"
res_cots = ADRIA.run_scenario(dom, scens[1, :])

# --- Extract Data ---
raw_data_no_cots = hasproperty(res_no_cots.raw, :data) ? res_no_cots.raw.data : res_no_cots.raw
total_coral_no_cots = dropdims(sum(raw_data_no_cots, dims=(2, 3)), dims=(2, 3))

raw_data_cots = hasproperty(res_cots.raw, :data) ? res_cots.raw.data : res_cots.raw
total_coral_cots = dropdims(sum(raw_data_cots, dims=(2, 3)), dims=(2, 3))

cots_data = hasproperty(res_cots.cots_log, :data) ? res_cots.cots_log.data : res_cots.cots_log
adult_cots = cots_data[:, 3, :]

# --- Plot 1: Key Reefs Comparison ---
println("Generating Key Reefs plot...")
key_plot_list = []
for name in target_names
    locs = get(key_reef_locs, name, Int[])
    if isempty(locs)
        continue
    end
    
    reef_coral_no_cots = vec(mean(total_coral_no_cots[:, locs], dims=2))
    reef_coral_cots = vec(mean(total_coral_cots[:, locs], dims=2))
    
    max_val = max(maximum(reef_coral_no_cots), maximum(reef_coral_cots)) * 100
    p = plot(title="$name (w/ vs w/o COTS)", ylabel="Coral %", legend=:topleft, ylims=(0, max_val*1.2))
    
    for loc in locs
        plot!(p, 1:tf, total_coral_no_cots[:, loc] .* 100, color=:green, alpha=0.15, label="")
        plot!(p, 1:tf, total_coral_cots[:, loc] .* 100, color=:blue, alpha=0.15, label="")
    end
    
    plot!(p, 1:tf, reef_coral_no_cots .* 100, label="No COTS Mean", color=:green, linewidth=3)
    plot!(p, 1:tf, reef_coral_cots .* 100, label="With COTS Mean", color=:blue, linewidth=3)
    
    push!(key_plot_list, p)
end

cols = 2
rows = ceil(Int, length(key_plot_list) / cols)
key_final_plot = plot(key_plot_list..., layout=(rows, cols), size=(1200, rows * 350), margin=5Plots.mm)
savefig(key_final_plot, joinpath(@__DIR__, "gbr_key_reefs_comparison.png"))
println("Saved Key Reefs plot to gbr_key_reefs_comparison.png")

# --- Plot 1.5: Key Reefs Coral vs COTS (With COTS Scenario) ---
println("Generating Key Reefs Coral vs COTS plot...")
cots_plot_list = []
for name in target_names
    locs = get(key_reef_locs, name, Int[])
    if isempty(locs)
        continue
    end
    
    reef_coral = vec(mean(total_coral_cots[:, locs], dims=2))
    reef_cots = vec(mean(adult_cots[:, locs], dims=2))
    
    p2 = plot(title="$name (Coral vs COTS)", ylabel="Percentage / Scaled Density", legend=:topleft, ylims=(0, maximum(reef_coral)*120))
    
    for loc in locs
        loc_coral = total_coral_cots[:, loc] .* 100
        loc_cots = adult_cots[:, loc] .* 20
        plot!(p2, 1:tf, loc_coral, color=:blue, alpha=0.15, label="")
        plot!(p2, 1:tf, loc_cots, color=:red, alpha=0.15, label="")
    end
    
    plot!(p2, 1:tf, reef_coral .* 100, label="Coral % Mean", color=:blue, linewidth=3)
    plot!(p2, 1:tf, reef_cots .* 20, label="Adult COTS (x20)", color=:red, linewidth=3)
    
    push!(cots_plot_list, p2)
end

cots_final_plot = plot(cots_plot_list..., layout=(rows, cols), size=(1200, rows * 350), margin=5Plots.mm)
savefig(cots_final_plot, joinpath(@__DIR__, "gbr_key_reefs_coral_vs_cots.png"))
println("Saved Coral vs COTS plot to gbr_key_reefs_coral_vs_cots.png")

# --- Plot 2: Regional Overviews ---
println("Generating Regional Overview plot...")
if :management_area in propertynames(dom.loc_data)
    regions = unique(dom.loc_data.management_area)
    reg_plot_list = []
    
    for reg in regions
        # Skip empty regions or "NaN" regions
        if ismissing(reg) || string(reg) == "NaN" || string(reg) == ""
            continue
        end
        
        # Get all locations in this region
        all_locs = findall(dom.loc_data.management_area .== reg)
        
        # Only aggregate over the locations we actually seeded, otherwise the effect is completely diluted
        locs = intersect(all_locs, seed_locs)
        
        if isempty(locs)
            continue
        end
        
        reg_coral_no_cots = vec(mean(total_coral_no_cots[:, locs], dims=2))
        reg_coral_cots = vec(mean(total_coral_cots[:, locs], dims=2))
        
        max_val = max(maximum(reg_coral_no_cots), maximum(reg_coral_cots)) * 100
        p_reg = plot(title="Region: $reg (Seeded Reefs Only)", ylabel="Mean Coral %", legend=:topleft, ylims=(0, max_val*1.2))
        
        plot!(p_reg, 1:tf, reg_coral_no_cots .* 100, label="No COTS", color=:green, linewidth=3)
        plot!(p_reg, 1:tf, reg_coral_cots .* 100, label="With COTS", color=:blue, linewidth=3)
        
        # Add adult COTS density on a secondary axis if possible, or just plot it scaled
        reg_cots = vec(mean(adult_cots[:, locs], dims=2))
        plot!(p_reg, 1:tf, reg_cots .* 20, label="Adult COTS (x20)", color=:red, linewidth=2, linestyle=:dash)
        
        push!(reg_plot_list, p_reg)
    end
    
    if !isempty(reg_plot_list)
        cols_reg = 2
        rows_reg = ceil(Int, length(reg_plot_list) / cols_reg)
        reg_final_plot = plot(reg_plot_list..., layout=(rows_reg, cols_reg), size=(1200, rows_reg * 350), margin=5Plots.mm)
        savefig(reg_final_plot, joinpath(@__DIR__, "gbr_regional_overview.png"))
        println("Saved Regional Overview plot to gbr_regional_overview.png")
    end
else
    println("Warning: :management_area not found. Cannot generate regional overview.")
end

# --- Plot 3: GBR COTS Heatmap ---
println("Generating GBR COTS Heatmap...")
p_cots_heatmap = heatmap(1:tf, 1:size(adult_cots, 2), adult_cots'; 
                         xlabel = "Time Step", 
                         ylabel = "Site ID", 
                         title = "COTS Density Across GBR", 
                         c = :viridis, 
                         yflip = true,
                         size=(1000, 800), margin=5Plots.mm)
savefig(p_cots_heatmap, joinpath(@__DIR__, "gbr_cots_heatmap.png"))
println("Saved COTS Heatmap to gbr_cots_heatmap.png")
