using Pkg; Pkg.activate(".")
using ADRIA
using YAXArrays

# 1. Load the domain
println("Loading domain...")
dom = ADRIA.load_domain(LizardDomain, "sandbox/data/Lizard_Island_Cluster_v0.1", "45")

# 2. Setup Parameters for Baseline (No COTS)
println("Setting up baseline scenario...")
p_baseline = ADRIA.param_table(dom)
# Disable COTS completely by setting immigration and initial density to 0 if we can, 
# or just setting ADRIA_COTS_ENABLED="false" in ENV which we can't do per scenario natively yet without ENV.
# But wait, we can just set COTS parameters in p_baseline to zero out their effect:
p_baseline.a_F .= 0.0
p_baseline.a_S .= 0.0
p_baseline.IMM .= 0.0

# 3. Setup Parameters for COTS Enabled
println("Setting up COTS scenario...")
p_cots = ADRIA.param_table(dom)
# Leave COTS enabled with default params

# 4. Run Scenarios
println("Running Baseline...")
ENV["ADRIA_COTS_ENABLED"] = "true"  # Ensure it runs the module
rs_base = ADRIA.run_scenario(dom, p_baseline[1, :])

println("Running COTS...")
rs_cots = ADRIA.run_scenario(dom, p_cots[1, :])

# 5. Extract results
tc_base = rs_base.raw
tc_cots = rs_cots.raw

println("Baseline cover array size: ", size(tc_base))
println("COTS cover array size: ", size(tc_cots))

# Just sum everything to get a single number as a sanity check
println("Baseline total raw cover sum: ", sum(tc_base))
println("COTS total raw cover sum: ", sum(tc_cots))

# If COTS is populated, let's look at the COTS results
cots_pop = rs_cots.cots_log
println("Max COTS population across time: ", maximum(cots_pop))
println("Total COTS population across time: ", sum(cots_pop))

println("Done!")
