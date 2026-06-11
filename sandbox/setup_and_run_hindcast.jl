using Pkg
Pkg.activate(".")
# Install missing dependencies
Pkg.instantiate()
println("Setting up sandbox environment...")
cd("sandbox")
Pkg.activate(".")
# Develop ADRIA locally using the sandbox architecture (from inside sandbox folder)
Pkg.develop(path="../")

# Force the registry version of CoralBlox
# Pkg.add("CoralBlox")



println("Environment is ready. Loading ADRIA...")
using ADRIA

# --- Run a Hindcast Scenario (2008 - 2024) ---
# We load the ReefMod domain explicitly for the hindcast period
# (Note: make sure the path to the ReefMod domain matches your local system)
# Check the default RME data package location:
RME_PATH = "C:/Users/smatthew/DataPackages/rme_ml_2025_06_05"
# If not there, please update RME_PATH to the correct path on your machine.

if isdir(RME_PATH)
    println("Loading domain from $RME_PATH for 2008-2024...")
    dom = ADRIA.load_domain(RMEDomain, RME_PATH, "45", timeframe=(2008, 2024))

    println("Running baseline hindcast scenario...")
    scens = ADRIA.sample(dom, 16) # Generate 1 scenario
    rs = ADRIA.run_scenario(dom, scens[1, :])

    println("Hindcast run complete! Results shape: ", size(rs.raw))
else
    println("WARNING: ReefMod domain not found at $RME_PATH. Please run this script again after updating the path.")
end
