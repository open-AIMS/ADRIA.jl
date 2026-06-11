using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using ADRIA

# We need to run a small batch of scenarios to see if COTS logs are saved to Zarr.
DOMAIN_PATH = joinpath(@__DIR__, "data", "rme_ml_2025_06_05")
dom = ADRIA.load_domain(ADRIA.RMEDomain, DOMAIN_PATH, "45")

# Enable COTS
ENV["ADRIA_COTS_ENABLED"] = "true"

scens = ADRIA.sample(dom, 2)
# Disable MC and Seeding to avoid unrelated `mc_target_locations` bug in this domain
scens.guided .= 0
scens.seed_TA .= 0
scens.seed_CA .= 0
scens.fogging .= 0
scens.mc_freq .= 0
try
    rs = ADRIA.run_scenarios(dom, scens, "45")

    println("Checking if COTS logs exist in ResultSet:")
    println("cots_pop_log: ", !isnothing(rs.cots_pop_log) ? size(rs.cots_pop_log) : "nothing")
    println("cots_condition_log: ", !isnothing(rs.cots_condition_log) ? size(rs.cots_condition_log) : "nothing")
catch e
    open(joinpath(@__DIR__, "error.log"), "w") do f
        for (exc, bt) in Base.catch_stack()
            showerror(f, exc, bt)
            println(f)
        end
    end
end
