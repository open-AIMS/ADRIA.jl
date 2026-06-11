using Pkg; Pkg.activate(".")
using ADRIA
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "historical")
p_df = ADRIA.param_table(dom)
rs = ADRIA.run_scenario(dom, p_df[1, :])
println("rs.raw size: ", size(rs.raw))
println("rs.cots_log size: ", size(rs.cots_log))
