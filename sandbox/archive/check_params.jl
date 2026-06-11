using Pkg; Pkg.activate(".")
using ADRIA
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "historical")
p_df = ADRIA.param_table(dom)
println(names(p_df))
