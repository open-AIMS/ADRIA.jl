using Pkg
Pkg.activate(".")
using ADRIA
try
    dom = ADRIA.load_domain(LizardDomain, "sandbox/data/Lizard_Island_Cluster_v0.1", "45")
    println("Success! Loaded ", length(dom.loc_ids), " sites.")
catch e
    showerror(stdout, e, catch_backtrace())
end
