using Pkg; Pkg.activate(".")
using ADRIA
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "45")
p = ADRIA.param_table(dom)
open("sandbox/param_names.txt", "w") do io
    println(io, join(names(p), "\n"))
end
