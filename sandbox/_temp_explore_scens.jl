using ADRIA, DataFrames
dom = ADRIA.load_domain(ADRIA.RMEDomain, "sandbox/data/rme_ml_2025_06_05", "45")
scens = ADRIA.sample(dom, 2)
println("Columns of scens: ")
println(propertynames(scens))
if :guided in propertynames(scens)
    println("guided type: ", typeof(scens.guided))
    println("guided values: ", scens.guided)
end
