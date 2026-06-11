using ADRIA, DataFrames
dom = ADRIA.load_domain(ADRIA.RMEDomain, "sandbox/data/rme_ml_2025_06_05", "45")
println(propertynames(dom.loc_data))
