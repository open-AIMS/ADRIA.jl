using ADRIA, DataFrames
dom = ADRIA.load_domain(ADRIA.RMEDomain, "sandbox/data/rme_ml_2025_06_05", "45")
println(typeof(dom.conn))
println(size(dom.conn))
