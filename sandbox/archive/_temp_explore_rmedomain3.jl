using ADRIA, DataFrames
dom = ADRIA.load_domain(ADRIA.RMEDomain, "sandbox/data/rme_ml_2025_06_05", "45")

# Let's see if UNIQUE_ID or GBRMPA_ID contains names
println("First 10 UNIQUE_ID:")
println(dom.loc_data.UNIQUE_ID[1:10])

println("First 10 GBRMPA_ID:")
println(dom.loc_data.GBRMPA_ID[1:10])

println("First 10 cscape_region:")
println(dom.loc_data.cscape_region[1:10])

println("Searching for target reefs in UNIQUE_ID...")
targets = ["lizard", "brewer", "batt", "gannet"]
for id in dom.loc_data.UNIQUE_ID
    id_str = lowercase(String(id))
    if any(occursin.(targets, id_str))
        println("Found in UNIQUE_ID: ", id)
    end
end
