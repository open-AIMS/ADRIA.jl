using ADRIA
using DataFrames

println("Loading RMEDomain...")
dom = ADRIA.load_domain(ADRIA.RMEDomain, "sandbox/data/rme_ml_2025_06_05", "45")
println("Loaded domain successfully!")

reef_col = :Reef
if !(:Reef in propertynames(dom.loc_data))
    reef_col = :reef_siteid
end

println("Unique reefs: ", length(unique(dom.loc_data[!, reef_col])))

targets = ["lizard", "brewer", "batt", "gannet"]
found_targets = []
for r in unique(dom.loc_data[!, reef_col])
    r_str = lowercase(String(r))
    if any(occursin.(targets, r_str))
        push!(found_targets, r)
        println("Found target reef: ", r)
    end
end
println("Total found: ", length(found_targets))
