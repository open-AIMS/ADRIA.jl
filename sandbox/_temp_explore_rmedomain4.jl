using ADRIA, DataFrames
dom = ADRIA.load_domain(ADRIA.RMEDomain, "sandbox/data/rme_ml_2025_06_05", "45")

targets = ["18-075", "16-029", "14-116", "21-566"]
found_targets = []
for (i, id) in enumerate(dom.loc_data.GBRMPA_ID)
    if String(id) in targets || startswith(String(id), "14-116")
        push!(found_targets, (i, id))
        println("Found target reef: ", id, " at index ", i)
    end
end
println("Total found: ", length(found_targets))
