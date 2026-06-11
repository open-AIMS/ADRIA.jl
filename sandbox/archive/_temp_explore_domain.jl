using ADRIA
using DataFrames

dom = ADRIA.load_domain("sandbox/data/rme_ml_2025_06_05", "45")
println(propertynames(dom.loc_data))
println(unique(dom.loc_data.Reef))

# Let's find Lizard Island, John Brewer Reef, Batt Reef, Gannet Cay
for r in unique(dom.loc_data.Reef)
    r_str = lowercase(String(r))
    if occursin("lizard", r_str) || occursin("brewer", r_str) || occursin("batt", r_str) || occursin("gannet", r_str)
        println("Found target reef: ", r)
    end
end
