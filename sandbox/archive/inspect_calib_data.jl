using Pkg; Pkg.activate(".")
using CSV, DataFrames
using Parquet2

println("=== COTS Data ===")
df_cots = CSV.read("sandbox/data/reef_cots.csv", DataFrame)
println(first(df_cots, 5))
println("Columns: ", names(df_cots))

println("\n=== Coral Data ===")
df_coral = CSV.read("sandbox/data/reef_manta.csv", DataFrame)
println(first(df_coral, 5))
println("Columns: ", names(df_coral))

println("\n=== DHW Data ===")
ds_dhw = Parquet2.Dataset("sandbox/data/dhw_scens_full.parquet")
# print schema or first few rows
df_dhw = DataFrame(ds_dhw)
println("Columns: ", names(df_dhw))
println("Number of rows: ", nrow(df_dhw))
println("First row sample: ")
println(first(df_dhw, 1))

