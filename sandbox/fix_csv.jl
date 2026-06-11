using Pkg; Pkg.activate(".")
using CSV, DataFrames
fn = "sandbox/data/Lizard_Island_Cluster_v0.1/connectivity/Lizard_Connectivity.csv"
df = CSV.read(fn, DataFrame)
if names(df)[1] != "Source"
    insertcols!(df, 1, :Source => names(df))
    CSV.write(fn, df)
end
