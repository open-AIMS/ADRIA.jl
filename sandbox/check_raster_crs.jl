using Pkg; Pkg.activate(".")
using ArchGDAL

ArchGDAL.read("sandbox/data/COTS_prob_0.02_cpue_year2025_clean.tif") do ds
    println("Raster CRS: ", ArchGDAL.getproj(ds))
end
