using CpuId, Pkg

if occursin("Intel", cpubrand())
    # Add MKL library to leverage faster performance for matrix operations
    Pkg.add("MKL")
end
