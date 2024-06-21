using CUDA
using BenchmarkTools

#=
Compare simple operations on CPU vs GPU to get an ideal of the performance
gains settler_cover may have.
=#

n1 = 6
n2 = 3806


@info "GPU"
begin
    x = CUDA.rand(n2, n2)
    y = CUDA.rand(n1, n2)
    z = CUDA.zeros(n2, n1)

    gpu_bm = @benchmark begin
        x .= CUDA.rand(n2, n2)
        y .= CUDA.rand(n1, n2)
        z .= x * y'
    end

    display(gpu_bm)
end


@info "CPU"
begin
    x = rand(n2, n2)
    y = rand(n1, n2)
    z = zeros(n2, n1)

    cpu_bm = @benchmark begin
        x .= rand(n2, n2)
        y .= rand(n1, n2)
        z .= x * y'
    end
    display(cpu_bm)
end

@info "ratio CPU vs GPU"
ratio(median(cpu_bm), median(gpu_bm))
