# Measures the *cold* JIT cost of the domain-model assembly + factor write-back path,
# i.e. the same work as the heavy block of `@compile_workload` in `src/ADRIA.jl`.
#
# Must be run in a process where ADRIA was precompiled with `precompile_workload=false`
# (see `precompile_budget.jl`); otherwise the workload has already cached this work and
# the measured time is meaningless.
#
# Prints a single line `COLD_MODEL_ASSEMBLY_SECONDS=<float>` and exits non-zero if the
# budget (ARGS[1], seconds) is exceeded.

using ADRIA
const A = ADRIA

budget = parse(Float64, get(ARGS, 1, "120"))

t = @elapsed begin
    arr = ones(Float32, 1, 1, 1)
    model = A._assemble_domain_model(
        A.EnvironmentalLayer(arr, arr, arr),
        A.Intervention(),
        A.CAqCriteriaWeights(),
        A.FogCriteriaWeights(),
        A.LvMCriteriaWeights(),
        A.DepthThresholds(),
        A.Coral(),
        A.GrowthAcceleration(),
        A.DepthAttenuation(),
        A.LowDHWBleaching(),
        A.CotsParams()
    )
    for _ in 1:2
        ms = A.model_spec(model)
        model[:val] = Float64.(ms.val)
        model[:dist_params] = Vector{Float64}[collect(Float64, d) for d in ms.dist_params]
    end
end

println("COLD_MODEL_ASSEMBLY_SECONDS=", round(t; digits=1))
exit(t <= budget ? 0 : 1)
