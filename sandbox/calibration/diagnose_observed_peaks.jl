using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPO_ROOT, "sandbox"))

using CSV, DataFrames, Statistics
include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))

observations = CSV.read(joinpath(REPO_ROOT, "sandbox", "data", "reef_cots.csv"), DataFrame)
reef_names = ["Lizard Isles", "Macgillivray Reef", "North Direction Island", "Eyrie Reef"]

for reef_name in reef_names
    reef = observations[observations.reef_name .== reef_name, :]
    yearly = combine(groupby(reef, :year), :cotsptow => mean => :value)
    println(reef_name)
    for threshold in (0.25, 0.4, 0.5)
        peaks = detect_cots_peaks(Int.(yearly.year), Float64.(yearly.value); min_height=threshold)
        description = join(["$(peak.year):$(round(peak.value; digits=3))" for peak in peaks], ", ")
        println("  relative threshold $threshold -> $description")
    end
end
