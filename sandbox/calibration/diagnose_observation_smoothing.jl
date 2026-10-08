using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using CSV, DataFrames
include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))

const OBSERVATIONS = CSV.read(joinpath(@__DIR__, "..", "data", "reef_cots.csv"), DataFrame)
const TARGETS = ["Lizard Isles", "Macgillivray Reef", "North Direction Island", "Eyrie Reef"]
const MODEL_YEARS = 1985:2024

for reef_name in TARGETS
    rows = sort(OBSERVATIONS[(OBSERVATIONS.reef_name .== reef_name) .&
                             in.(OBSERVATIONS.year, Ref(MODEL_YEARS)), :], :year)
    years = Int.(rows.year)
    cpue = Float64.(rows.cotsptow)
    println("\n", reef_name, ": ", length(years), " survey years; gaps=",
            sort(unique(diff(years))))
    for window_years in (1, 3, 5, 7)
        smooth = calendar_moving_average(years, cpue; window_years=window_years)
        peaks = detect_cots_peaks(years, smooth)
        println("  window=", window_years, " years: peaks=",
                [(p.year, round(p.value; digits=3)) for p in peaks],
                "; max=", round(maximum(smooth); digits=3))
    end
end
