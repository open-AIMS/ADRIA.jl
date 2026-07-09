$ErrorActionPreference = "Stop"

$durations = @(1, 2, 3)

foreach ($duration in $durations) {
    $env:COTS_OUTPUT_TAG = "pulse_start20_dur$duration"
    $env:COTS_EXTERNAL_PULSE = "true"
    $env:COTS_PULSE_START = "20"
    $env:COTS_PULSE_DURATION = "$duration"
    $env:COTS_PULSE_REPEAT_INTERVAL = "0"
    $env:COTS_PULSE_RELATIVE_MAGNITUDE = "0.5"

    Write-Host "Running pulse scenario: start=20 duration=$duration relative_magnitude=0.5"
    julia sandbox\calibration\simulate_best_stochastic.jl

    Write-Host "Plotting pulse scenario: $env:COTS_OUTPUT_TAG"
    python sandbox\plotting\plot_best_stochastic.py
}

