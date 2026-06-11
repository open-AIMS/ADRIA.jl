# toy_cots_cycles.jl — Minimal toy model to diagnose COTS predator-prey cycles
#
# Strips away all ADRIA/CoralBlox complexity and tests whether different
# formulations of the COTS population dynamics can produce 15-year
# boom-bust cycles characteristic of real GBR COTS outbreaks.
#
# Key biological assumptions (from user feedback):
# 1. COTS produce MILLIONS of larvae → the larval explosion is the
#    primary overshoot mechanism
# 2. Adult COTS mortality is near-zero while prey is available, then
#    DRASTIC once preferred prey runs out (threshold, not smooth)
# 3. Maternal body condition tracks nutrition → determines larval output
# 4. The 2-year age structure (age-0 → age-1 → age-2+ adults) provides
#    a natural delay in the numerical response
#
# We test 5 model variants to see which produces 15-year cycles:
#   A: Current model (Type I, instant smooth starvation)
#   B: Type II functional response + instant starvation
#   C: Type I + threshold starvation (step-function mortality)
#   D: Type II + threshold starvation + body condition → larval gating
#   E: Variant D with Ricker recruitment (overcompensation)
#
# Usage: julia --project=sandbox sandbox/toy_cots_cycles.jl

using Pkg
Pkg.activate(dirname(@__FILE__))

using Plots
using Statistics
using StaticArrays
using Printf

# ================================================================
# Coral Growth (simple logistic, shared by all variants)
# ================================================================
struct CoralState
    F::Float64   # Fast-growing coral cover (0-1 relative)
    S::Float64   # Slow-growing coral cover (0-1 relative)
end

function coral_step(c::CoralState, Cons_F::Float64, Cons_S::Float64;
                    r_F=0.12, r_S=0.06, K_F=0.45, K_S=0.30)
    total = c.F + c.S
    K_total = K_F + K_S

    # Logistic growth minus consumption
    F_new = c.F + r_F * c.F * (1.0 - total / K_total) - Cons_F
    S_new = c.S + r_S * c.S * (1.0 - total / K_total) - Cons_S

    return CoralState(clamp(F_new, 0.0, 1.0), clamp(S_new, 0.0, 1.0))
end

# ================================================================
# COTS State (shared structure for all variants)
# ================================================================
mutable struct COTSState
    N::MVector{3, Float64}        # Age 0 (larvae/recruits), Age 1 (juveniles), Age 2+ (adults)
    body_condition::Float64       # Maternal nutrition index (0=starving, 1=well-fed)
end

# ================================================================
# Variant A: Current model (Type I, instant smooth starvation)
# ================================================================
function cots_step_A!(cots::COTSState, F::Float64, S::Float64, p::NamedTuple)
    N = cots.N
    total_coral = max(0.0, F + S)

    # Smooth starvation (current model) — survival drops smoothly with coral
    rho = min(1.0, (total_coral / p.C_max)^2)
    f = (1.0 - p.p_tilde) + p.p_tilde * rho  # 0.15 when coral=0, 1.0 when coral>=C_max

    # Beverton-Holt recruitment
    R = (p.a * N[3]) / (1.0 + p.b * N[3])

    # Survival rates
    s1 = 1.0 - p.m1
    s2 = 1.0 - p.m2
    s3 = 1.0 - p.m3

    # State update (age-structured)
    N_next = MVector{3, Float64}(0.0, 0.0, 0.0)
    N_next[1] = R + p.IMM
    N_next[2] = N[1] * s1
    N_next[3] = N[2] * (s2 * f) + N[3] * (s3 * f)

    # Type I consumption (linear, no saturation)
    Cons_F = N[3] * p.a_F * F
    Cons_S = N[3] * p.a_S * S
    Cons_F = min(Cons_F, F)
    Cons_S = min(Cons_S, S)

    cots.N .= max.(0.0, N_next)
    return Cons_F, Cons_S
end

# ================================================================
# Variant B: Type II functional response + instant starvation
# ================================================================
function cots_step_B!(cots::COTSState, F::Float64, S::Float64, p::NamedTuple)
    N = cots.N
    total_coral = max(0.0, F + S)

    # Same smooth starvation as A
    rho = min(1.0, (total_coral / p.C_max)^2)
    f = (1.0 - p.p_tilde) + p.p_tilde * rho

    R = (p.a * N[3]) / (1.0 + p.b * N[3])

    s1 = 1.0 - p.m1
    s2 = 1.0 - p.m2
    s3 = 1.0 - p.m3

    N_next = MVector{3, Float64}(0.0, 0.0, 0.0)
    N_next[1] = R + p.IMM
    N_next[2] = N[1] * s1
    N_next[3] = N[2] * (s2 * f) + N[3] * (s3 * f)

    # Type II: consumption SATURATES at high prey density
    # This means at high coral, each COTS eats at a fixed max rate
    # → coral can grow faster than COTS can eat it → overshoot
    h = p.h  # handling time
    denom = 1.0 + h * (p.a_F * F + p.a_S * S)
    Cons_F = N[3] * (p.a_F * F) / denom
    Cons_S = N[3] * (p.a_S * S) / denom
    Cons_F = min(Cons_F, F)
    Cons_S = min(Cons_S, S)

    cots.N .= max.(0.0, N_next)
    return Cons_F, Cons_S
end

# ================================================================
# Variant C: Type I + threshold starvation (step-function mortality)
# ================================================================
function cots_step_C!(cots::COTSState, F::Float64, S::Float64, p::NamedTuple)
    N = cots.N
    total_coral = max(0.0, F + S)

    # THRESHOLD starvation: mortality stays LOW while prey is available,
    # then jumps DRASTICALLY once prey drops below a critical threshold.
    # This is biologically accurate: COTS are fine until they literally
    # have nothing to eat, then they starve quickly.
    starve_threshold = p.C_max * 0.15  # below 15% of C_max → starvation kicks in
    if total_coral > starve_threshold
        f = 1.0   # no extra mortality while food is available
    else
        # Sharp crash: survival drops as a steep function near zero coral
        frac = total_coral / starve_threshold
        f = (1.0 - p.p_tilde) + p.p_tilde * frac^3  # cubic for sharpness
    end

    R = (p.a * N[3]) / (1.0 + p.b * N[3])

    s1 = 1.0 - p.m1
    s2 = 1.0 - p.m2
    s3 = 1.0 - p.m3

    N_next = MVector{3, Float64}(0.0, 0.0, 0.0)
    N_next[1] = R + p.IMM
    N_next[2] = N[1] * s1
    N_next[3] = N[2] * (s2 * f) + N[3] * (s3 * f)

    # Type I consumption
    Cons_F = N[3] * p.a_F * F
    Cons_S = N[3] * p.a_S * S
    Cons_F = min(Cons_F, F)
    Cons_S = min(Cons_S, S)

    cots.N .= max.(0.0, N_next)
    return Cons_F, Cons_S
end

# ================================================================
# Variant D: Type II + threshold starvation + body condition
#            Body condition gates LARVAL PRODUCTION (maternal nutrition)
# ================================================================
function cots_step_D!(cots::COTSState, F::Float64, S::Float64, p::NamedTuple)
    N = cots.N
    total_coral = max(0.0, F + S)

    # --- Body condition: exponential moving average of food availability ---
    # Represents maternal nutrition — determines how many viable larvae
    # are produced. Responds SLOWLY to changes in food.
    # tau = memory timescale in years (higher = more lag)
    tau = get(p, :tau_condition, 3.0)  # 3-year memory
    alpha = 1.0 / tau  # smoothing factor
    food_signal = min(1.0, total_coral / (p.C_max * 0.5))  # normalised food
    new_condition = (1.0 - alpha) * cots.body_condition + alpha * food_signal
    cots.body_condition = clamp(new_condition, 0.0, 1.0)

    # --- THRESHOLD starvation: mortality stays low while prey available ---
    starve_threshold = p.C_max * 0.15
    if total_coral > starve_threshold
        f = 1.0
    else
        frac = total_coral / starve_threshold
        f = (1.0 - p.p_tilde) + p.p_tilde * frac^3
    end

    # --- Recruitment gated by body condition (maternal nutrition) ---
    # COTS produce millions of larvae when well-fed.
    # Body condition determines what fraction of potential fecundity is realized.
    # This creates a LAGGED positive feedback: good conditions now →
    # larvae explosion 2-3 years later (when they mature to adults).
    fecundity = p.a * cots.body_condition^2  # squared: emphasise boom when well-fed
    R = (fecundity * N[3]) / (1.0 + p.b * N[3])

    s1 = 1.0 - p.m1
    s2 = 1.0 - p.m2
    s3 = 1.0 - p.m3

    N_next = MVector{3, Float64}(0.0, 0.0, 0.0)
    N_next[1] = R + p.IMM
    N_next[2] = N[1] * s1
    N_next[3] = N[2] * (s2 * f) + N[3] * (s3 * f)

    # Type II functional response (handling time saturation)
    h = p.h
    denom = 1.0 + h * (p.a_F * F + p.a_S * S)
    Cons_F = N[3] * (p.a_F * F) / denom
    Cons_S = N[3] * (p.a_S * S) / denom
    Cons_F = min(Cons_F, F)
    Cons_S = min(Cons_S, S)

    cots.N .= max.(0.0, N_next)
    return Cons_F, Cons_S
end

# ================================================================
# Variant E: Same as D but with Ricker recruitment
#            Ricker allows OVERCOMPENSATION → dramatic spikes
# ================================================================
function cots_step_E!(cots::COTSState, F::Float64, S::Float64, p::NamedTuple)
    N = cots.N
    total_coral = max(0.0, F + S)

    # Body condition (same as D)
    tau = get(p, :tau_condition, 3.0)
    alpha = 1.0 / tau
    food_signal = min(1.0, total_coral / (p.C_max * 0.5))
    new_condition = (1.0 - alpha) * cots.body_condition + alpha * food_signal
    cots.body_condition = clamp(new_condition, 0.0, 1.0)

    # Threshold starvation (same as D)
    starve_threshold = p.C_max * 0.15
    if total_coral > starve_threshold
        f = 1.0
    else
        frac = total_coral / starve_threshold
        f = (1.0 - p.p_tilde) + p.p_tilde * frac^3
    end

    # --- RICKER recruitment: allows overcompensation ---
    # R = a * N * exp(-b * N) * body_condition
    # At low N: recruitment scales with N (more adults → more larvae)
    # At high N: density dependence REDUCES recruitment (overcrowding)
    # The peak at N = 1/b creates a massive pulse of recruits
    # followed by a crash — exactly the boom-bust pattern we want.
    fecundity = p.a_ricker * cots.body_condition^2
    R = fecundity * N[3] * exp(-p.b_ricker * N[3])

    s1 = 1.0 - p.m1
    s2 = 1.0 - p.m2
    s3 = 1.0 - p.m3

    N_next = MVector{3, Float64}(0.0, 0.0, 0.0)
    N_next[1] = R + p.IMM
    N_next[2] = N[1] * s1
    N_next[3] = N[2] * (s2 * f) + N[3] * (s3 * f)

    # Type II functional response
    h = p.h
    denom = 1.0 + h * (p.a_F * F + p.a_S * S)
    Cons_F = N[3] * (p.a_F * F) / denom
    Cons_S = N[3] * (p.a_S * S) / denom
    Cons_F = min(Cons_F, F)
    Cons_S = min(Cons_S, S)

    cots.N .= max.(0.0, N_next)
    return Cons_F, Cons_S
end

# ================================================================
# Run a simulation for any variant
# ================================================================
function run_toy_simulation(step_fn!::Function, params::NamedTuple;
                            F_init=0.35, S_init=0.25,
                            N_init=(0.01, 0.01, 0.05),
                            body_init=0.8,
                            steps=120,
                            n_reefs=1)
    # Support n_reefs for future multi-reef extension
    # For now just single-reef
    results = []

    for reef in 1:n_reefs
        coral = CoralState(F_init, S_init)
        cots = COTSState(MVector{3, Float64}(N_init...), body_init)

        F_hist = zeros(steps)
        S_hist = zeros(steps)
        cots_hist = zeros(steps)    # adult density
        larvae_hist = zeros(steps)  # age-0 recruits
        bc_hist = zeros(steps)      # body condition

        for t in 1:steps
            F_hist[t] = coral.F
            S_hist[t] = coral.S
            cots_hist[t] = cots.N[3]
            larvae_hist[t] = cots.N[1]
            bc_hist[t] = cots.body_condition

            Cons_F, Cons_S = step_fn!(cots, coral.F, coral.S, params)
            coral = coral_step(coral, Cons_F, Cons_S)
        end

        push!(results, (F=F_hist, S=S_hist, cots=cots_hist,
                        larvae=larvae_hist, body_condition=bc_hist))
    end

    return n_reefs == 1 ? results[1] : results
end

# ================================================================
# Cycle detection
# ================================================================
function count_cycles(series; min_peak_ratio=3.0, min_separation=4)
    n = length(series)
    peaks = Int[]
    troughs = Int[]

    for t in 3:(n-2)
        if series[t] > series[t-1] && series[t] > series[t-2] &&
           series[t] >= series[t+1] && series[t] >= series[t+2]
            if isempty(peaks) || (t - peaks[end]) > min_separation
                push!(peaks, t)
            elseif series[t] > series[peaks[end]]
                peaks[end] = t
            end
        end
    end

    if length(peaks) < 2
        return (n_cycles=0, period=0.0, peaks=peaks, amp_ratio=0.0)
    end

    # Find troughs between consecutive peaks
    for i in 1:(length(peaks)-1)
        segment = series[peaks[i]:peaks[i+1]]
        _, idx = findmin(segment)
        push!(troughs, peaks[i] + idx - 1)
    end

    # Check amplitude ratios
    real_cycles = 0
    for i in 1:length(troughs)
        peak_val = max(series[peaks[i]], series[peaks[i+1]])
        trough_val = series[troughs[i]]
        ratio = trough_val > 0 ? peak_val / trough_val : (peak_val > 0 ? Inf : 0.0)
        if ratio >= min_peak_ratio
            real_cycles += 1
        end
    end

    periods = diff(peaks)
    avg_period = mean(periods)
    peak_vals = series[peaks]
    trough_vals = length(troughs) > 0 ? series[troughs] : [0.0]
    amp_ratio = mean(trough_vals) > 0 ? mean(peak_vals) / mean(trough_vals) : Inf

    return (n_cycles=real_cycles, period=avg_period, peaks=peaks,
            troughs=troughs, amp_ratio=amp_ratio)
end

# ================================================================
# Parameter sets
# ================================================================
# Shared base parameters
base_params = (
    a = 2.5,          # Beverton-Holt: max recruitment rate
    b = 0.3,          # Beverton-Holt: density dependence
    IMM = 0.001,      # Background immigration (low — allows re-seeding)
    p_tilde = 0.92,   # Max starvation mortality fraction (92% of base can be lost)
    C_max = 0.75,     # Total coral cover at which starvation = 0
    m1 = 0.4,         # Age-0 mortality
    m2 = 0.2,         # Age-1 mortality
    m3 = 0.08,        # Adult base mortality (LOW — adults survive well when fed)
    a_F = 0.5,        # Fast coral consumption rate
    a_S = 0.12,       # Slow coral consumption rate
    h = 2.5,          # Handling time for Type II (variants B,D,E)
    tau_condition = 3.0,  # Body condition memory (years) for variants D,E
    a_ricker = 4.0,       # Ricker recruitment rate for variant E
    b_ricker = 0.8,       # Ricker density dependence for variant E
)

STEPS = 120  # 120 years for multiple cycles

# ================================================================
# Run all variants
# ================================================================
println("=" ^ 70)
println("TOY MODEL: COTS Predator-Prey Cycle Diagnostics")
println("=" ^ 70)

variant_names = ["A: Type I + Smooth Starvation (current)",
                 "B: Type II + Smooth Starvation",
                 "C: Type I + Threshold Starvation",
                 "D: Type II + Threshold + Body Condition",
                 "E: Type II + Threshold + Body Cond + Ricker"]

variant_fns = [cots_step_A!, cots_step_B!, cots_step_C!, cots_step_D!, cots_step_E!]

results = []
for (i, (name, fn)) in enumerate(zip(variant_names, variant_fns))
    println("\nRunning Variant $name ...")
    r = run_toy_simulation(fn, base_params; steps=STEPS)
    coral_total = r.F .+ r.S
    cycles = count_cycles(r.cots)

    push!(results, (name=name, result=r, coral=coral_total, cycles=cycles))

    println("  COTS: peak=$(round(maximum(r.cots), digits=3)), final=$(round(r.cots[end], digits=4))")
    println("  Coral: min=$(round(minimum(coral_total)*100, digits=1))%, final=$(round(coral_total[end]*100, digits=1))%")
    println("  Cycles: $(cycles.n_cycles), period=$(round(cycles.period, digits=1))yr, amp_ratio=$(round(cycles.amp_ratio, digits=1))")
    if length(cycles.peaks) > 0
        println("  Peak times: $(cycles.peaks)")
    end
end

# ================================================================
# Print detailed time series for the best variant
# ================================================================
# Find the variant with the most cycles and closest to 15yr period
best_idx = argmax([r.cycles.n_cycles > 0 ?
    r.cycles.n_cycles * exp(-0.5*((r.cycles.period - 15.0)/5.0)^2) : 0.0
    for r in results])

println("\n" * "=" ^ 70)
println("BEST VARIANT: $(results[best_idx].name)")
println("=" ^ 70)
r = results[best_idx]
println("\n  t | COTS Adults | Larvae (Age0) | Coral %  | Body Cond")
println("  " * "-"^65)
for t in 1:STEPS
    @printf("  %3d | %10.4f | %13.4f | %6.1f%% | %.3f\n",
        t, r.result.cots[t], r.result.larvae[t],
        r.coral[t]*100, r.result.body_condition[t])
end

# ================================================================
# Plotting
# ================================================================
Plots.gr()

# --- Main comparison: 5 panels ---
variant_plots = []
for (i, r) in enumerate(results)
    c = r.cycles
    title_str = "$(r.name)\nCycles=$(c.n_cycles) Period=$(round(c.period, digits=1))yr"

    pl = plot(1:STEPS, r.coral .* 100,
        label="Coral Cover (%)", color=:blue, lw=2,
        ylabel="Coral %", ylims=(0, 100),
        title=title_str, titlefontsize=7,
        legend=:topleft, legendfontsize=6)

    pl2 = twinx(pl)
    plot!(pl2, 1:STEPS, r.result.cots,
        label="COTS Adults", color=:red, lw=2,
        ylabel="COTS Density", legend=:topright, legendfontsize=6)

    # Mark peaks
    if length(c.peaks) > 0
        scatter!(pl2, c.peaks, r.result.cots[c.peaks],
            color=:black, markersize=4, label="", markershape=:dtriangle)
    end

    push!(variant_plots, pl)
end

p_comparison = plot(variant_plots...,
    layout=(length(variant_plots), 1),
    size=(900, 300 * length(variant_plots)),
    margin=6Plots.mm)

savefig(p_comparison, joinpath(@__DIR__, "toy_cots_variants.png"))
println("\nSaved variant comparison: sandbox/toy_cots_variants.png")

# --- Phase portrait for best variant ---
r_best = results[best_idx]
p_phase = plot(r_best.coral .* 100, r_best.result.cots,
    xlabel="Total Coral Cover (%)", ylabel="COTS Adult Density",
    title="Phase Portrait: $(r_best.name)",
    lw=1.5, color=:purple, legend=false, arrow=true)
scatter!([r_best.coral[1]*100], [r_best.result.cots[1]],
    color=:green, markersize=8, label="Start", markershape=:circle)
scatter!([r_best.coral[end]*100], [r_best.result.cots[end]],
    color=:red, markersize=8, label="End", markershape=:square)

savefig(p_phase, joinpath(@__DIR__, "toy_cots_phase.png"))
println("Saved phase portrait: sandbox/toy_cots_phase.png")

# --- Body condition plot for variants D and E ---
p_bc = plot(ylabel="Body Condition", xlabel="Year",
    title="Body Condition (Maternal Nutrition) Over Time",
    legend=:topright)
for i in [4, 5]  # D and E
    r = results[i]
    plot!(p_bc, 1:STEPS, r.result.body_condition,
        label=r.name[1:1], lw=2)
end
savefig(p_bc, joinpath(@__DIR__, "toy_cots_body_condition.png"))
println("Saved body condition: sandbox/toy_cots_body_condition.png")

# --- Summary table ---
println("\n" * "=" ^ 70)
println("VARIANT COMPARISON SUMMARY")
println("=" ^ 70)
println(@sprintf("%-45s | %6s | %6s | %8s | %s",
    "Variant", "Cycles", "Period", "AmpRatio", "Verdict"))
println("-"^90)
for r in results
    c = r.cycles
    verdict = if c.n_cycles >= 3 && 10 <= c.period <= 20
        "✓ TARGET MET"
    elseif c.n_cycles >= 2
        "~ Partial"
    elseif c.n_cycles == 1
        "○ Single bump"
    else
        "✗ No cycles"
    end
    println(@sprintf("%-45s | %6d | %5.1fyr | %8.1f | %s",
        r.name, c.n_cycles, c.period, c.amp_ratio, verdict))
end

println("\nDone. Check sandbox/ for plots.")
