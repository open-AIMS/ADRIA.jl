using Pkg; Pkg.activate(".")
using ADRIA
using CSV, DataFrames, Statistics, Random

println("Loading Lizard Island historical domain...")
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "historical")
site_to_reef = CSV.read("sandbox/data/Lizard_Historical_v0.1/site_to_reef.csv", DataFrame)

clean_reef_names = [split(r, " (")[1] for r in site_to_reef.reef_name]
site_to_reef.reef_name_clean = clean_reef_names
unique_reefs_clean = unique(clean_reef_names)

# Best parameter candidate from Phase 9b (run_id = 2, loss = 23.0897):
best_a_F = 1.34378
best_a_S = 0.0957831
best_IMM = 0.0947791
best_seed_mult = 2.23695
best_pulse_relative_magnitude = parse(Float64, get(ENV, "COTS_PULSE_RELATIVE_MAGNITUDE", "0.0"))

# Stochastic diagnostics configuration.
# Modes:
# - "environmental": only ADRIA environmental/sample fields vary.
# - "demographic": environmental fields plus controlled COTS demographic jitter.
stochastic_mode = lowercase(get(ENV, "COTS_STOCHASTIC_MODE", "demographic"))
stochastic_seed = parse(Int, get(ENV, "COTS_STOCHASTIC_SEED", "20260709"))
cots_demographic_cv = parse(Float64, get(ENV, "COTS_DEMOGRAPHIC_CV", "0.08"))
cots_seed_cv = parse(Float64, get(ENV, "COTS_SEED_MULT_CV", "0.05"))
cots_pulse_relative_cv = parse(Float64, get(ENV, "COTS_PULSE_RELATIVE_MAGNITUDE_CV", "0.0"))
seed_first_n = parse(Int, get(ENV, "COTS_SEED_FIRST_N", "10"))
pulse_start = parse(Int, get(ENV, "COTS_PULSE_START", "1"))
pulse_duration = max(1, parse(Int, get(ENV, "COTS_PULSE_DURATION", "1")))
pulse_repeat_interval = parse(Int, get(ENV, "COTS_PULSE_REPEAT_INTERVAL", "0"))
output_tag = strip(get(ENV, "COTS_OUTPUT_TAG", ""))
output_prefix = output_tag == "" ? "best_stochastic" : "best_stochastic_$(output_tag)"

# Set up environment variables for Cairns Initiation Box pulse & initial density.
# Per-run seed multiplier and pulse amplitude are set immediately before run_scenario.
ENV["COTS_EXTERNAL_PULSE"] = get(ENV, "COTS_EXTERNAL_PULSE", "false")
ENV["COTS_PULSE_START"] = string(pulse_start)
ENV["COTS_PULSE_DURATION"] = string(pulse_duration)
ENV["COTS_PULSE_REPEAT_INTERVAL"] = string(pulse_repeat_interval)
ENV["COTS_PULSE_RELATIVE_MAGNITUDE"] = string(best_pulse_relative_magnitude)
# The pulse code injects into the explicit COTS seed locations. Use the first
# sites as a simple upstream initiation box for this diagnostic workflow.
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = string(seed_first_n)

# Run stochastic simulation runs of our single best parameter candidate.
# Sobol sampling requires a power of 2.
N_scens = parse(Int, get(ENV, "COTS_N_STOCHASTIC_SCENS", "16"))
println("Sampling $N_scens stochastic scenarios...")
scens = ADRIA.sample(dom, N_scens)

p_df = ADRIA.param_table(dom)
for col in names(scens)
    if !(col in ["dhw_scenario", "wave_scenario", "cyclone_mortality_scenario"]) && !startswith(col, "surv_") && !in(col, ["fecundity", "coral_recruitment"])
        scens[!, col] .= p_df[1, col]
    end
end

# Lock calibrated COTS parameters across all scenarios to our single best candidate.
scens.a_F .= best_a_F
scens.a_S .= best_a_S
scens.IMM .= best_IMM

rng = MersenneTwister(stochastic_seed)
lognormal_multiplier(rng::AbstractRNG, cv::Float64) = exp(randn(rng) * cv)
clamp_param(x::Float64, lo::Float64, hi::Float64) = min(max(x, lo), hi)

meta_df = DataFrame(
    sim_id = Int[],
    stochastic_mode = String[],
    stochastic_seed = Int[],
    seed_mult = Float64[],
    external_pulse_enabled = String[],
    pulse_start = Int[],
    pulse_duration = Int[],
    pulse_repeat_interval = Int[],
    pulse_relative_magnitude = Float64[],
    a_ricker = Float64[],
    b_ricker = Float64[],
    m1 = Float64[],
    m2 = Float64[],
    m3 = Float64[],
    p_tilde = Float64[],
    C_max = Float64[],
    tau_condition = Float64[],
    allee_threshold = Float64[],
    imm_threshold = Float64[],
    eta_imm = Float64[]
)

for s in 1:N_scens
    seed_mult = best_seed_mult
    pulse_relative_magnitude = best_pulse_relative_magnitude

    if stochastic_mode == "demographic"
        seed_mult = clamp_param(best_seed_mult * lognormal_multiplier(rng, cots_seed_cv), 0.5, 3.0)
        if get(ENV, "COTS_EXTERNAL_PULSE", "false") == "true" && best_pulse_relative_magnitude > 0.0
            pulse_relative_magnitude = clamp_param(best_pulse_relative_magnitude * lognormal_multiplier(rng, cots_pulse_relative_cv), 0.0, 10.0)
        end
        scens.a_ricker[s] = clamp_param(p_df[1, "a_ricker"] * lognormal_multiplier(rng, cots_demographic_cv), 2.0, 10.0)
        scens.b_ricker[s] = clamp_param(p_df[1, "b_ricker"] * lognormal_multiplier(rng, cots_demographic_cv), 0.01, 0.5)
        scens.m1[s] = clamp_param(p_df[1, "m1"] * lognormal_multiplier(rng, cots_demographic_cv), 0.1, 0.9)
        scens.m2[s] = clamp_param(p_df[1, "m2"] * lognormal_multiplier(rng, cots_demographic_cv), 0.05, 0.5)
        scens.m3[s] = clamp_param(p_df[1, "m3"] * lognormal_multiplier(rng, cots_demographic_cv), 0.05, 0.3)
        scens.p_tilde[s] = clamp_param(p_df[1, "p_tilde"] * lognormal_multiplier(rng, cots_demographic_cv), 0.8, 1.0)
        scens.C_max[s] = clamp_param(p_df[1, "C_max"] * lognormal_multiplier(rng, cots_demographic_cv), 0.4, 1.0)
        scens.tau_condition[s] = clamp_param(p_df[1, "tau_condition"] * lognormal_multiplier(rng, cots_demographic_cv), 1.0, 10.0)
        scens.allee_threshold[s] = clamp_param(p_df[1, "allee_threshold"] * lognormal_multiplier(rng, cots_demographic_cv), 0.1, 5.0)
        scens.imm_threshold[s] = clamp_param(p_df[1, "imm_threshold"] * lognormal_multiplier(rng, cots_demographic_cv), 0.1, 0.8)
        scens.eta_imm[s] = clamp_param(p_df[1, "eta_imm"] * lognormal_multiplier(rng, cots_demographic_cv), 1.0, 5.0)
    elseif stochastic_mode != "environmental"
        error("Unknown COTS_STOCHASTIC_MODE='$stochastic_mode'. Use 'environmental' or 'demographic'.")
    end

    push!(
        meta_df,
        (
            s,
            stochastic_mode,
            stochastic_seed,
            seed_mult,
            get(ENV, "COTS_EXTERNAL_PULSE", "false"),
            pulse_start,
            pulse_duration,
            pulse_repeat_interval,
            pulse_relative_magnitude,
            scens.a_ricker[s],
            scens.b_ricker[s],
            scens.m1[s],
            scens.m2[s],
            scens.m3[s],
            scens.p_tilde[s],
            scens.C_max[s],
            scens.tau_condition[s],
            scens.allee_threshold[s],
            scens.imm_threshold[s],
            scens.eta_imm[s]
        )
    )
end

println("Running $N_scens stochastic simulation runs sequentially in '$stochastic_mode' mode...")
sim_df = DataFrame(
    sim_id = Int[],
    year = Int[],
    reef_name = String[],
    sim_cots_adult = Float64[],
    sim_cots_norm = Float64[],
    sim_coral = Float64[]
)
site_df = DataFrame(
    sim_id = Int[],
    year = Int[],
    reef_name = String[],
    site_index = Int[],
    sim_cots_adult = Float64[],
    sim_coral = Float64[]
)

for s in 1:N_scens
    println("Simulating stochastic run $s / $N_scens ...")
    ENV["COTS_INITIAL_MULTIPLIER"] = string(meta_df.seed_mult[s])
    ENV["COTS_PULSE_RELATIVE_MAGNITUDE"] = string(meta_df.pulse_relative_magnitude[s])

    p = scens[s, :]
    rs = ADRIA.run_scenario(dom, p)
    adult_cots_site = rs.cots_log[:, 3, :]
    total_cover_site = dropdims(sum(rs.raw, dims=(2, 3)), dims=(2, 3))
    n_timesteps = size(adult_cots_site, 1)
    years = 1985:(1984 + n_timesteps)

    for reef in unique_reefs_clean
        site_indices = findall(site_to_reef.reef_name_clean .== reef)
        if isempty(site_indices) continue end

        reef_sim_cots = [mean(adult_cots_site[t, site_indices]) for t in 1:n_timesteps]
        reef_sim_coral = [mean(total_cover_site[t, site_indices]) for t in 1:n_timesteps]

        for t in 1:n_timesteps
            push!(sim_df, (s, years[t], reef, reef_sim_cots[t], 0.0, reef_sim_coral[t]))
            for site_idx in site_indices
                push!(
                    site_df,
                    (
                        s,
                        years[t],
                        reef,
                        site_idx,
                        adult_cots_site[t, site_idx],
                        total_cover_site[t, site_idx]
                    )
                )
            end
        end
    end
end

# Normalize COTS with one shared maximum per reef across all stochastic runs.
for reef in unique(sim_df.reef_name)
    reef_rows = sim_df.reef_name .== reef
    max_sim = maximum(sim_df[reef_rows, :sim_cots_adult])
    if max_sim > 0
        sim_df[reef_rows, :sim_cots_norm] .= sim_df[reef_rows, :sim_cots_adult] ./ max_sim
    end
end

CSV.write("sandbox/data/$(output_prefix)_trajectories.csv", sim_df)
CSV.write("sandbox/data/$(output_prefix)_site_trajectories.csv", site_df)
CSV.write("sandbox/data/$(output_prefix)_metadata.csv", meta_df)
println("Saved stochastic trajectories to sandbox/data/$(output_prefix)_trajectories.csv")
println("Saved site-level stochastic trajectories to sandbox/data/$(output_prefix)_site_trajectories.csv")
println("Saved stochastic metadata to sandbox/data/$(output_prefix)_metadata.csv")