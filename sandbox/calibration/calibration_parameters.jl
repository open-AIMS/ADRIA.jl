function calibration_search_space(mode::AbstractString)
    normalized_mode = uppercase(mode)
    focused_names = ["a_F", "a_S", "IMM", "seed_mult"]
    focused_bounds = [(0.1, 2.0), (0.01, 0.9), (0.0, 0.05), (0.5, 4.0)]
    normalized_mode == "FOCUSED" && return focused_names, focused_bounds

    expanded_names = vcat(focused_names, [
        "a_ricker", "b_ricker", "m1", "m2", "m3", "p_tilde", "C_max",
        "tau_condition", "imm_threshold", "eta_imm"
    ])
    expanded_bounds = vcat(focused_bounds, [
        (2.0, 10.0), (0.01, 0.5), (0.1, 0.9), (0.05, 0.5), (0.05, 0.3),
        (0.8, 1.0), (0.4, 1.0), (1.0, 10.0), (0.1, 0.8), (1.0, 5.0)
    ])
    normalized_mode == "EXPANDED" && return expanded_names, expanded_bounds
    if normalized_mode == "EXPANDED_PULSE"
        return vcat(expanded_names, ["pulse_start", "pulse_duration", "pulse_relative_magnitude"]),
            vcat(expanded_bounds, [(15.0, 30.0), (1.0, 4.0), (0.0, 1.0)])
    end
    throw(ArgumentError("Unknown calibration mode: $mode"))
end
