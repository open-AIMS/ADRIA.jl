# Helpers declaring each factor as a `ParamSpec` (constant at its default unless made
# otherwise, with its full range in `dist_params`); `spec_rows` maps them to table rows.
_cots_real(name, val, (lo, hi), label, description) = ParamSpec(
    name, val, (lo, hi); kind=:real, dist=:triang, dist_params=[lo, hi, val],
    label=label, description=description
)
_cots_switch(name, val, label, description) = ParamSpec(
    name, val, (0, 1); kind=:bool, label=label, description=description
)
_cots_levels(name, val, levels, label, description) = ParamSpec(
    name, val, (minimum(levels), maximum(levels)); kind=:categorical, levels=levels,
    label=label, description=description
)
_cots_ordered(name, val, (lo, hi, step), label, description) = ParamSpec(
    name, val, (lo, hi); kind=:integer, dist=:discrete_ordered,
    dist_params=[lo, hi, step], label=label, description=description
)
_cots_discrete(name, val, (lo, hi), label, description) = ParamSpec(
    name, val, (lo, hi); kind=:integer, dist=:discrete_uniform, dist_params=[lo, hi],
    label=label, description=description
)

"""
Model factors of the crown-of-thorns starfish (COTS) submodel, run by COTSMod.jl.

`cots_on` switches the submodel on (1) or off (0, the default). With it off, ADRIA runs the
coral-only model and every other factor here is ignored.

The remaining factors are COTSMod parameters and opt-in COTSMod mechanisms. Mechanism
switches are categorical factors (0 = off, the default; 1 = on). Every factor is constant
at its value by default, so `sample()` does not vary it; the full bounds (broad
placeholders, not evidence-derived ranges) are kept for `set_factor_bounds!`. Calibrated
values can be supplied through a calibration NetCDF (`load_calib_params`) or by fixing the
factors directly.

Units follow COTSMod: population densities and the Allee threshold are COTS ha^-1, and
immigration is COTS ha^-1 per timestep.
"""
const COTS_PARAM_SPECS = ParamSpec[
    _cots_switch(:cots_on, 0,
        "COTS On",
        "Run the COTS submodel (1) or the coral-only model (0)."),

    # Core demographic parameters
    _cots_real(:cots_bh_a, 1.5, (0.5, 2.5),
        "COTS Beverton-Holt a",
        "Beverton-Holt recruitment parameter 'a'."),
    _cots_real(:cots_bh_b, 0.5, (0.1, 1.0),
        "COTS Beverton-Holt b",
        "Beverton-Holt density dependence 'b'."),
    _cots_real(:cots_ricker_a, 6.0, (2.0, 10.0),
        "COTS Ricker a",
        "Ricker recruitment maximum parameter."),
    _cots_real(:cots_ricker_b, 0.1, (0.01, 0.5),
        "COTS Ricker b",
        "Ricker density dependence parameter."),
    _cots_real(:cots_allee_threshold, 1.0, (0.1, 5.0),
        "Allee Threshold",
        "Allee half-saturation density (COTS ha^-1)."),
    _cots_switch(:cots_fecundity_gate, 0,
        "Fecundity Gate",
        "Gate fecundity on maternal condition (0 = off, 1 = on)."),

    # Mortality and starvation
    _cots_real(:cots_mortality_recruit, 0.4, (0.1, 0.9),
        "COTS Recruit Mortality",
        "Background mortality of age 0 (recruits)."),
    _cots_real(:cots_mortality_juvenile, 0.2, (0.05, 0.5),
        "COTS Juvenile Mortality",
        "Background mortality of age 1 (juveniles)."),
    _cots_real(:cots_mortality_adult, 0.1, (0.05, 0.3),
        "COTS Adult Mortality",
        "Background mortality of age 2+ (adults)."),
    _cots_real(:cots_starvation_max_fraction, 0.97, (0.8, 1.0),
        "Starvation Max Fraction",
        "Maximum fraction of mortality attributable to starvation."),
    _cots_real(:cots_max_food_cover, 0.8, (0.4, 1.0),
        "COTS Max Food Cover",
        "Coral cover giving maximum food availability."),
    _cots_real(:cots_starvation_exponent, 2.0, (1.0, 5.0),
        "Starvation Exponent",
        "Shape of the starvation threshold function."),
    _cots_real(:cots_condition_timescale, 5.0, (1.0, 10.0),
        "Condition Timescale",
        "Timescale (years) of the maternal-condition moving average."),

    # Predation (functional response)
    _cots_real(:cots_consumption_rate_fast, 0.6, (0.1, 1.0),
        "Fast Consumption Rate",
        "Consumption rate coefficient for fast-growing corals."),
    _cots_real(:cots_consumption_rate_slow, 0.15, (0.01, 0.5),
        "Slow Consumption Rate",
        "Consumption rate coefficient for slow-growing corals."),
    _cots_real(:cots_handling_time, 0.0, (0.0, 1.0),
        "Handling Time",
        "Predator handling time (Holling type II/III)."),
    _cots_real(:cots_consumption_exponent_fast, 1.0, (1.0, 3.0),
        "Fast Consumption Exponent",
        "Shape exponent for fast-coral consumption (type II/III " *
        "switch)."),
    _cots_real(:cots_consumption_exponent_slow, 1.0, (1.0, 3.0),
        "Slow Consumption Exponent",
        "Shape exponent for slow-coral consumption (type II/III " *
        "switch)."),

    # Immigration
    _cots_real(:cots_immigration_base, 0.002, (0.0, 0.01),
        "Base Immigration",
        "Background recruit immigration (COTS ha^-1 per timestep)."),
    _cots_real(:cots_immigration_threshold, 0.35, (0.1, 0.8),
        "Immigration Threshold",
        "Coral cover fraction below which immigration drops off."),
    _cots_real(:cots_immigration_exponent, 2.0, (1.0, 5.0),
        "Immigration Exponent",
        "Shape of the immigration gating function."),

    # Opt-in COTSMod mechanisms and their parameters
    _cots_switch(:cots_adult_senescence, 0,
        "Adult Senescence",
        "Extra mortality of adults older than the senescence age (0 = " *
        "off, 1 = on)."),
    _cots_ordered(:cots_senescence_age, 6, (3, 10, 1),
        "Senescence Age",
        "Adult age (years) from which senescence mortality applies."),
    _cots_real(:cots_senescence_mortality, 0.8, (0.5, 1.0),
        "Senescence Mortality",
        "Annual mortality of senescent adults."),
    _cots_switch(:cots_allee_stretched, 0,
        "Gradual Allee Effect",
        "Use the stretched (gradual) Allee function (0 = off, 1 = on)."),
    _cots_real(:cots_allee_exponent, 0.5, (0.25, 1.0),
        "Allee Exponent",
        "Exponent of the stretched Allee function."),
    _cots_switch(:cots_juvenile_predation, 0,
        "Juvenile Predation",
        "Saturating predation on juveniles, a predator pit (0 = off, 1 " *
        "= on)."),
    _cots_real(:cots_juvenile_predation_capacity, 20.0, (5.0, 50.0),
        "Predator Pit Capacity",
        "Juveniles predators can remove (COTS ha^-1 yr^-1)."),
    _cots_switch(:cots_per_capita_consumption, 0,
        "Per-capita Consumption",
        "Graze coral at a fixed area per adult per year (0 = off, 1 = " *
        "on)."),
    _cots_real(:cots_consumption_m2_per_adult_year, 10.0, (5.0, 20.0),
        "Consumption per Adult",
        "Coral area consumed per adult per year (m^2)."),
    _cots_switch(:cots_preferred_condition, 0,
        "Preferred-coral Condition",
        "Drive body condition by preferred (fast) coral (0 = off, 1 = " *
        "on)."),
    _cots_real(:cots_condition_saturation, 0.1, (0.05, 0.5),
        "Condition Saturation",
        "Preferred-coral cover fraction at which condition saturates."),
    _cots_real(:cots_condition_memory_years, 1.0, (0.5, 3.0),
        "Condition Memory",
        "Memory (years) of preferred-coral availability for condition."),
    _cots_switch(:cots_preferred_prey_starvation, 0,
        "Preferred-prey Starvation",
        "Drive starvation by remembered preferred coral (0 = off, 1 = " *
        "on)."),
    _cots_real(:cots_starvation_preferred_threshold, 0.05, (0.01, 0.2),
        "Starvation Threshold",
        "Preferred-coral cover fraction below which starvation applies."),
    _cots_real(:cots_starvation_memory_years, 0.0, (0.0, 3.0),
        "Starvation Memory",
        "Memory (years) of preferred-coral cover for starvation; 0 = " *
        "current year."),

    # Inputs from the data package (see `CotsInputs`); each is neutral by default
    _cots_switch(:cots_separated_recruitment, 0,
        "Separated Recruitment",
        "Settle larvae only through reef-scale dispersal, not locally " *
        "(0 = off, 1 = on). Needs annual reef connectivity in the data " *
        "package."),
    _cots_real(:cots_larval_fecundity, 6.0, (1.0, 60.0),
        "Larval Fecundity",
        "Larvae produced per adult before dispersal, with separated " *
        "recruitment."),
    _cots_levels(:cots_larval_conn_mode, 0, [0, 1, 2],
        "Larval Connectivity Mode",
        "COTS larval connectivity: 0 = mean, 1 = annual matrices in " *
        "turn, 2 = an annual matrix drawn each year."),
    _cots_discrete(:cots_larval_conn_seed, 0, (0, 1.0e9),
        "Larval Connectivity Seed",
        "Seed of the annual-matrix draws in mode 2."),
    _cots_real(:cots_external_outbreak_production, 0.0, (0.0, 50000.0),
        "External Outbreak Production",
        "Larval production of outbreaks outside the domain, scaling the " *
        "outbreak boundary in the data package. Needs separated " *
        "recruitment."),
    _cots_real(:cots_initial_state_scale, 1.0, (0.0, 10.0),
        "Initial State Scale",
        "Multiplier on the starting COTS densities in the data package.")
]

# Table-row view of `COTS_PARAM_SPECS`: per factor, the value, `ptype`, `dist`, full-range
# `dist_params`, `name` and `description` (what `_cots_spec` collapses to constants).
const COTS_FACTOR_SPECS = NamedTuple(
    fn => (;
        val=sp.val, ptype=sp.ptype, dist=sp.dist,
        dist_params=Tuple(sp.default_dist_params),
        name=sp.name, description=sp.description
    ) for (fn, sp) in spec_rows(COTS_PARAM_SPECS)
)

# Each factor is constant at its value by default (its sampling bounds collapse to the
# value, as `fix_factor!` does), so `sample()` leaves COTS off and adds no sampled
# dimensions. The full bounds are kept as `default_dist_params`; widen a factor's bounds
# with `set_factor_bounds!` to sample it.
function _cots_spec(fieldname::Symbol, val)::NamedTuple
    spec = getproperty(COTS_FACTOR_SPECS, fieldname)
    return _factor_spec(
        val;
        ptype=spec.ptype,
        dist=spec.dist,
        dist_params=Tuple(fill(Float64(val), length(spec.dist_params))),
        default_dist_params=spec.dist_params,
        name=spec.name,
        description=spec.description
    )
end
_cots_factor(fieldname::Symbol, val)::Param = Param(_cots_spec(fieldname, val))
_cots_factor(fieldname::Symbol)::Param = _cots_factor(
    fieldname, getproperty(COTS_FACTOR_SPECS, fieldname).val
)

Base.@kwdef struct CotsParams <: EcoModel
    cots_on::Param = _cots_factor(:cots_on)
    cots_bh_a::Param = _cots_factor(:cots_bh_a)
    cots_bh_b::Param = _cots_factor(:cots_bh_b)
    cots_ricker_a::Param = _cots_factor(:cots_ricker_a)
    cots_ricker_b::Param = _cots_factor(:cots_ricker_b)
    cots_allee_threshold::Param = _cots_factor(:cots_allee_threshold)
    cots_fecundity_gate::Param = _cots_factor(:cots_fecundity_gate)
    cots_mortality_recruit::Param = _cots_factor(:cots_mortality_recruit)
    cots_mortality_juvenile::Param = _cots_factor(:cots_mortality_juvenile)
    cots_mortality_adult::Param = _cots_factor(:cots_mortality_adult)
    cots_starvation_max_fraction::Param = _cots_factor(:cots_starvation_max_fraction)
    cots_max_food_cover::Param = _cots_factor(:cots_max_food_cover)
    cots_starvation_exponent::Param = _cots_factor(:cots_starvation_exponent)
    cots_condition_timescale::Param = _cots_factor(:cots_condition_timescale)
    cots_consumption_rate_fast::Param = _cots_factor(:cots_consumption_rate_fast)
    cots_consumption_rate_slow::Param = _cots_factor(:cots_consumption_rate_slow)
    cots_handling_time::Param = _cots_factor(:cots_handling_time)
    cots_consumption_exponent_fast::Param = _cots_factor(:cots_consumption_exponent_fast)
    cots_consumption_exponent_slow::Param = _cots_factor(:cots_consumption_exponent_slow)
    cots_immigration_base::Param = _cots_factor(:cots_immigration_base)
    cots_immigration_threshold::Param = _cots_factor(:cots_immigration_threshold)
    cots_immigration_exponent::Param = _cots_factor(:cots_immigration_exponent)
    cots_adult_senescence::Param = _cots_factor(:cots_adult_senescence)
    cots_senescence_age::Param = _cots_factor(:cots_senescence_age)
    cots_senescence_mortality::Param = _cots_factor(:cots_senescence_mortality)
    cots_allee_stretched::Param = _cots_factor(:cots_allee_stretched)
    cots_allee_exponent::Param = _cots_factor(:cots_allee_exponent)
    cots_juvenile_predation::Param = _cots_factor(:cots_juvenile_predation)
    cots_juvenile_predation_capacity::Param = _cots_factor(
        :cots_juvenile_predation_capacity
    )
    cots_per_capita_consumption::Param = _cots_factor(:cots_per_capita_consumption)
    cots_consumption_m2_per_adult_year::Param = _cots_factor(
        :cots_consumption_m2_per_adult_year
    )
    cots_preferred_condition::Param = _cots_factor(:cots_preferred_condition)
    cots_condition_saturation::Param = _cots_factor(:cots_condition_saturation)
    cots_condition_memory_years::Param = _cots_factor(:cots_condition_memory_years)
    cots_preferred_prey_starvation::Param = _cots_factor(:cots_preferred_prey_starvation)
    cots_starvation_preferred_threshold::Param = _cots_factor(
        :cots_starvation_preferred_threshold
    )
    cots_starvation_memory_years::Param = _cots_factor(:cots_starvation_memory_years)
    cots_separated_recruitment::Param = _cots_factor(:cots_separated_recruitment)
    cots_larval_fecundity::Param = _cots_factor(:cots_larval_fecundity)
    cots_larval_conn_mode::Param = _cots_factor(:cots_larval_conn_mode)
    cots_larval_conn_seed::Param = _cots_factor(:cots_larval_conn_seed)
    cots_external_outbreak_production::Param = _cots_factor(
        :cots_external_outbreak_production
    )
    cots_initial_state_scale::Param = _cots_factor(:cots_initial_state_scale)
end

# Names of the COTS factors. The derived per-scenario RNG seed (`set_random_seed`) excludes
# them, so adding the COTS component does not change the seeds of existing coral scenarios.
const COTS_FACTOR_NAMES = Set(string.(fieldnames(CotsParams)))

"""
    create_cots_params_instance(; overrides::Dict{String,Float64}=Dict{String,Float64}())::CotsParams

Build a `CotsParams` instance, replacing default factor values with any matching entry in
`overrides` (keyed by field name). Bounds and metadata are unaffected.
"""
function create_cots_params_instance(;
    overrides::Dict{String,Float64}=Dict{String,Float64}()
)::CotsParams
    return CotsParams(;
        (
            fn => _cots_factor(
                fn, get(overrides, string(fn), getproperty(COTS_FACTOR_SPECS, fn).val)
            )
            for fn in fieldnames(CotsParams)
        )...
    )
end

# ADRIA factor name => COTSMod.COTSParams field. Factors not listed here (`cots_on`, the
# larval connectivity mode and seed, external outbreak production and the initial state
# scale) are read by `run_model` directly. Switch factors (categorical) become `Bool`s, and
# ordered discrete factors become `Int`s, when COTSMod's parameters are built.
const COTS_PARAM_MAP = (
    cots_bh_a=:a,
    cots_bh_b=:b,
    cots_ricker_a=:a_ricker,
    cots_ricker_b=:b_ricker,
    cots_allee_threshold=:allee_threshold,
    cots_fecundity_gate=:fecundity_gate,
    cots_mortality_recruit=:m1,
    cots_mortality_juvenile=:m2,
    cots_mortality_adult=:m3,
    cots_starvation_max_fraction=:p_tilde,
    cots_max_food_cover=:C_max,
    cots_starvation_exponent=:eta_starve,
    cots_condition_timescale=:tau_condition,
    cots_consumption_rate_fast=:a_F,
    cots_consumption_rate_slow=:a_S,
    cots_handling_time=:h,
    cots_consumption_exponent_fast=:eta_F,
    cots_consumption_exponent_slow=:eta_S,
    cots_immigration_base=:IMM,
    cots_immigration_threshold=:imm_threshold,
    cots_immigration_exponent=:eta_imm,
    cots_adult_senescence=:adult_senescence,
    cots_senescence_age=:senescence_age,
    cots_senescence_mortality=:senescence_mortality,
    cots_allee_stretched=:allee_stretched,
    cots_allee_exponent=:allee_exponent,
    cots_juvenile_predation=:juvenile_predation,
    cots_juvenile_predation_capacity=:juvenile_predation_capacity,
    cots_per_capita_consumption=:per_capita_consumption,
    cots_consumption_m2_per_adult_year=:consumption_m2_per_adult_year,
    cots_preferred_condition=:preferred_condition,
    cots_condition_saturation=:condition_saturation,
    cots_condition_memory_years=:condition_memory_years,
    cots_preferred_prey_starvation=:preferred_prey_starvation,
    cots_starvation_preferred_threshold=:starvation_preferred_threshold,
    cots_starvation_memory_years=:starvation_memory_years,
    cots_separated_recruitment=:separated_recruitment,
    cots_larval_fecundity=:larval_fecundity
)

"""
    _cots_specs(; overrides)::Vector{Pair{Symbol,NamedTuple}}

Factor specs of every `CotsParams` field (in field order), with `overrides` (keyed by field
name) replacing the default values. The source of the `CotsParams` table rows.
"""
function _cots_specs(;
    overrides::Dict{String,Float64}=Dict{String,Float64}()
)::Vector{Pair{Symbol,NamedTuple}}
    return Pair{Symbol,NamedTuple}[
        fn => _cots_spec(
            fn, get(overrides, string(fn), getproperty(COTS_FACTOR_SPECS, fn).val)
        ) for fn in fieldnames(CotsParams)
    ]
end
