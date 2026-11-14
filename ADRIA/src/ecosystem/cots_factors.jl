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
const COTS_FACTOR_SPECS = (
    cots_on=(val=0, ptype="unordered categorical", dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="COTS On",
        description="Run the COTS submodel (1) or the coral-only model (0)."),

    # Core demographic parameters
    cots_bh_a=(val=1.5, ptype="continuous", dist=TriangularDist,
        dist_params=(0.5, 2.5, 1.5),
        name="COTS Beverton-Holt a", description="Beverton-Holt recruitment parameter 'a'."
    ),
    cots_bh_b=(val=0.5, ptype="continuous", dist=TriangularDist,
        dist_params=(0.1, 1.0, 0.5),
        name="COTS Beverton-Holt b", description="Beverton-Holt density dependence 'b'."),
    cots_ricker_a=(val=6.0, ptype="continuous", dist=TriangularDist,
        dist_params=(2.0, 10.0, 6.0),
        name="COTS Ricker a", description="Ricker recruitment maximum parameter."),
    cots_ricker_b=(val=0.1, ptype="continuous", dist=TriangularDist,
        dist_params=(0.01, 0.5, 0.1),
        name="COTS Ricker b", description="Ricker density dependence parameter."),
    cots_allee_threshold=(val=1.0, ptype="continuous", dist=TriangularDist,
        dist_params=(0.1, 5.0, 1.0), name="Allee Threshold",
        description="Allee half-saturation density (COTS ha^-1)."),
    cots_fecundity_gate=(val=0, ptype="unordered categorical", dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Fecundity Gate",
        description="Gate fecundity on maternal condition (0 = off, 1 = on)."),

    # Mortality and starvation
    cots_mortality_recruit=(val=0.4, ptype="continuous", dist=TriangularDist,
        dist_params=(0.1, 0.9, 0.4),
        name="COTS Recruit Mortality",
        description="Background mortality of age 0 (recruits)."),
    cots_mortality_juvenile=(val=0.2, ptype="continuous", dist=TriangularDist,
        dist_params=(0.05, 0.5, 0.2),
        name="COTS Juvenile Mortality",
        description="Background mortality of age 1 (juveniles)."),
    cots_mortality_adult=(val=0.1, ptype="continuous", dist=TriangularDist,
        dist_params=(0.05, 0.3, 0.1),
        name="COTS Adult Mortality", description="Background mortality of age 2+ (adults)."
    ),
    cots_starvation_max_fraction=(val=0.97, ptype="continuous", dist=TriangularDist,
        dist_params=(0.8, 1.0, 0.97),
        name="Starvation Max Fraction",
        description="Maximum fraction of mortality attributable to starvation."),
    cots_max_food_cover=(val=0.8, ptype="continuous", dist=TriangularDist,
        dist_params=(0.4, 1.0, 0.8),
        name="COTS Max Food Cover",
        description="Coral cover giving maximum food availability."),
    cots_starvation_exponent=(val=2.0, ptype="continuous", dist=TriangularDist,
        dist_params=(1.0, 5.0, 2.0),
        name="Starvation Exponent",
        description="Shape of the starvation threshold function."),
    cots_condition_timescale=(val=5.0, ptype="continuous", dist=TriangularDist,
        dist_params=(1.0, 10.0, 5.0), name="Condition Timescale",
        description="Timescale (years) of the maternal-condition moving average."),

    # Predation (functional response)
    cots_consumption_rate_fast=(val=0.6, ptype="continuous", dist=TriangularDist,
        dist_params=(0.1, 1.0, 0.6),
        name="Fast Consumption Rate",
        description="Consumption rate coefficient for fast-growing corals."),
    cots_consumption_rate_slow=(val=0.15, ptype="continuous", dist=TriangularDist,
        dist_params=(0.01, 0.5, 0.15),
        name="Slow Consumption Rate",
        description="Consumption rate coefficient for slow-growing corals."),
    cots_handling_time=(val=0.0, ptype="continuous", dist=TriangularDist,
        dist_params=(0.0, 1.0, 0.0),
        name="Handling Time", description="Predator handling time (Holling type II/III)."),
    cots_consumption_exponent_fast=(val=1.0, ptype="continuous", dist=TriangularDist,
        dist_params=(1.0, 3.0, 1.0),
        name="Fast Consumption Exponent",
        description="Shape exponent for fast-coral consumption (type II/III switch)."),
    cots_consumption_exponent_slow=(val=1.0, ptype="continuous", dist=TriangularDist,
        dist_params=(1.0, 3.0, 1.0),
        name="Slow Consumption Exponent",
        description="Shape exponent for slow-coral consumption (type II/III switch)."),

    # Immigration
    cots_immigration_base=(val=0.002, ptype="continuous", dist=TriangularDist,
        dist_params=(0.0, 0.01, 0.002),
        name="Base Immigration",
        description="Background recruit immigration (COTS ha^-1 per timestep)."),
    cots_immigration_threshold=(val=0.35, ptype="continuous", dist=TriangularDist,
        dist_params=(0.1, 0.8, 0.35), name="Immigration Threshold",
        description="Coral cover fraction below which immigration drops off."),
    cots_immigration_exponent=(val=2.0, ptype="continuous", dist=TriangularDist,
        dist_params=(1.0, 5.0, 2.0),
        name="Immigration Exponent", description="Shape of the immigration gating function."
    ),

    # Opt-in COTSMod mechanisms and their parameters
    cots_adult_senescence=(val=0, ptype="unordered categorical",
        dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Adult Senescence",
        description="Extra mortality of adults older than the senescence age (0 = off, 1 = on)."
    ),
    cots_senescence_age=(val=6, ptype="ordered discrete", dist=DiscreteOrderedUniformDist,
        dist_params=(3.0, 10.0, 1.0), name="Senescence Age",
        description="Adult age (years) from which senescence mortality applies."),
    cots_senescence_mortality=(val=0.8, ptype="continuous", dist=TriangularDist,
        dist_params=(0.5, 1.0, 0.8), name="Senescence Mortality",
        description="Annual mortality of senescent adults."),
    cots_allee_stretched=(val=0, ptype="unordered categorical",
        dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Gradual Allee Effect",
        description="Use the stretched (gradual) Allee function (0 = off, 1 = on)."),
    cots_allee_exponent=(val=0.5, ptype="continuous", dist=TriangularDist,
        dist_params=(0.25, 1.0, 0.5), name="Allee Exponent",
        description="Exponent of the stretched Allee function."),
    cots_juvenile_predation=(val=0, ptype="unordered categorical",
        dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Juvenile Predation",
        description="Saturating predation on juveniles, a predator pit (0 = off, 1 = on)."),
    cots_juvenile_predation_capacity=(val=20.0, ptype="continuous", dist=TriangularDist,
        dist_params=(5.0, 50.0, 20.0), name="Predator Pit Capacity",
        description="Juveniles predators can remove (COTS ha^-1 yr^-1)."),
    cots_per_capita_consumption=(val=0, ptype="unordered categorical",
        dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Per-capita Consumption",
        description="Graze coral at a fixed area per adult per year (0 = off, 1 = on)."),
    cots_consumption_m2_per_adult_year=(val=10.0, ptype="continuous", dist=TriangularDist,
        dist_params=(5.0, 20.0, 10.0), name="Consumption per Adult",
        description="Coral area consumed per adult per year (m^2)."),
    cots_preferred_condition=(val=0, ptype="unordered categorical",
        dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Preferred-coral Condition",
        description="Drive body condition by preferred (fast) coral (0 = off, 1 = on)."),
    cots_condition_saturation=(val=0.10, ptype="continuous", dist=TriangularDist,
        dist_params=(0.05, 0.5, 0.10), name="Condition Saturation",
        description="Preferred-coral cover fraction at which condition saturates."),
    cots_condition_memory_years=(val=1.0, ptype="continuous", dist=TriangularDist,
        dist_params=(0.5, 3.0, 1.0), name="Condition Memory",
        description="Memory (years) of preferred-coral availability for condition."),
    cots_preferred_prey_starvation=(val=0, ptype="unordered categorical",
        dist=CategoricalDistribution, dist_params=(0.0, 1.0),
        name="Preferred-prey Starvation",
        description="Drive starvation by remembered preferred coral (0 = off, 1 = on)."),
    cots_starvation_preferred_threshold=(val=0.05, ptype="continuous", dist=TriangularDist,
        dist_params=(0.01, 0.2, 0.05), name="Starvation Threshold",
        description="Preferred-coral cover fraction below which starvation applies."),
    cots_starvation_memory_years=(val=0.0, ptype="continuous", dist=TriangularDist,
        dist_params=(0.0, 3.0, 0.0), name="Starvation Memory",
        description="Memory (years) of preferred-coral cover for starvation; 0 = current year."
    ),

    # Inputs from the data package (see `CotsInputs`); each is neutral by default
    cots_separated_recruitment=(val=0, ptype="unordered categorical",
        dist=CategoricalDistribution, dist_params=(0.0, 1.0),
        name="Separated Recruitment",
        description="Settle larvae only through reef-scale dispersal, not locally " *
                    "(0 = off, 1 = on). Needs annual reef connectivity in the data package."
    ),
    cots_larval_fecundity=(val=6.0, ptype="continuous", dist=TriangularDist,
        dist_params=(1.0, 60.0, 6.0), name="Larval Fecundity",
        description="Larvae produced per adult before dispersal, with separated recruitment."
    ),
    cots_larval_conn_mode=(val=0, ptype="unordered categorical",
        dist=CategoricalDistribution,
        dist_params=(0.0, 1.0, 2.0), name="Larval Connectivity Mode",
        description="COTS larval connectivity: 0 = mean, 1 = annual matrices in turn, " *
                    "2 = an annual matrix drawn each year."),
    cots_larval_conn_seed=(val=0, ptype="ordered discrete", dist=DiscreteUniform,
        dist_params=(0.0, 1.0e9), name="Larval Connectivity Seed",
        description="Seed of the annual-matrix draws in mode 2."),
    cots_external_outbreak_production=(val=0.0, ptype="continuous", dist=TriangularDist,
        dist_params=(0.0, 50000.0, 0.0), name="External Outbreak Production",
        description="Larval production of outbreaks outside the domain, scaling the " *
                    "outbreak boundary in the data package. Needs separated recruitment."),
    cots_initial_state_scale=(val=1.0, ptype="continuous", dist=TriangularDist,
        dist_params=(0.0, 10.0, 1.0), name="Initial State Scale",
        description="Multiplier on the starting COTS densities in the data package.")
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
