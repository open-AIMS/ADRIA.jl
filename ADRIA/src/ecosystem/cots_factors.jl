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
    a=(val=1.5, ptype="continuous", dist=TriangularDist, dist_params=(0.5, 2.5, 1.5),
        name="BH a", description="Beverton-Holt recruitment parameter 'a'."),
    b=(val=0.5, ptype="continuous", dist=TriangularDist, dist_params=(0.1, 1.0, 0.5),
        name="BH b", description="Beverton-Holt density dependence 'b'."),
    a_ricker=(val=6.0, ptype="continuous", dist=TriangularDist,
        dist_params=(2.0, 10.0, 6.0),
        name="Ricker a", description="Ricker recruitment maximum parameter."),
    b_ricker=(val=0.1, ptype="continuous", dist=TriangularDist,
        dist_params=(0.01, 0.5, 0.1),
        name="Ricker b", description="Ricker density dependence parameter."),
    allee_threshold=(val=1.0, ptype="continuous", dist=TriangularDist,
        dist_params=(0.1, 5.0, 1.0), name="Allee Threshold",
        description="Allee half-saturation density (COTS ha^-1)."),
    fecundity_gate=(val=0, ptype="unordered categorical", dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Fecundity Gate",
        description="Gate fecundity on maternal condition (0 = off, 1 = on)."),

    # Mortality and starvation
    m1=(val=0.4, ptype="continuous", dist=TriangularDist, dist_params=(0.1, 0.9, 0.4),
        name="m1", description="Background mortality of age 0 (recruits)."),
    m2=(val=0.2, ptype="continuous", dist=TriangularDist, dist_params=(0.05, 0.5, 0.2),
        name="m2", description="Background mortality of age 1 (juveniles)."),
    m3=(val=0.1, ptype="continuous", dist=TriangularDist, dist_params=(0.05, 0.3, 0.1),
        name="m3", description="Background mortality of age 2+ (adults)."),
    p_tilde=(val=0.97, ptype="continuous", dist=TriangularDist,
        dist_params=(0.8, 1.0, 0.97),
        name="Starvation Max Fraction",
        description="Maximum fraction of mortality attributable to starvation."),
    C_max=(val=0.8, ptype="continuous", dist=TriangularDist, dist_params=(0.4, 1.0, 0.8),
        name="C_max", description="Coral cover giving maximum food availability."),
    eta_starve=(val=2.0, ptype="continuous", dist=TriangularDist,
        dist_params=(1.0, 5.0, 2.0),
        name="Starvation Exponent",
        description="Shape of the starvation threshold function."),
    tau_condition=(val=5.0, ptype="continuous", dist=TriangularDist,
        dist_params=(1.0, 10.0, 5.0), name="Condition Timescale",
        description="Timescale (years) of the maternal-condition moving average."),

    # Predation (functional response)
    a_F=(val=0.6, ptype="continuous", dist=TriangularDist, dist_params=(0.1, 1.0, 0.6),
        name="Fast Consumption Rate",
        description="Consumption rate coefficient for fast-growing corals."),
    a_S=(val=0.15, ptype="continuous", dist=TriangularDist, dist_params=(0.01, 0.5, 0.15),
        name="Slow Consumption Rate",
        description="Consumption rate coefficient for slow-growing corals."),
    h=(val=0.0, ptype="continuous", dist=TriangularDist, dist_params=(0.0, 1.0, 0.0),
        name="Handling Time", description="Predator handling time (Holling type II/III)."),
    eta_F=(val=1.0, ptype="continuous", dist=TriangularDist, dist_params=(1.0, 3.0, 1.0),
        name="Fast Consumption Exponent",
        description="Shape exponent for fast-coral consumption (type II/III switch)."),
    eta_S=(val=1.0, ptype="continuous", dist=TriangularDist, dist_params=(1.0, 3.0, 1.0),
        name="Slow Consumption Exponent",
        description="Shape exponent for slow-coral consumption (type II/III switch)."),

    # Immigration
    IMM=(val=0.002, ptype="continuous", dist=TriangularDist, dist_params=(0.0, 0.01, 0.002),
        name="Base Immigration",
        description="Background recruit immigration (COTS ha^-1 per timestep)."),
    imm_threshold=(val=0.35, ptype="continuous", dist=TriangularDist,
        dist_params=(0.1, 0.8, 0.35), name="Immigration Threshold",
        description="Coral cover fraction below which immigration drops off."),
    eta_imm=(val=2.0, ptype="continuous", dist=TriangularDist, dist_params=(1.0, 5.0, 2.0),
        name="Immigration Exponent", description="Shape of the immigration gating function."
    ),

    # Opt-in COTSMod mechanisms and their parameters
    adult_senescence=(val=0, ptype="unordered categorical", dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Adult Senescence",
        description="Extra mortality of adults older than the senescence age (0 = off, 1 = on)."
    ),
    senescence_age=(val=6, ptype="ordered discrete", dist=DiscreteOrderedUniformDist,
        dist_params=(3.0, 10.0, 1.0), name="Senescence Age",
        description="Adult age (years) from which senescence mortality applies."),
    senescence_mortality=(val=0.8, ptype="continuous", dist=TriangularDist,
        dist_params=(0.5, 1.0, 0.8), name="Senescence Mortality",
        description="Annual mortality of senescent adults."),
    allee_stretched=(val=0, ptype="unordered categorical", dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Gradual Allee Effect",
        description="Use the stretched (gradual) Allee function (0 = off, 1 = on)."),
    allee_exponent=(val=0.5, ptype="continuous", dist=TriangularDist,
        dist_params=(0.25, 1.0, 0.5), name="Allee Exponent",
        description="Exponent of the stretched Allee function."),
    juvenile_predation=(val=0, ptype="unordered categorical", dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Juvenile Predation",
        description="Saturating predation on juveniles, a predator pit (0 = off, 1 = on)."),
    juvenile_predation_capacity=(val=20.0, ptype="continuous", dist=TriangularDist,
        dist_params=(5.0, 50.0, 20.0), name="Predator Pit Capacity",
        description="Juveniles predators can remove (COTS ha^-1 yr^-1)."),
    per_capita_consumption=(val=0, ptype="unordered categorical",
        dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Per-capita Consumption",
        description="Graze coral at a fixed area per adult per year (0 = off, 1 = on)."),
    consumption_m2_per_adult_year=(val=10.0, ptype="continuous", dist=TriangularDist,
        dist_params=(5.0, 20.0, 10.0), name="Consumption per Adult",
        description="Coral area consumed per adult per year (m^2)."),
    preferred_condition=(val=0, ptype="unordered categorical", dist=CategoricalDistribution,
        dist_params=(0.0, 1.0), name="Preferred-coral Condition",
        description="Drive body condition by preferred (fast) coral (0 = off, 1 = on)."),
    condition_saturation=(val=0.10, ptype="continuous", dist=TriangularDist,
        dist_params=(0.05, 0.5, 0.10), name="Condition Saturation",
        description="Preferred-coral cover fraction at which condition saturates."),
    condition_memory_years=(val=1.0, ptype="continuous", dist=TriangularDist,
        dist_params=(0.5, 3.0, 1.0), name="Condition Memory",
        description="Memory (years) of preferred-coral availability for condition."),
    preferred_prey_starvation=(val=0, ptype="unordered categorical",
        dist=CategoricalDistribution, dist_params=(0.0, 1.0),
        name="Preferred-prey Starvation",
        description="Drive starvation by remembered preferred coral (0 = off, 1 = on)."),
    starvation_preferred_threshold=(val=0.05, ptype="continuous", dist=TriangularDist,
        dist_params=(0.01, 0.2, 0.05), name="Starvation Threshold",
        description="Preferred-coral cover fraction below which starvation applies."),
    starvation_memory_years=(val=0.0, ptype="continuous", dist=TriangularDist,
        dist_params=(0.0, 3.0, 0.0), name="Starvation Memory",
        description="Memory (years) of preferred-coral cover for starvation; 0 = current year."
    )
)

# Each factor is constant at its value by default (its sampling bounds collapse to the
# value, as `fix_factor!` does), so `sample()` leaves COTS off and adds no sampled
# dimensions. The full bounds are kept as `default_dist_params`; widen a factor's bounds
# with `set_factor_bounds!` to sample it.
function _cots_factor(fieldname::Symbol, val)::Param
    spec = getproperty(COTS_FACTOR_SPECS, fieldname)
    return Factor(
        val;
        ptype=spec.ptype,
        dist=spec.dist,
        dist_params=Tuple(fill(Float64(val), length(spec.dist_params))),
        default_dist_params=spec.dist_params,
        name=spec.name,
        description=spec.description
    )
end
_cots_factor(fieldname::Symbol)::Param = _cots_factor(
    fieldname, getproperty(COTS_FACTOR_SPECS, fieldname).val
)

Base.@kwdef struct CotsParams <: EcoModel
    cots_on::Param = _cots_factor(:cots_on)
    a::Param = _cots_factor(:a)
    b::Param = _cots_factor(:b)
    a_ricker::Param = _cots_factor(:a_ricker)
    b_ricker::Param = _cots_factor(:b_ricker)
    allee_threshold::Param = _cots_factor(:allee_threshold)
    fecundity_gate::Param = _cots_factor(:fecundity_gate)
    m1::Param = _cots_factor(:m1)
    m2::Param = _cots_factor(:m2)
    m3::Param = _cots_factor(:m3)
    p_tilde::Param = _cots_factor(:p_tilde)
    C_max::Param = _cots_factor(:C_max)
    eta_starve::Param = _cots_factor(:eta_starve)
    tau_condition::Param = _cots_factor(:tau_condition)
    a_F::Param = _cots_factor(:a_F)
    a_S::Param = _cots_factor(:a_S)
    h::Param = _cots_factor(:h)
    eta_F::Param = _cots_factor(:eta_F)
    eta_S::Param = _cots_factor(:eta_S)
    IMM::Param = _cots_factor(:IMM)
    imm_threshold::Param = _cots_factor(:imm_threshold)
    eta_imm::Param = _cots_factor(:eta_imm)
    adult_senescence::Param = _cots_factor(:adult_senescence)
    senescence_age::Param = _cots_factor(:senescence_age)
    senescence_mortality::Param = _cots_factor(:senescence_mortality)
    allee_stretched::Param = _cots_factor(:allee_stretched)
    allee_exponent::Param = _cots_factor(:allee_exponent)
    juvenile_predation::Param = _cots_factor(:juvenile_predation)
    juvenile_predation_capacity::Param = _cots_factor(:juvenile_predation_capacity)
    per_capita_consumption::Param = _cots_factor(:per_capita_consumption)
    consumption_m2_per_adult_year::Param = _cots_factor(:consumption_m2_per_adult_year)
    preferred_condition::Param = _cots_factor(:preferred_condition)
    condition_saturation::Param = _cots_factor(:condition_saturation)
    condition_memory_years::Param = _cots_factor(:condition_memory_years)
    preferred_prey_starvation::Param = _cots_factor(:preferred_prey_starvation)
    starvation_preferred_threshold::Param = _cots_factor(:starvation_preferred_threshold)
    starvation_memory_years::Param = _cots_factor(:starvation_memory_years)
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
