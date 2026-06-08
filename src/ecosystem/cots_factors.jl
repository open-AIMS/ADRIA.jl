Base.@kwdef struct CotsParams <: EcoModel
    # Core demographic parameters
    a::Param = Factor(1.5; ptype="continuous", dist=TriangularDist, dist_params=(0.5, 2.5, 1.5), name="BH a", description="Beverton-Holt recruitment parameter 'a'")
    b::Param = Factor(0.5; ptype="continuous", dist=TriangularDist, dist_params=(0.1, 1.0, 0.5), name="BH b", description="Beverton-Holt density dependence 'b'")
    
    # Ricker recruitment parameters
    a_ricker::Param = Factor(6.0; ptype="continuous", dist=TriangularDist, dist_params=(2.0, 10.0, 6.0), name="Ricker a", description="Ricker recruitment max parameter")
    b_ricker::Param = Factor(0.1; ptype="continuous", dist=TriangularDist, dist_params=(0.01, 0.5, 0.1), name="Ricker b", description="Ricker density dependence parameter")
    allee_threshold::Param = Factor(1.0; ptype="continuous", dist=TriangularDist, dist_params=(0.1, 5.0, 1.0), name="Allee Threshold", description="Allee effect threshold (population size where fertilization halves)")
    fecundity_gate::Param = Factor(0; ptype="unordered categorical", dist=CategoricalDistribution, dist_params=(0.0, 1.0), name="Fecundity Gate", description="Toggle maternal condition gating on fecundity (0=off, 1=on)")
    
    # Mortality and starvation
    m1::Param = Factor(0.4; ptype="continuous", dist=TriangularDist, dist_params=(0.1, 0.9, 0.4), name="m1", description="Background mortality of Age 0 (recruits)")
    m2::Param = Factor(0.2; ptype="continuous", dist=TriangularDist, dist_params=(0.05, 0.5, 0.2), name="m2", description="Background mortality of Age 1 (juveniles)")
    m3::Param = Factor(0.1; ptype="continuous", dist=TriangularDist, dist_params=(0.05, 0.3, 0.1), name="m3", description="Background mortality of Age 2+ (adults)")
    p_tilde::Param = Factor(0.97; ptype="continuous", dist=TriangularDist, dist_params=(0.8, 1.0, 0.97), name="Starvation max fraction", description="Maximum fraction of mortality attributable to starvation")
    C_max::Param = Factor(0.8; ptype="continuous", dist=TriangularDist, dist_params=(0.4, 1.0, 0.8), name="C_max", description="Coral cover level corresponding to maximum food availability")
    eta_starve::Param = Factor(2.0; ptype="continuous", dist=TriangularDist, dist_params=(1.0, 5.0, 2.0), name="Starvation Exponent", description="Shape parameter for starvation threshold function")
    tau_condition::Param = Factor(5.0; ptype="continuous", dist=TriangularDist, dist_params=(1.0, 10.0, 5.0), name="Condition Timescale", description="Timescale (years) for maternal condition exponential moving average")

    # Predation (Functional Response)
    a_F::Param = Factor(0.6; ptype="continuous", dist=TriangularDist, dist_params=(0.1, 1.0, 0.6), name="Fast Consumption Rate", description="Base consumption rate coefficient for fast-growing corals")
    a_S::Param = Factor(0.15; ptype="continuous", dist=TriangularDist, dist_params=(0.01, 0.5, 0.15), name="Slow Consumption Rate", description="Base consumption rate coefficient for slow-growing corals")
    h::Param = Factor(0.0; ptype="continuous", dist=TriangularDist, dist_params=(0.0, 1.0, 0.0), name="Handling Time", description="Predator handling time (for Holling Type II/III)")
    eta_F::Param = Factor(1.0; ptype="continuous", dist=TriangularDist, dist_params=(1.0, 3.0, 1.0), name="Fast Consumption Exponent", description="Shape exponent for fast coral consumption (Type II/III switch)")
    eta_S::Param = Factor(1.0; ptype="continuous", dist=TriangularDist, dist_params=(1.0, 3.0, 1.0), name="Slow Consumption Exponent", description="Shape exponent for slow coral consumption (Type II/III switch)")

    # Dispersal / Immigration
    IMM::Param = Factor(0.002; ptype="continuous", dist=TriangularDist, dist_params=(0.0, 0.01, 0.002), name="Base Immigration", description="Base background immigration rate per timestep")
    imm_threshold::Param = Factor(0.35; ptype="continuous", dist=TriangularDist, dist_params=(0.1, 0.8, 0.35), name="Immigration Threshold", description="Coral cover fraction below which immigration drops off")
    eta_imm::Param = Factor(2.0; ptype="continuous", dist=TriangularDist, dist_params=(1.0, 5.0, 2.0), name="Immigration Exponent", description="Shape parameter for immigration gating function")
end
