"""
    shading_intervention!(ctx::SimulationContext, tstep::Int64)

Process the shading intervention phase.
"""
function shading_intervention!(ctx::SimulationContext, tstep::Int64)
    if ctx.apply_shading && ctx.shade_decision_years[tstep]
        # Apply regional cooling effect
        ctx.shading_log[tstep, :] .= ctx.srm

        # Apply reduction in DHW due to SRM
        ctx.dhw_t .= max.(0.0, ctx.dhw_t .- ctx.srm)
    end
end