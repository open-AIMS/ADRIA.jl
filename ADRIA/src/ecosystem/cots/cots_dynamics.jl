# ADRIA adapter for the crown-of-thorns starfish (COTS) submodel in COTSMod.jl.
#
# ADRIA owns the parameter table, domain wiring, coral-cover tensors, result logging and
# scenario timing. COTSMod owns the COTS state transition, predation, larval dispersal and
# initialisation. Units follow COTSMod: COTS densities are COTS ha^-1.

"""
    cots_runtime_params(param_set::YAXArray)::COTSMod.COTSParams

Build COTSMod's runtime parameters from a scenario's `CotsParams` factors, following
`COTS_PARAM_MAP`. Mechanisms that are not exposed as factors keep COTSMod's defaults (off).
"""
function cots_runtime_params(param_set::YAXArray)::COTSMod.COTSParams
    kwargs = (
        field => _cots_param_value(factor, param_set[At(String(factor))])
        for (factor, field) in pairs(COTS_PARAM_MAP)
    )
    return COTSMod.COTSParams(; kwargs...)
end

# Convert a factor value to the type COTSMod expects, from the factor's `ptype`
function _cots_param_value(factor::Symbol, val)
    ptype = getproperty(COTS_FACTOR_SPECS, factor).ptype
    ptype == "unordered categorical" && return val == 1
    ptype == "ordered discrete" && return round(Int, val)
    return Float64(val)
end

"""
    cots_inputs(domain::Domain)::CotsInputs

The optional COTS inputs of `domain`. Domain types that do not carry any (e.g. ReefMod)
have none; `ADRIADomain` returns those loaded from its data package.
"""
cots_inputs(::Domain)::CotsInputs = CotsInputs()

"""
    cots_connectivity(domain::Domain)

The source-to-sink connectivity used for COTS larval dispersal: the data package's
`cots_connectivity/` when it has one (see `CotsInputs`), else the coral connectivity.
"""
function cots_connectivity(domain::Domain)
    conn = cots_inputs(domain).conn
    return isnothing(conn) ? domain.conn : conn
end

"""
    check_cots_settings(inputs::CotsInputs, separated::Bool, conn_mode::Integer, outbreak_production::Real)

Throw an `ErrorException` when a scenario with COTS on asks for something the domain's
`inputs` cannot supply. Call this for COTS-on scenarios only.
"""
function check_cots_settings(
    inputs::CotsInputs, separated::Bool, conn_mode::Integer, outbreak_production::Real
)::Nothing
    forcing = inputs.reef_forcing
    (separated || conn_mode != 0) && isnothing(forcing) &&
        error(
            "Separated recruitment and annual larval connectivity need annual reef " *
            "connectivity in the data package (cots_connectivity/annual_forcing/)"
        )
    outbreak_production > 0.0 && !separated &&
        error("External outbreak supply needs separated recruitment")
    outbreak_production > 0.0 && isempty(forcing.outbreak_years) &&
        error("External outbreak supply needs an outbreak boundary in the data package")
    return nothing
end

"""
    validate_cots_settings(domain::Domain, scens::DataFrame)::Nothing

Check every COTS-on scenario in `scens` against the inputs `domain` provides (see
[`check_cots_settings`](@ref)), so a bad setting fails before any scenario runs.
"""
function validate_cots_settings(domain::Domain, scens::DataFrame)::Nothing
    "cots_on" in names(scens) || return nothing
    inputs = cots_inputs(domain)
    for (i, row) in enumerate(eachrow(scens))
        row.cots_on == 1 || continue
        try
            check_cots_settings(
                inputs,
                row.cots_separated_recruitment == 1,
                round(Int64, row.cots_larval_conn_mode),
                row.cots_external_outbreak_production
            )
        catch err
            err isa ErrorException || rethrow()
            error("Scenario $i: $(err.msg)")
        end
    end
    return nothing
end

"""
Per-scenario COTS state and settings for `run_model`: the COTSMod population of every
location, the prey map and the settings read from the scenario's `CotsParams` factors.
Scenarios with COTS off use `nothing` in its place, for which [`cots_step!`](@ref) and
[`log_cots!`](@ref) do nothing.
"""
struct CotsRun
    models::COTSMod.COTSState
    prey_map::COTSMod.COTSPreyMap
    conn::SparseMatrixCSC{Float64,Int64}
    inputs::CotsInputs
    separated::Bool
    conn_mode::Int64
    conn_seed::Int64
    outbreak_production::Float64
    mean_reef_conn::Union{Nothing,Matrix{Float64}}
    habitable_fraction::Union{Nothing,Vector{Float64}}
end

"""
    CotsRun(domain::Domain, param_set::YAXArray, rng::AbstractRNG)

Initialise the COTS state for one scenario.
"""
function CotsRun(domain::Domain, param_set::YAXArray, rng::AbstractRNG)::CotsRun
    n_locs = domain.coral_growth.n_locs
    inputs = cots_inputs(domain)
    separated::Bool = param_set[At("cots_separated_recruitment")] == 1
    conn_mode = round(Int64, param_set[At("cots_larval_conn_mode")])
    outbreak_production::Float64 = param_set[At("cots_external_outbreak_production")]
    check_cots_settings(inputs, separated, conn_mode, outbreak_production)

    models = COTSMod.initialize_cots(n_locs, cots_runtime_params(param_set); rng=rng)
    if !isnothing(inputs.initial_state)
        COTSMod.set_cots_state!(
            models, param_set[At("cots_initial_state_scale")] .* inputs.initial_state
        )
    end

    forcing = inputs.reef_forcing
    # Reef connectivity averaged over the annual matrices (larval connectivity mode 0)
    mean_reef_conn =
        separated && conn_mode == 0 ?
        reduce(+, forcing.reef_connectivity) ./ length(forcing.reef_connectivity) :
        nothing
    # Per-capita grazing converts adults per hectare into habitable-cover fractions
    habitable_fraction =
        param_set[At("cots_per_capita_consumption")] == 1 ?
        Float64.(location_k(domain)) : nothing

    return CotsRun(
        models,
        # Fast prey: the three Acropora groups; slow prey: the two massive groups
        COTSMod.COTSPreyMap([1, 2, 3], [4, 5]),
        sparse(cots_connectivity(domain).data),
        inputs,
        separated,
        conn_mode,
        round(Int64, param_set[At("cots_larval_conn_seed")]),
        outbreak_production,
        mean_reef_conn,
        habitable_fraction
    )
end

"""
    cots_step!(cots, coral_cover, tstep::Int64, year::Int64)

Advance COTS one year: grazing on `coral_cover` (groups x sizes x locations), then larval
dispersal between locations and larvae arriving from outside the domain. `tstep` is the
timestep index and `year` its calendar year. Does nothing when `cots` is `nothing`.
"""
cots_step!(::Nothing, coral_cover, tstep::Int64, year::Int64) = nothing

function cots_step!(
    cots::CotsRun, coral_cover::AbstractArray{Float64,3}, tstep::Int64, year::Int64
)::Nothing
    models = cots.models
    forcing = cots.inputs.reef_forcing

    COTSMod.apply_predation!(
        coral_cover, models, cots.prey_map; habitable_fraction=cots.habitable_fraction
    )
    if cots.separated
        isnothing(cots.inputs.larval_survival) || COTSMod.scale_larval_production!(
            models, @view(cots.inputs.larval_survival[tstep, :])
        )
        isnothing(cots.inputs.settlement_modifier) ||
            COTSMod.scale_settlement!(models, cots.inputs.settlement_modifier)
    end

    forcing_idx = 0  # 0 = no annual matrix in use
    if cots.conn_mode == 0 && !cots.separated
        COTSMod.disperse_larvae!(models, cots.conn)
    elseif cots.conn_mode == 0
        COTSMod.disperse_larvae_by_reef!(
            models, cots.mean_reef_conn, forcing.site_to_reef;
            reef_areas=forcing.reef_habitable_area_ha
        )
    else
        forcing_idx = larval_conn_index(
            forcing, tstep - 1, cots.conn_mode, cots.conn_seed
        )
        COTSMod.disperse_larvae_by_reef!(
            models, forcing.reef_connectivity[forcing_idx], forcing.site_to_reef;
            reef_areas=cots.separated ? forcing.reef_habitable_area_ha : nothing
        )
    end

    if cots.outbreak_production > 0.0
        # The supply is passed as both the potential and the pelagic (post-survival)
        # supply, so pelagic survival is not applied again
        supply = outbreak_supply_by_location(
            forcing, forcing_idx, year, cots.outbreak_production
        )
        COTSMod.apply_external_larval_supply!(models, supply, supply)
    end
    return nothing
end

"""
    log_cots!(density_log, condition_log, cots, tstep::Int64)

Record each location's COTS densities (`density_log`, timesteps x 3 x locations) and body
condition (`condition_log`, timesteps x locations) for `tstep`. Does nothing when `cots` is
`nothing`.
"""
log_cots!(density_log, condition_log, ::Nothing, tstep::Int64) = nothing

function log_cots!(
    density_log::AbstractArray{<:Real,3}, condition_log::AbstractMatrix{<:Real},
    cots::CotsRun, tstep::Int64
)::Nothing
    for (loc, model) in enumerate(cots.models)
        density_log[tstep, :, loc] .= model.N
        condition_log[tstep, loc] = model.body_condition
    end
    return nothing
end

"""
    cots_log_specs(tf::Int, n_locs::Int)::Vector{LogSpec}
    cots_log_specs(domain::Domain)::Vector{LogSpec}

The output logs of the COTS submodel, as `LogSpec`s: `cots_pop_log` (COTS density by age
class, per hectare) and `cots_condition_log` (body condition). The result store names are
unchanged. Core function until COTS is a submodel that declares its own `log_specs`.
"""
function cots_log_specs(tf::Int, n_locs::Int)::Vector{ADRIAinterfaces.LogSpec}
    return [
        ADRIAinterfaces.LogSpec(
            :cots_pop_log, [:timesteps, :age_class, :locations], [tf, 3, n_locs];
            labels=Dict(:age_class => ["recruits", "juveniles", "adults"]),
            units="COTS ha^-1"
        ),
        ADRIAinterfaces.LogSpec(
            :cots_condition_log, [:timesteps, :locations], [tf, n_locs]
        )
    ]
end
function cots_log_specs(domain::Domain)::Vector{ADRIAinterfaces.LogSpec}
    return cots_log_specs(length(timesteps(domain)), domain.coral_growth.n_locs)
end

"""
    larval_conn_index(forcing::CotsReefForcing, timestep::Int, mode::Int, seed::Int)::Int

Index of the annual reef matrix used at `timestep`: in turn (`mode` 1) or drawn with a
generator seeded by `seed + timestep` (`mode` 2), so the draws do not depend on other
random numbers in the scenario.
"""
function larval_conn_index(
    forcing::CotsReefForcing, timestep::Int, mode::Int, seed::Int
)::Int
    mode == 1 && return mod1(timestep, length(forcing.years))
    mode == 2 && return rand(MersenneTwister(seed + timestep), eachindex(forcing.years))
    throw(ArgumentError("Larval connectivity mode must be 1 or 2 for annual matrices"))
end

"""
    outbreak_supply_by_location(forcing, forcing_idx, calendar_year, production)::Vector{Float64}

Larvae arriving at each location from outbreaks outside the domain in `calendar_year`:
`production` times the reef's outbreak-boundary coefficient for that year, under annual
matrix `forcing_idx` (0 = the mean over the annual matrices). Zero in years the boundary
does not cover.
"""
function outbreak_supply_by_location(
    forcing::CotsReefForcing, forcing_idx::Int, calendar_year::Int, production::Float64
)::Vector{Float64}
    year_idx = findfirst(==(calendar_year), forcing.outbreak_years)
    isnothing(year_idx) && return zeros(Float64, length(forcing.site_to_reef))
    coefficients = if forcing_idx == 0
        vec(mean(@view(forcing.outbreak_supply[:, year_idx, :]); dims=2))
    else
        @view(forcing.outbreak_supply[:, year_idx, forcing_idx])
    end
    return [production * coefficients[reef] for reef in forcing.site_to_reef]
end
