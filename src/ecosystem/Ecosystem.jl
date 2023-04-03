using Setfield

using DataFrames
using ModelParameters
import ModelParameters: update!, Model


abstract type EcoModel end


"""Set a model parameter value directly."""
function set(p::Param, val::Union{Int64,Float64})
    if hasproperty(p, :ptype) && p.ptype == "integer" && !isinteger(val)
        val = map_to_discrete(val, p.bounds[2])
    end

    return val
end


"""
    map_to_discrete(v::Union{Int64,Float64}, u::Union{Int64,Float64})::Int64

For integer/categorical parameters, take floor of `v`, capping to `u - 1`
"""
function map_to_discrete(v::Union{Int64,Float64}, u::Union{Int64,Float64})::Int64
    return Int64(min(floor(v), u - 1))
end

"""
    map_to_discrete!(df::DataFrame, u::AbstractArray)::Nothing

Update a dataframe of parameters.
Length of `u` is expected to match number of columns in `df`.
"""
function map_to_discrete!(df::DataFrame, u::Union{AbstractVector,Tuple})::Nothing
    for (idx, b) in enumerate(u)
        df[!, idx] .= map_to_discrete.(df[!, idx], b)
    end
end


Base.@kwdef struct Intervention{N,P,N2,P2} <: EcoModel
    # Intervention Parameters
    # Integer values have a +1 offset to allow for discrete value mapping
    # (see `set()` and `map_to_discrete()` methods)
    # Bounds are defined as floats to maintain type stability
    guided::N = Param(0, ptype="integer", bounds=(-1.0, 3.0 + 1.0), dists="unif",
        name="Guided", description="Choice of MCDA approach.")
    seed_TA::N = Param(0, ptype="integer", bounds=(0.0, 1000000.0 + 1.0), dists="unif",
        name="Seeded Tabular Acropora", description="Number of enhanced Tabular Acropora to seed per deployment year.")
    seed_CA::N = Param(0, ptype="integer", bounds=(0.0, 1000000.0 + 1.0), dists="unif",
        name="Seeded Corymbose Acropora", description="Number of enhanced Corymbose Acropora to seed per deployment year.")
    fogging::P = Param(0.16, ptype="real", bounds=(0.0, 0.3, 0.16 / 0.3), dists="triang",
        name="Fogging", description="Assumed reduction in bleaching mortality.")
    SRM::P = Param(0.0, ptype="real", bounds=(0.0, 7.0, 0.0), dists="triang",
        name="SRM", description="Reduction in DHWs due to shading.")
    a_adapt::P = Param(0.0, ptype="real", bounds=(0.0, 8.0, 0.0), dists="triang",
        name="Assisted Adaptation", description="Assisted adaptation in terms of DHW resistance.")
    n_adapt::N2 = Param(0.0, ptype="real", bounds=(0.0, 0.05), dists="unif",
        name="Natural Adaptation", description="Natural adaptation rate (yearly increase).")
    seed_years::P2 = Param(10, ptype="integer", bounds=(5.0, 74.0 + 1.0, 5 / 70), dists="triang",
        name="Years to Seed", description="Number of years to seed for.")
    shade_years::P2 = Param(10, ptype="integer", bounds=(5.0, 74.0 + 1.0, 5 / 70), dists="triang",
        name="Years to Shade", description="Number of years to shade for.")
    seed_freq::N = Param(5, ptype="integer", bounds=(0.0, 5.0 + 1.0), dists="unif",
        name="Seeding Frequency", description="Frequency of seeding site selection (0 is set and forget).")
    shade_freq::N = Param(1, ptype="integer", bounds=(0.0, 5.0 + 1.0), dists="unif",
        name="Shading Frequency", description="Frequency of shading site selection (0 is set and forget).")
    seed_year_start::N = Param(2, ptype="integer", bounds=(2.0, 25.0 + 1.0), dists="unif",
        name="Seeding Start Year", description="Start seeding deployments after this number of years has elapsed.")
    shade_year_start::N = Param(2, ptype="integer", bounds=(2.0, 25.0 + 1.0), dists="unif",
        name="Shading Start Year", description="Start of shading deployments after this number of years has elapsed.")
end


Base.@kwdef struct Criteria{P,N} <: EcoModel
    wave_stress::P = Param(1.0, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Wave Stress", description="Importance of avoiding wave stress. Higher values places more weight on areas with low wave stress.")
    heat_stress::P = Param(1.0, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Heat Stress", description="Importance of avoiding heat stress. Higher values places more weight on areas with low heat stress.")
    shade_connectivity::P = Param(0.0, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Shade Connectivity", description="Higher values give preference to locations with high connectivity for shading deployments.")
    in_seed_connectivity::P = Param(1.0, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Incoming Connectivity (Seed)", description="Higher values give preference to locations with high incoming connectivity (i.e., receives larvae from other sites) for enhanced coral deployments.")
    out_seed_connectivity::P = Param(1.0, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Outgoing Connectivity (Seed)", description="Higher values give preference to locations with high outgoing connectivity (i.e., provides larvae to other sites) for enhanced coral deployments.")
    coral_cover_low::P = Param(0.0, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Low Coral Cover", description="Higher values give greater preference to sites with low coral cover for seeding deployments.")
    coral_cover_high::P = Param(0.0, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="High Coral Cover", description="Higher values give preference to sites with high coral cover for shading deployments.")
    seed_priority::P = Param(1.0, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Predecessor Priority (Seed)", description="Importance of seeding sites that provide larvae to priority reefs.")
    shade_priority::P = Param(0.0, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Predecessor Priority (Shade)", description="Importance of shading sites that provide larvae to priority reefs.")
    zone_seed::P = Param(0.0, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Zone Predecessor (Seed)", description="Importance of seeding sites that provide larvae to priority (target) zones.")
    zone_shade::P = Param(0.0, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Zone Predecessor (Shade)", description="Importance of shading sites that provide larvae to priority (target) zones.")
    coral_cover_tol::P = Param(0.2, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Low Area Tolerance", description="Tolerance for low proportional space for seeding deployments.")
    deployed_coral_risk_tol::P = Param(1.0, ptype="real", bounds=(0.75, 1.0), dists="unif",
        name="Risk Tolerance", description="Filters out sites with heat/wave stress above threshold.")
    use_dist::N = Param(1, ptype="integer", bounds=(0.0, 1.0 + 1.0), dists="unif",
        name="Use Distance Threshold", description="Turns distance sorting on or off.")
    dist_thresh::P = Param(0.1, ptype="real", bounds=(0.0, 1.0), dists="unif",
        name="Distance Threshold", description="Sites selected by MCDA must be further apart than median(dist)-dist_thresh*median(dist).")
    top_n::N = Param(10, ptype="integer", bounds=(5.0, 50.0 + 1.0), dists="unif",
        name="Top N", description="Replaces a given deployment site with a top-ranked site if it does not satisfy the minimum distance threshold.")
    depth_min::P = Param(5.0, ptype="real", bounds=(3.0, 5.0), dists="unif",
        name="Minimum Depth", description="Minimum depth for a site to be included for consideration.\nNote: This value will be replaced with the shallowest depth value found if all sites are found to be deeper than `depth_min + depth_offset`.")
    depth_offset::P = Param(10.0, ptype="real", bounds=(10.0, 25.0), dists="unif", name="Depth Offset",
        description="Offset from minimum depth, used to indicate maximum depth.")
end


struct EnvironmentalLayer{P} <: EcoModel
    dhw_scenario::P
    wave_scenario::P
end

function EnvironmentalLayer(dhw::AbstractArray, wave::AbstractArray)
    return EnvironmentalLayer(
        Param(1, bounds=(1.0, Float64(size(dhw, 3)) + 1.0), ptype="integer", dists="unif", name="DHW Scenario", description="DHW scenario member identifier."),
        Param(1, bounds=(1.0, Float64(size(wave, 3)) + 1.0), ptype="integer", dists="unif", name="Wave Scenario", description="Wave scenario member identifier.")
    )
end


"""
_coral_struct(field_defs::Dict)::Nothing

Helper function to dynamically create structs

https://stackoverflow.com/a/27084705/2694952
https://stackoverflow.com/questions/27083816/is-it-possible-to-create-types-in-julia-at-runtime


# Example

```julia
# Generate Coral struct with single attribute "a"
_coral_struct(Dict("a"=>200))

y = Coral()
# Coral{Int64}(200)

y.a
# 200
```
"""
function _coral_struct(field_defs::Dict)::Nothing
    s = IOBuffer()
    write(s, "Base.@kwdef struct Coral{P} <: EcoModel\n")

    for (f, v) in field_defs
        write(s, "$(f)::P = $(v)\n")
    end

    write(s, "end")
    eval(Meta.parse(String(take!(s))))

    return
end


"""
    create_coral_struct(bounds=(0.9, 1.1))

Generates Coral struct using the default parameter spec.

# Example
```julia
# Define coral struct with auto-generated parameter ranges
# (default in ADRIA is ± 10%, triangular distribution with peak at 0.5)
create_coral_struct()
coral = Coral()

# Recreate coral spec ± 50% from nominal values
create_coral_struct((0.5, 1.5))
coral = Coral()
```
"""
function create_coral_struct(bounds::Tuple{Float64,Float64}=(0.9, 1.1))::Nothing
    _, base_coral_params, x = coral_spec()

    coral_ids = x.coral_id

    struct_fields = Dict()
    for c_id in coral_ids
        for p in base_coral_params
            f_name = c_id * "_" * p
            f_val = x[x.coral_id.==c_id, p][1]
            struct_fields[f_name] = Param(f_val, ptype="real", bounds=(f_val * bounds[1], f_val * bounds[2], 0.5), dists="triang",
                name=human_readable_name(f_name, title_case=true), description="")
        end
    end

    _coral_struct(struct_fields)

    return
end


"""
    colony_areas()

Generate colony area data based on Bozec et al., [1].

# Returns
- `colony_area_mean_cm2` : mean colony areas in cm²
- `colony_diam_means_m` : mean colony diameter (in meters)


# References
1. Bozec, Y.-M., Rowell, D., Harrison, L., Gaskell, J., Hock, K.,
     Callaghan, D., Gorton, R., Kovacs, E. M., Lyons, M., Mumby, P.,
     & Roelfsema, C. (2021).
       Baseline mapping to support reef restoration and resilience-based
       management in the Whitsundays.
     https://doi.org/10.13140/RG.2.2.26976.20482
"""
function colony_areas()
    bin_edges = [0, 2, 5, 10, 20, 40, 80]

    # Diameters in cm
    mean_cm_diameters = bin_edges[1:end-1] + (bin_edges[2:end] - bin_edges[1:end-1]) / 2

    nclasses::Int64 = length(mean_cm_diameters)

    # The coral colony diameter bin edges (cm) are: 0, 2, 5, 10, 20, 40, 80
    # To convert to cover we locate bin means and calculate bin mean areas
    colony_diam_means = repeat(mean_cm_diameters', nclasses, 1)
    colony_area_mean_cm2 = @. pi * ((colony_diam_means / 2)^2)

    return colony_area_mean_cm2, (colony_diam_means ./ 100.0)
end


"""
    coral_spec()

Template for coral parameter values for ADRIA.
Includes "vital" bio/ecological parameters, to be filled with
sampled or user-specified values.

Any parameter added to the `params` DataFrame defined here will automatically be
made available to the ADRIA model.

Notes:
Values for the historical, temporal patterns of degree heating weeks
between bleaching years come from [1].

# Returns
- `params` : NamedTuple[taxa_names, param_names, params], taxa names, parameter
           names, and parameter values for each coral taxa, group and size class

# References
1. Lough, J. M., Anderson, K. D., & Hughes, T. P. (2018).
    Increasing thermal stress for tropical coral reefs: 1871-2017.
    Scientific Reports, 8(1), 6079.
    https://doi.org/10.1038/s41598-018-24530-9

2. Hall, V.R. & Hughes, T.P. 1996.
   Reproductive strategies of modular organisms:
     comparative studies of reef-building corals.
   Ecology, 77: 950 - 963.
   https://dx.doi.org/10.2307/2265514

3. Bozec, Y.-M., Rowell, D., Harrison, L., Gaskell, J., Hock, K.,
    Callaghan, D., Gorton, R., Kovacs, E. M., Lyons, M., Mumby, P.,
    & Roelfsema, C. (2021).
   Baseline mapping to support reef restoration and
     resilience-based management in the Whitsundays.
   https://doi.org/10.13140/RG.2.2.26976.20482

4. Bozec, Y.-M., Hock, K., Mason, R. A. B., Baird, M. E., Castro-Sanguino, C.,
    Condie, S. A., Puotinen, M., Thompson, A., & Mumby, P. J. (2022).
   Cumulative impacts across Australia's Great Barrier Reef: A mechanistic evaluation.
   Ecological Monographs, 92(1), e01494.
   https://doi.org/10.1002/ecm.1494
"""
function coral_spec()::NamedTuple
    # Below parameters pertaining to species are new. We now add size classes
    # to each coral species, but treat each coral size class as a 'species'.
    # Need a better word than 'species' to signal that we are using different
    # sizes within groups and real taxonomic species.

    params = DataFrame()

    # Coral species are divided into taxa and size classes
    taxa_names = String[
        "tabular_acropora_enhanced"
        "tabular_acropora_unenhanced"
        "corymbose_acropora_enhanced"
        "corymbose_acropora_unenhanced"
        "small_massives"
        "large_massives"
    ]

    # total number of "species" modelled in the current version.
    n_classes::Int64 = 6
    n_species::Int64 = length(taxa_names) * n_classes

    tn = repeat(taxa_names; inner=n_classes)

    # Create combinations of taxa names and size classes
    params.name = human_readable_name(tn; title_case=true)
    params.taxa_id = repeat(1:n_classes; inner=n_classes)

    params.class_id = repeat(1:n_classes, n_classes)
    params.coral_id = String[join(x, "_") for x in zip(tn, params.taxa_id, params.class_id)]

    # Ecological parameters
    # To be more consistent with parameters in ReefMod, IPMF and RRAP
    # interventions, we express coral abundance as colony numbers in different
    # size classes and growth rates as linear extension (in cm per year).

    colony_area_mean_cm², mean_colony_diameter_m = colony_areas()
    params.colony_area_cm2 = reshape(colony_area_mean_cm²', n_species)[:]

    ## Coral growth rates as linear extensions (Bozec et al 2021 S2, Table 1)
    # we assume similar growth rates for enhanced and unenhanced corals
    # all values in cm/year
    linear_extension = Array{Float64,2}([
        1.0 3.0 3.0 4.4 4.4 4.4     # Tabular Acropora Enhanced
        1.0 3.0 3.0 4.4 4.4 4.4     # Tabular Acropora Unenhanced
        1.0 3.0 3.0 3.0 3.0 3.0     # Corymbose Acropora Enhanced
        1.0 3.0 3.0 3.0 3.0 3.0     # Corymbose Acropora Unenhanced
        1.0 1.0 1.0 1.0 0.8 0.8     # small massives
        1.0 1.0 1.0 1.0 1.2 1.2])   # large massives

    # Convert linear extensions to delta coral in two steps.
    # First calculate what proportion of coral numbers that change size class
    # given linear extensions. This is based on the simple assumption that
    # coral sizes are evenly distributed within each bin

    bin_widths = Float64[2, 3, 5, 10, 20, 40]  # These bin widths have to line up with values in colony_areas()

    # Second, growth as transitions of cover to higher bins is estimated as
    # rate of growth per year
    params.growth_rate .= growth_rate(linear_extension, bin_widths)

    # Adjust growth rate for size class 6 to 20% of assumed value.
    params.growth_rate[params.class_id.==6] .= params.growth_rate[params.class_id.==6] .* 0.2

    # Scope for fecundity as a function of colony area (Hall and Hughes 1996)
    fec_par_a = Float64[1.03; 1.03; 1.69; 1.69; 0.86; 0.86]  # fecundity parameter a
    fec_par_b = Float64[1.28; 1.28; 1.05; 1.05; 1.21; 1.21]  # fecundity parameter b
    min_size_full_fec_cm2 = Float64[123.0; 123.0; 134.0; 134.0; 38.0; 38.0]

    # fecundity as a function of colony basal area (cm2) from Hall and Hughes 1996
    # unit is number of larvae per colony
    colony_area_cm2 = pi .* ((mean_colony_diameter_m .* 100.0) ./ 2.0) .^ 2
    fec = exp.(log.(fec_par_a) .+ fec_par_b .* log.(colony_area_cm2)) ./ 0.1
    fec[colony_area_cm2.<min_size_full_fec_cm2] .= 0.0

    # Smallest size class do not reproduce
    fec[:, 1:2] .= 0.0

    # then convert to number of larvae produced per m2
    fec_m² = fec ./ (pi .* (mean_colony_diameter_m ./ 2.0) .^ 2)  # convert from per colony area to per m2
    params.fecundity = fec_m²'[:]

    ## Mortality
    # Wave mortality risk : wave damage for the 90 percentile of routine wave stress
    wavemort90 = Array{Float64,2}([
        0.0 0.0 0.0 0.0 0.05 0.1  # Tabular Acropora Enhanced
        0.0 0.0 0.0 0.0 0.05 0.1  # Tabular Acropora Unenhanced
        0.0 0.0 0.0 0.0 0.02 0.05  # Corymbose Acropora Enhanced
        0.0 0.0 0.0 0.0 0.02 0.05  # Corymbose Acropora Unenhanced
        0.0 0.0 0.0 0.0 0.0 0.0  # small massives and encrusting
        0.0 0.0 0.0 0.0 0.0 0.0])  # large massives

    params.wavemort90 = wavemort90'[:]

    # Background mortality taken from Bozec et al. 2022 (Supplementary 2, Table S1)
    # Using values for:
    # - juvenile mortality (first two columns)
    # - < 5cm² (Columns 1 and 2)
    # - < 250cm² (Columns 3 and 4)
    # - > 250cm² (Columns 5 and 6)
    # Values for size class 4 are then interpolated by K.A
    mb = Array{Float64,2}([
        0.2 0.2 0.19 0.125 0.05 0.0    # Tabular Acropora Enhanced
        0.2 0.2 0.19 0.125 0.05 0.0    # Tabular Acropora Unenhanced
        0.2 0.2 0.172 0.113 0.06 0.04    # Corymbose Acropora Enhanced
        0.2 0.2 0.172 0.113 0.06 0.04    # Corymbose Acropora Unenhanced
        0.2 0.2 0.04 0.026 0.02 0.02    # Small massives and encrusting
        0.2 0.2 0.04 0.026 0.02 0.02])   # Large massives
    params.mb_rate = mb'[:]

    # Bleaching sensitivity of each coral group
    # Bozec et al., (2022)
    bleaching_sensitivity = Float64[
        1.50 1.50 1.50 1.50 1.50 1.50  # Tabular Acropora Enhanced (Arborescent staghorn corals)
        1.50 1.50 1.50 1.50 1.50 1.50  # Tabular Acropora Unenhanced
        1.40 1.40 1.40 1.40 1.40 1.40  # Corymbose Acropora Enhanced
        1.40 1.40 1.40 1.40 1.40 1.40  # Corymbose Acropora Unenhanced
        0.25 0.25 0.25 0.25 0.25 0.25  # Small massives and encrusting
        0.25 0.25 0.25 0.25 0.25 0.25] # Large massives
    params.bleaching_sensitivity = bleaching_sensitivity'[:]

    # Get perturbable coral parameters
    # i.e., the parameter names not defined in the second list
    param_names = setdiff(names(params), ["name", "taxa_id", "class_id", "coral_id"])

    return (taxa_names=taxa_names, param_names=param_names, params=params)
end
