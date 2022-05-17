using Setfield

using DataFrames
using ModelParameters
import ModelParameters: update!, Model


abstract type EcoModel end


"""Set a model parameter value directly."""
function set!(p::Param, val)
    if hasproperty(p, :ptype) && p.ptype == "integer" && !isinteger(val)
        # For integer/categorical parameters, take floor of v + 1, capping to the upper bound
        val = floor(val + 1)
        val = min(val, p.bounds[2])
    end

    # @set! p.val = val
    return val
end


"""Update a given model with new uncertain parameter values."""
function update!(m::Model, vals::Union{Vector,Tuple,Array})::Nothing
    m[:val] = map((x) -> set!(x...), zip(params(m), vals))

    return
end


Base.@kwdef struct Intervention{N,P} <: EcoModel
    # Intervention Parameters
    guided::N = Param(0, ptype="integer", bounds=(0, 5)) # Guided, choice of MCDA approach
    seed_TA::N = Param(0, ptype="integer", bounds=(0, 500000)) # Seed1, integer, number of Enhanced TA to seed
    seed_CA::N = Param(0, ptype="integer", bounds=(0, 500000)) # Seed2, integer, number of Enhanced CA to seed
    fogging::P = Param(0.2, ptype="real", bounds=(0.0, 0.3)) # fogging, float, assumed percent reduction in bleaching mortality
    SRM::P = Param(0.0, ptype="real", bounds=(0.0, 12.0)) # SRM, float, reduction in DHWs due to shading
    a_adapt::P = Param(0.0, ptype="real", bounds=(0.0, 12.0)) # Aadpt, float, float, increased adaptation rate
    n_adapt::P = Param(0.0, ptype="real", bounds=(0.0, 0.05)) # Natad, float, natural adaptation rate
    seed_years::N = Param(10, ptype="integer", bounds=(5, 16)) # Seedyrs, integer, years into simulation during which seeding is considered
    shade_years::N = Param(10, ptype="integer", bounds=(5, 74)) # Shadeyrs, integer, years into simulation during which shading is considered
    seed_freq::N = Param(5, ptype="integer", bounds=(0, 6)) # Seedfreq, integer, yearly intervals to adjust seeding site selection (0 is set and forget)
    shade_freq::N = Param(1, ptype="integer", bounds=(0, 6)) # Shadefreq, integer, yearly intervals to adjust shading (fogging) site selection (0 is set and forget)
    seed_year_start::N = Param(2, ptype="integer", bounds=(2, 26)) # Seedyr_start, integer, seed intervention start offset from simulation start
    shade_year_start::N = Param(2, ptype="integer", bounds=(2, 26)) # Shadeyr_start, integer, shade intervention start offset from simulation start
end


Base.@kwdef struct Criteria{P} <: EcoModel
    wave_stress::P = Param(1.0, ptype="real", bounds=(0.0, 1.0))
    heat_stress::P = Param(1.0, ptype="real", bounds=(0.0, 1.0))
    shade_connectivity::P = Param(0.0, ptype="real", bounds=(0.0, 1.0))
    seed_connectivity::P = Param(1.0, ptype="real", bounds=(0.0, 1.0))
    coral_cover_high::P = Param(0.0, ptype="real", bounds=(0.0, 1.0))
    coral_cover_low::P = Param(1.0, ptype="real", bounds=(0.0, 1.0))
    seed_priority::P = Param(1.0, ptype="real", bounds=(0.0, 1.0))
    shade_priority::P = Param(0.0, ptype="real", bounds=(0.0, 1.0))
    deployed_coral_risk_tol::P = Param(5.0, ptype="real", bounds=(0.0, 1.0))
    depth_min::P = Param(5.0, ptype="real", bounds=(3.0, 5.0))     # minimum depth
    depth_offset::P = Param(5.0, ptype="real", bounds=(5.0, 6.0))  # offset from minimum depth to indicate maximum depth**
end
# **This is simply to avoid parameterization/implementation
#   that requires one parameter to be greater than another.


"""
_coral_struct(field_defs::Dict)::Nothing

Helper function to dynamically create structs

https://stackoverflow.com/a/27084705/2694952
https://stackoverflow.com/questions/27083816/is-it-possible-to-create-types-in-julia-at-runtime


Example
-------
julia> _coral_struct(Dict("a"=>200))

julia> y = Coral()
Coral{Int64}(200)

julia> y.a
    200
"""
function _coral_struct(field_defs::Dict)::Nothing
    s::String = "Base.@kwdef struct Coral{P} <: EcoModel\n"

    for (f, v) in field_defs
        s = s * "$(f)::P = $(v)\n"
    end

    eval(Meta.parse(s * "end"))

    return
end


"""
    create_coral_struct(bounds=(0.9, 1.1))

Generates Coral struct

Example
-------
    # Define coral struct with auto-generated parameter ranges (default in ADRIA is ± 10%)
    create_coral_struct()
    coral = Coral()

    # Recreate coral spec ± 50% from nominal values
    create_coral_struct((0.5, 1.5))
    coral = Coral()
"""
function create_coral_struct(bounds=(0.9, 1.1))
    _, base_coral_params, x = coral_spec()

    coral_ids = x.coral_id

    struct_fields = Dict()
    for c_id in coral_ids
        for p in base_coral_params
            f_name = c_id * "_" * p
            f_val = x[x.coral_id.==c_id, p][1]
            struct_fields[f_name] = Param(f_val, ptype="real", bounds=(f_val * bounds[1], f_val * bounds[2]))
        end
    end

    _coral_struct(struct_fields)

    return
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

Returns
-------
params : NamedTuple[taxa_names, param_names, params], taxa names, parameter 
         names, and parameter values for each coral taxa, group and size class

References
----------
1. Lough, J.M., Anderson, K.D. and Hughes, T.P. (2018)
       'Increasing thermal stress for tropical coral reefs: 1871-2017',
       Scientific Reports, 8(1), p. 6079.
       doi: 10.1038/s41598-018-24530-9.

2. Bozec, Y.-M., Rowell, D., Harrison, L., Gaskell, J., Hock, K.,
     Callaghan, D., Gorton, R., Kovacs, E. M., Lyons, M., Mumby, P.,
     & Roelfsema, C. (2021).
       Baseline mapping to support reef restoration and resilience-based
       management in the Whitsundays.
     https://doi.org/10.13140/RG.2.2.26976.20482

3. Hall,V.R. & Hughes, T.P. 1996. Reproductive strategies of modular
      organisms: comparative studies of reef-building corals. Ecology,
      77: 950 - 963.
"""
function coral_spec()::NamedTuple
    # Below parameters pertaining to species are new. We now add size classes
    # to each coral species, but treat each coral size class as a 'species'.
    # Need a better word than 'species' to signal that we are using different
    # sizes within groups and real taxonomic species.

    params = DataFrame();

    param_names = ["growth_rate", "fecundity", "wavemort90", "mb_rate", "bleach_resist"]

    # Coral species are divided into taxa and size classes
    taxa_names = [
        "tabular_acropora_enhanced";
        "tabular_acropora_unenhanced";
        "corymbose_acropora_enhanced";
        "corymbose_acropora_unenhanced";
        "small_massives";
        "large_massives"
        ];

    tn = repeat(taxa_names, 6, 1)

    size_cm = [2; 5; 10; 20; 40; 80]
    size_class_means_from = [1; 3.5; 7.5; 15; 30; 60]
    size_class_means_to = [size_class_means_from[2:end]; 100.0];

    # total number of "species" modelled in the current version.
    nclasses = length(size_cm);
    nspecies = length(taxa_names) * nclasses;

    # Create combinations of taxa names and size classes
    params.name = human_readable_name(tn, true);
    params.taxa_id = repeat(1:nclasses, inner=nclasses);

    params.class_id = repeat(1:nclasses, nclasses);
    params.size_cm = repeat(size_cm, nclasses);

    params.coral_id = ["$(x[1])_$(x[2])_$(x[3])" for x in zip(tn, params.taxa_id, params.class_id)]

    # rec = [0.00, 0.01, 0.00, 0.01, 0.01, 0.01];
    # params.recruitment_factor = repmat(rec, 1, nclasses)';

    ## Ecological parameters

    # To be more consistent with parameters in ReefMod, IPMF and RRAP
    # interventions, we express coral abundance as colony numbers in different
    # size classes and growth rates as linear extention (in cm per year).

    ### Base covers
    #First express as number of colonies per size class per 100m2 of reef

    # NOTE: These values are currently being overwritten by init_coral_cover
    # in the ADRIA class (see init_coral_cover method in ADRIA.m)
    base_coral_numbers =
        [0 0 0 0 0 0;           # Tabular Acropora Enhanced
         0 0 0 0 0 0;           # Tabular Acropora Unenhanced
         0 0 0 0 0 0;           # Corymbose Acropora Enhanced
         200 100 100 50 30 10;  # Corymbose Acropora Unenhanced
         200 100 200 30 0 0;    # small massives
         0 0 0 0 0 0];             # large massives

    # To convert to covers we need to first calculate the area of colonies,
    # multiply by how many corals in each bin, and divide by reef area

    # The coral colony diameter bin edges (cm) are: 0, 2, 5, 10, 20, 40, 80
    # To convert to cover we locate bin means and calculate bin mean areas
    colony_diam_means_from = repeat(size_class_means_from', nclasses, 1)
    colony_diam_means_to = repeat(size_class_means_to', nclasses', 1)

    colony_area_m2_from = @. pi * ((colony_diam_means_from / 2)^2) / (10^4)
    colony_area = @. pi * ((colony_diam_means_to / 2)^2)
    params.colony_area_cm2 = reshape(colony_area', nspecies)
    colony_area_m2_to = colony_area ./ (10^4)

    a_arena = 100 # m2 of reef arena where corals grow, survive and reproduce

    # convert to coral covers (proportions) and convert to vector
    params.basecov = @. base_coral_numbers'[:] * colony_area_m2_from'[:] / a_arena

    ## Coral growth rates as linear extensions (Bozec et al 2021 Table S2)
    # we assume similar growth rates for enhanced and unenhanced corals
    linear_extension =
       [1 3 3 4.4 4.4 4.4;  # Tabular Acropora Enhanced
        1 3 3 4.4 4.4 4.4;   # Tabular Acropora Unenhanced
        1 3 3 3 3 3;         # Corymbose Acropora Enhanced
        1 3 3 3 3 3;         # Corymbose Acropora Unenhanced
        1 1 1 1 0.8 0.8;     # small massives
        1 1 1 1 1.2 1.2];       # large massives

    # Convert linear extensions to delta coral in two steps.
    # First calculate what proportion of coral numbers that change size class
    # given linear extensions. This is based on the simple assumption that
    # coral sizes are evenly distributed within each bin
    bin_widths = [2, 3, 5, 10, 20, 40];
    diam_bin_widths = repeat(bin_widths, nclasses, 1)
    prop_change = @views linear_extension'[:] ./ diam_bin_widths

    # Second, growth as transitions of cover to higher bins is estimated as
    params.growth_rate = vec(prop_change .* (colony_area_m2_to'[:] ./ colony_area_m2_from'[:]))

    # note that we use proportion of bin widths and linear extension to estimate
    # number of corals changing size class, but we use the bin means to estimate
    # the cover equivalent because we assume coral sizes shift from edges to mean
    # over the year (used in 'growthODE4()').

    # Scope for fecundity as a function of colony area (Hall and Hughes 1996)
    fec_par_a = [1.02; 1.02; 1.69; 1.69; 0.86; 0.86]; # fecundity parameter a
    fec_par_b = [1.28; 1.28; 1.05; 1.05; 1.21; 1.21]; # fecundity parameter b

    # fecundity as a function of colony basal area (cm2) from Hall and Hughes 1996
    # unit is number of larvae per colony
    fec = exp.(fec_par_a .+ fec_par_b .* log.(colony_area_m2_from * 10^4))

    # then convert to number of larvae produced per m2
    fec_m2 = fec ./ colony_area_m2_from;  # convert from per colony area to per m2
    params.fecundity = fec_m2'[:];

    ## Mortality
    # Wave mortality risk : wave damage for the 90 percentile of routine wave stress
    wavemort90 =
       [0 0 0.00 0.00 0.00 0.00;   # Tabular Acropora Enhanced
        0 0 0.00 0.00 0.00 0.00;   # Tabular Acropora Unenhanced
        0 0 0.00 0.00 0.00 0.00;   # Corymbose Acropora Enhanced
        0 0 0.00 0.00 0.00 0.00;   # Corymbose Acropora Unenhanced
        0 0 0.00 0.00 0.00 0.00;   # Small massives
        0 0 0.00 0.00 0.00 0.00];  # Large massives

    params.wavemort90 = wavemort90'[:];

    # Background mortality taken from Bozec et al. 2021 (Table S2)
    mb = [0.20 0.19 0.15 0.098 0.098 0.098;    # Tabular Acropora Enhanced
          0.20 0.19 0.15 0.098 0.098 0.098;    # Tabular Acropora Unenhanced
          0.20 0.17 0.12 0.088 0.088 0.088;    # Corymbose Acropora Enhanced
          0.20 0.17 0.12 0.088 0.088 0.088;    # Corymbose Acropora Unenhanced
          0.20 0.10 0.04 0.030 0.020 0.020;    # Small massives and encrusting
          0.20 0.10 0.04 0.030 0.020 0.020];   # Large massives

    params.mb_rate = mb'[:];

    # Background rates of natural adaptation. User-defined natad rates will be
    # added to these
    # natad = zeros(36);
    # params.n_adapt = natad;

    # Estimated bleaching resistance (as DHW) relative to the assemblage
    # response for 2016 bleaching on the GBR (based on Hughes et al. 2018).
    bleach_resist = [
        0.0 0.0 0.0 0.0 0.0 0.0;  # Tabular Acropora Enhanced
        0.0 0.0 0.0 0.0 0.0 0.0;  # Tabular Acropora Unenhanced
        0.0 0.0 0.0 0.0 0.0 0.0;  # Corymbose Acropora Enhanced
        0.0 0.0 0.0 0.0 0.0 0.0;  # Corymbose Acropora Unenhanced
        1.5 1.5 1.5 1.5 1.5 1.5;  # Small massives and encrusting
        1.0 1.0 1.0 1.0 1.0 1.0]; # Large massives 

    params.bleach_resist = bleach_resist'[:];

    return (taxa_names=taxa_names, param_names=param_names, params=params)
end
