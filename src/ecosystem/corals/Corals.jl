using DataFrames
import ModelParameters: Model

# Upper bound offset to use when re-creating critical DHW distributions
const HEAT_UB = 10.0

"""
    functional_group_names()::Vector{Symbol}

Name of functional groups represented by ADRIAmod.
"""
function functional_group_names()::Vector{Symbol}
    return [
        #:arborescent_Acropora,
        :tabular_Acropora,
        :corymbose_Acropora,
        :corymbose_non_Acropora,  # and Pocillopora
        :small_massives,
        :large_massives
    ]
end

"""
    colony_mean_area(colony_diam_means::Array{T})::Array{T} where {T<:Real}

Generate mean colony areas for given colony diameter(s).

# Arguments
- `colony_diam_means` : mean colony diameter (in meters)
"""
function colony_mean_area(colony_diam_means::Array{T})::Array{T} where {T<:Float64}
    return pi .* ((colony_diam_means ./ 2.0) .^ 2)
end

"""
    linear_extensions()

Linear extensions. The values are converted from `cm` to the desired unit.
The default unit is `m`.
"""
function linear_extensions(; unit=:m)::Matrix{Float64}
    return [
        0.609456 1.07184 2.55149 5.07988 9.45091 16.8505 0.0;
        0.768556 1.22085 1.86447 2.82297 3.52938 3.00422 0.0;
        0.190455 0.343747 0.615467 0.97477 1.70079 2.91729 0.0;
        0.318034 0.47385 0.683729 0.710587 0.581085 0.581085 0.0;
        0.122478 0.217702 0.382098 0.718781 1.24172 2.08546 0.0
    ] .* linear_scale(:cm, unit)
end

"""
    bin_edges()

Helper function defining coral colony diameter bin edges. The values are converted from `cm`
to the desired unit. The default unit is `m`.
"""
function bin_edges(; unit=:m)
    return Matrix(
        [
            2.5 5.0 7.5 10.0 20.0 40.0 100.0 150.0;
            2.5 5.0 7.5 10.0 20.0 35.0 50.0 100.0;
            2.5 5.0 7.5 10.0 15.0 20.0 40.0 50.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0
        ]
    ) .* linear_scale(:cm, unit)
end
# function bin_edges()
#     return Matrix([
#         0.0 1.0 2.0 6.0 15.0 36.0 89.0 90.0;
#         0.0 1.0 2.0 4.0  9.0 18.0 38.0 39.0;
#         0.0 1.0 2.0 4.0  7.0 14.0 27.0 28.0;
#         0.0 1.0 2.0 5.0  8.0 12.0 26.0 27.0;
#         0.0 1.0 2.0 4.0  9.0 19.0 40.0 41.0
#     ])
# end

"""
    bin_widths()

Helper function defining coral colony diameter bin widths.
"""
function bin_widths()
    return bin_edges()[:, 2:end] .- bin_edges()[:, 1:(end - 1)]
end

"""
    planar_area_params()

Colony planar area parameters (see Fig 2B in Aston et al., [1])
First column is `b`, second column is `a`
log(S) = b + a * log(x)
"""
function planar_area_params()
    return Array{Float64,2}([
        # -8.97 3.14   # Abhorescent Acropora (using branching porites parameters as similar method of growing ever expanding colonies).
        -8.95 2.80   # Tabular Acropora
        -9.13 2.94   # Corymbose Acropora
        -8.90 2.94   # Corymbose non-Acropora (using branching pocillopora values from fig2B)
        -8.87 2.30   # Small massives
        -8.87 2.30   # Large massives
    ])
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
    # The coral colony diameter bin edges (cm) are: 0, 2, 5, 10, 20, 40, 80
    edges = bin_edges(; unit=:cm)

    # Diameters in cm
    mean_cm_diameters =
        edges[:, 1:(end - 1)] + (edges[:, 2:end] - edges[:, 1:(end - 1)]) / 2.0

    # To convert to cover we locate bin means and calculate bin mean areas
    colony_area_mean_cm2 = colony_mean_area(mean_cm_diameters)

    return colony_area_mean_cm2, linear_scale.(mean_cm_diameters, :cm, :m)
end

function bins_bounds(mean_diam::Matrix{Float64})::Matrix{Float64}
    bins::Matrix{Float64} = zeros(size(mean_diam)...)
    bins[:, 1] .= mean_diam[:, 1] .* 2
    for i in 2:(size(mean_diam)[2])
        bins[:, i] .= (mean_diam[:, i] .- bins[:, i - 1]) .* 2 .+ bins[:, i - 1]
    end
    return bins
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

5. Bairos-Novak, K.R., Hoogenboom, M.O., van Oppen, M.J.H., Connolly, S.R., 2021.
   Coral adaptation to climate change: Meta-analysis reveals high heritability across
     multiple traits.
   Global Change Biology 27, 5694-5710.
   https://doi.org/10.1111/gcb.15829
"""
function coral_spec()::NamedTuple
    # Below parameters pertaining to species are new. We now add size classes
    # to each coral species, but treat each coral size class as a 'species'.
    # Need a better word than 'species' to signal that we are using different
    # sizes within groups and real taxonomic species.

    params = DataFrame()

    # Coral species are divided into taxa and size classes
    group_names = string.(functional_group_names())

    # Coral growth rates as linear extensions.
    # All values in cm/year and are from (unpublished) ecoRRAP data.
    _linear_extensions::Matrix{Float64} = linear_extensions()

    # number of functional groups and size classes modelled in the current version.
    n_groups::Int64, n_sizes::Int64 = size(_linear_extensions)
    n_groups_and_sizes::Int64 = n_groups * n_sizes

    tn = repeat(group_names; inner=n_sizes)

    # Create combinations of taxa names and size classes
    params.name = human_readable_name(tn; title_case=true)
    params.taxa_id = repeat(1:n_groups; inner=n_sizes)

    params.class_id = repeat(1:n_sizes, n_groups)::Vector{Int64}
    params.coral_id = String[join(x, "_") for x in zip(tn, params.taxa_id, params.class_id)]

    # Ecological parameters
    # To be more consistent with parameters in ReefMod, C~Scape and RRAP
    # interventions, we express coral abundance as colony numbers in different
    # size classes and growth rates as linear extension (in cm per year).
    colony_area_mean_cm², mean_colony_diameter_m = colony_areas()
    params.mean_colony_diameter_m = reshape(mean_colony_diameter_m', n_groups_and_sizes)[:]
    params.linear_extension = reshape(_linear_extensions', n_groups_and_sizes)[:]

    # Convert linear extensions to delta coral in two steps.
    # First calculate what proportion of coral numbers that change size class
    #     given linear extensions. This is based on the simple assumption that
    #     coral sizes are evenly distributed within each bin.
    # Second, growth as transitions of cover to higher bins is estimated as
    #     rate of growth per year.
    params.growth_rate .= reshape(
        growth_rate(_linear_extensions, bin_widths()), n_groups_and_sizes
    )[:]

    # Scope for fecundity as a function of colony area (Hall and Hughes 1996)
    # Corymbose non-acropora uses the Stylophora data from Hall and Hughes with interpolation
    fec_par_a = Float64[1.03; 1.69; 0.02; 0.86; 0.86]  # fecundity parameter a
    fec_par_b = Float64[1.28; 1.05; 2.27; 1.21; 1.21]  # fecundity parameter b
    min_colony_area_full_fec = Float64[123.0; 134.0; 134.0; 38.0; 38.0]

    # fecundity as a function of colony basal area (cm2) from Hall and Hughes 1996
    # unit is number of larvae per colony
    colony_area_cm2 = colony_mean_area(mean_colony_diameter_m .* 100.0)
    fec = exp.(log.(fec_par_a) .+ fec_par_b .* log.(colony_area_cm2)) ./ 0.1

    # Colonies with area (in cm2) below indicated size are not fecund (reproductive)
    fec[colony_area_cm2 .< min_colony_area_full_fec] .= 0.0

    # then convert to number of larvae produced per m2
    fec_m² = fec ./ (colony_mean_area(mean_colony_diameter_m)) # convert from per colony area to per m2
    params.fecundity = fec_m²'[:]

    # Mortality
    # Survival rates taken from (unpublished) ecoRRAP data.
    # Background mortality is then the complement.
    # survival_rate::Matrix{Float64} = [
    #     0.859017851 0.858528906 0.857044217 0.856477498 0.856104353 0.855852241 0.855852241;    # Tabular Acropora
    #     0.865006527 0.87915437 0.892044073 0.905304164 0.915373252 0.925707536 0.925707536;     # Corymbose Acropora
    #     0.953069031 0.959152694 0.964460394 0.968306361 0.972598906 0.97621179 0.97621179;     # Corymbose non-Acropora
    #     0.869976692 0.938029324 0.977889252 0.987199004 0.99207702 0.996931548 0.996931548;     # Small massives and encrusting
    #     0.9782479 0.979496637 0.980850254 0.982178103 0.983568572 0.984667677 0.984667677       # Large massives
    # ]
    survival_rate::Matrix{Float64} = [
        0.6 0.76 0.805 0.76 0.85 0.86 0.86;    # Tabular Acropora
        0.6 0.76 0.77 0.875 0.83 0.90 0.90;    # Corymbose Acropora
        0.52 0.77 0.77 0.875 0.89 0.97621179 0.97621179;                # Corymbose non-Acropora
        0.72 0.87 0.77 0.98 0.996931548 0.996931548 0.996931548;        # Small massives and encrusting
        0.58 0.87 0.78 0.983568572 0.984667677 0.984667677 0.984667677  # Large massives
    ]
    mb = 1.0 .- survival_rate
    params.mb_rate = mb'[:]

    upper_bound::Matrix{Float64} = bin_edges()[:, 2:end]
    params.bin_ub = reshape(upper_bound', n_groups_and_sizes)[:]

    # Natural adaptation / heritability
    # Values here informed by Bairos-Novak et al., (2022) and (unpublished) data from
    # Hughes et al., (2018)
    # Mean and std for each species (row) and size class (cols)
    params.dist_mean = repeat(
        Float64[
            # 3.345484656,  # arborescent Acropora
            3.751612251,  # tabular Acropora
            4.081622683,  # corymbose Acropora
            4.487465256,  # Pocillopora + non-Acropora corymbose
            6.165751937,  # Small massives and encrusting
            7.153507902   # Large massives
        ]; inner=n_sizes)

    params.dist_std = repeat(
        Float64[
            # 2.590016677,  # arborescent Acropora
            2.904433676,  # tabular Acropora
            3.159922076,  # corymbose Acropora
            3.474118416,  # Pocillopora + non-Acropora corymbose
            4.773419097,  # Small massives and encrusting
            5.538122776   # Large massives
        ]; inner=n_sizes)

    # Get perturbable coral parameters
    # i.e., the parameter names not defined in the second list
    param_names = setdiff(names(params), ["name", "taxa_id", "class_id", "coral_id"])

    return (taxa_names=group_names, param_names=param_names, params=params)
end

"""
    _coral_struct(field_defs::Dict)::Nothing

Helper function to dynamically create coral struct.

https://stackoverflow.com/a/27084705/2694952
https://stackoverflow.com/questions/27083816/is-it-possible-to-create-types-in-julia-at-runtime


# Example

```julia
# Generate Coral struct with single attribute "a"
_coral_struct(OrderedDict("a"=>200))

y = Coral()
# Coral{Int64}(200)

y.a
# 200
```
"""
function _coral_struct(field_defs::OrderedDict)::Nothing
    s = IOBuffer()
    write(s, "Base.@kwdef struct Coral{P,P2} <: EcoModel\n")

    for (f, v) in field_defs
        if f == "heritability"
            write(s, "$(f)::P2 = $(v)\n")
            continue
        end

        write(s, "$(f)::P = $(v)\n")
    end

    write(s, "end")
    eval(Meta.parse(String(take!(s))))

    return nothing
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
    _, base_coral_params, p_vals = coral_spec()

    struct_fields = OrderedDict{String,Param}()
    struct_fields["heritability"] = Factor(
        0.3;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.25, 0.5),
        name="Heritability",
        description="Heritability of DHW tolerance."
    )

    for c_id in p_vals.coral_id
        for p in base_coral_params
            f_name::String = c_id * "_" * p
            f_val = p_vals[p_vals.coral_id .== c_id, p][1]
            struct_fields[f_name] = Factor(
                f_val;
                ptype="continuous",
                dist=TriangularDist,
                dist_params=(f_val * bounds[1], f_val * bounds[2], f_val),
                name=human_readable_name(f_name; title_case=true),
                description=""
            )
        end
    end

    _coral_struct(struct_fields)

    return nothing
end

# Generate base coral struct from default spec.
# Have to call this before including specification methods
create_coral_struct()

"""
    to_coral_spec(m::Model)::DataFrame

Convert Coral Model specification to a coral spec DataFrame
"""
function to_coral_spec(m::Coral)::DataFrame
    _, pnames, spec = coral_spec()
    val_df = DataFrame(Model(m))

    return _update_coral_spec(spec, pnames, val_df)
end

"""
    to_coral_spec(coral_df::DataFrame)::DataFrame

Convert dataframe of model parameters to a coral spec.
"""
function to_coral_spec(coral_df::DataFrame)::DataFrame
    _, pnames, spec = coral_spec()

    return _update_coral_spec(spec, pnames, coral_df)
end

function _update_coral_spec(
    spec::DataFrame, pnames::Vector{String}, coral_params::DataFrame
)::DataFrame
    fnames::Vector{String} = String.(coral_params[!, :fieldname])
    for p in pnames
        target::BitVector = occursin.(p, fnames)
        for tn in fnames[target]
            idx::BitVector = spec.coral_id .== rsplit(tn, "_$p"; keepempty=false)
            spec[idx, [p]] .= coral_params[coral_params.fieldname .== tn, :val]
        end
    end

    return spec
end

"""
    _update_coral_factors(spec::DataFrame, coral_params::DataFrame)::DataFrame

Update scenario specification with updated coral parameters/factors.
Allows changes to `coral_spec()` to be represented without requiring a session restart.
For debugging purposes only.

# Examples
```julia
using Revise

dom = ADRIA.load_domain("path to a Domain", "45")
default_scen = ADRIA.param_table(dom)

# ... modify the coral specification at some point
# this would normally require restarting the Julia session for changes to take effect.
# Instead, we update the specification with the updated values instead.
default_scen = ADRIA._update_coral_factors(default_scen, ADRIA.coral_spec().params)
```
"""
function _update_coral_factors(spec::DataFrame, coral_params::DataFrame)::DataFrame
    c_ids = coral_params.coral_id
    for factor in [
        "mean_colony_diameter_m",
        "growth_rate",
        "fecundity",
        "mb_rate",
        "dist_mean",
        "dist_std"
    ]
        for id in c_ids
            spec[!, "$(id)_$(factor)"] = coral_params[coral_params.coral_id .== id, factor]
        end
    end

    return spec
end

function to_coral_spec(inputs::YAXArray)::DataFrame
    _, pnames, spec = coral_spec()

    coral_ids::Vector{String} = spec[:, :coral_id]
    for p in pnames
        # wrapping `p` in an array is necessary so update of DF works
        spec[!, [p]] .= Array(inputs[At(coral_ids .* "_" .* p)])
    end

    return spec
end
function to_coral_spec(inputs::DataFrameRow)::DataFrame
    ins = DataCube(Vector(inputs); factors=names(inputs))
    return to_coral_spec(ins)
end

"""
    n_sizes(n_groups::Int64, n_group_sizes::Int64)::Int64

Number of size classes. `n_groups` must exactly divide `n_group_sizes`.

# Arguments
- `n_groups` : Number of functional groups.
- `n_group_sizes` : Number of functional groups multiplied by number of size classes.
"""
function n_sizes(n_groups::Int64, n_group_sizes::Int64)::Int64
    if n_group_sizes % n_groups != 0
        throw(
            ArgumentError(
                "Number of groups must divide n_group_sizes. " *
                "n_group_sizes: $(n_group_sizes), n_groups: $(n_groups)"
            )
        )
    end
    return Int64(n_group_sizes / n_groups)
end

function group_indices(n_sizes::Int64, n_group_sizes::Int64)::Vector{UnitRange{Int64}}
    return [i:(i + (n_sizes - 1)) for i in 1:n_sizes:n_group_sizes]
end
