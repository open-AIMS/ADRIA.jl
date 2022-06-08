module metrics

using Interpolations, Statistics

using DataFrames
import ADRIA: coral_spec, ResultSet


function relative_cover(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    return dropdims(sum(X, dims=2), dims=2)  # sum over all species and size classes
end

"""
    coral_cover(X)::NamedTuple

Converts outputs from scenario runs to relative cover of the four different coral taxa.

# Returns
NamedTuple
    - total_cover : total coral cover
    - enhanced_tab_acr : cover of enhanced tabular acropora
    - unenhanced_tab_acr : area covered by unenhanced tabular acropora
    - enhanced_cor_acr : area covered by enhanced corymbose acropora
    - unenhanced_cor_acr : area covered by unenhanced corymbose acropora
    - tab_acr : cover of tabular acropora
    - cor_acr : cover of corymbose acropora
    - small_enc : cover of small encrusting
    - large_mass : cover of large massives
    - juveniles : area covered by juveniles
    - large : area covered by large mature corals
"""
function coral_cover(X::AbstractArray{<:Real})::NamedTuple
    # Relative total coral cover
    TC::AbstractArray{<:Real} = relative_cover(X)  # sum over all species and size classes

    _, _, cs_p = coral_spec()

    screen = (x, idx) -> findall(x .== idx)

    sc1 = X[:, screen(cs_p.taxa_id, 1), :, :, :]
    sc2 = X[:, screen(cs_p.taxa_id, 2), :, :, :]
    sc3 = X[:, screen(cs_p.taxa_id, 3), :, :, :]
    sc4 = X[:, screen(cs_p.taxa_id, 4), :, :, :]

    C1 = sc1 .+ sc2  # enhanced to unenhanced tabular Acropora
    C2 = sc3 .+ sc4  # enhanced to unenhanced corymbose Acropora
    C3 = X[:, screen(cs_p.taxa_id, 5), :, :, :]  # Encrusting and small massives
    C4 = X[:, screen(cs_p.taxa_id, 6), :, :, :]  # Large massives

    # Cover of juvenile corals (< 5cm diameter)
    juv_groups = X[:, screen(cs_p.class_id, 1), :, :, :] .+ X[:, screen(cs_p.class_id, 2), :, :, :]
    juv_all = dropdims(sum(juv_groups, dims=2), dims=2)

    large_corals = X[:, screen(cs_p.class_id, 5), :, :, :] + X[:, screen(cs_p.class_id, 6), :, :, :]
    large_all = dropdims(sum(large_corals, dims=2), dims=2)

    covers = (total_cover = TC,
              enhanced_tab_acr = sc1,
              unenhanced_tab_acr = sc2,
              enhanced_cor_acr = sc3,
              unenhanced_cor_acr = sc4,
              tab_acr = C1, cor_acr = C2, 
              small_enc = C3, large_mass = C4,
              juveniles=juv_all, large=large_all)

    return covers
end
function coral_cover(rs::ResultSet)::NamedTuple
    return coral_cover(rs.raw)
end


"""
    coral_evenness(rs::ResultSet)

Calculates evenness across functional coral groups in ADRIA.
Inverse Simpsons diversity indicator.

# Notes
Number of taxa (distinct groups with enhanced lumped with unenhanced) is hardcoded in this function.
"""
function coral_evenness(rs::ResultSet)::AbstractArray{<:Real}
    X = rs.raw
    return coral_evenness(X)
end
function coral_evenness(X::AbstractArray{<:Real})::AbstractArray{<:Real}
    x = min.(max.(X, 0.0), 1.0)
    covers = relative_cover(x)

    # Evenness as a functional diversity metric
    n = 4;  # number of taxa
    p1 = dropdims(sum(covers.tab_acr, dims=2), dims=2) ./ covers;
    p2 = dropdims(sum(covers.cor_acr, dims=2), dims=2) ./ covers;
    p3 = dropdims(sum(covers.small_enc, dims=2), dims=2) ./ covers;
    p4 = dropdims(sum(covers.large_mass, dims=2), dims=2) ./ covers;

    sum_psqr = p1.^2 + p2.^2 + p3.^2 + p4.^2;  # functional diversity
    simpson_D = 1 ./ sum_psqr;  # Hill 1973, Ecology 54:427-432
    return simpson_D ./ n;  # Group evenness
end


"""
    shelter_volume(rs::ResultSet)
    shelter_volume(X::AbstractArray, inputs::DataFrame)

Provide indication of shelter volume.
"""
function shelter_volume(X::AbstractArray{<:Real}, inputs::DataFrame)::AbstractArray{<:Real}
    _, _, cs_p = coral_spec()
    n_corals = length(unique(cs_p.taxa_id))

    colony_area_cm2 = Array(inputs[:, contains.(names(inputs), "colony_area_cm2")])'

    sheltervolume_parameters = [
        -8.32 1.50;   # tabular from Urbina-Barretto 2021
        -8.32 1.50;   # tabular from Urbina-Barretto 2021
        -7.37 1.34;   # columnar from Urbina-Barretto 2021, assumed similar for corymbose Acropora
        -7.37 1.34;   # columnar from Urbina-Barretto 2021, assumed similar for corymbose Acropora
        -9.69 1.49;   # massives from Urbina-Barretto 2021, assumed similar for encrusting and small massives
        -9.69 1.49]  # massives from Urbina-Barretto 2021,  assumed similar for large massives

    sheltervolume_parameters = repeat(sheltervolume_parameters, n_corals, 1)

    ntsteps, nspecies, nsites, nint, nreps = size(X);

    #  Estimate log colony volume (litres) based on relationship
    #  established by Urbina-Barretto 2021
    logcolony_sheltervolume = sheltervolume_parameters[:,1] .+ sheltervolume_parameters[:,2] .* log10.(colony_area_cm2);
    maxlogcolony_sheltervolume = sheltervolume_parameters[:,1] .+ sheltervolume_parameters[:,2] .* log10.(maximum(colony_area_cm2, dims=1));

    shelter_volume_colony_litres_per_cm2 = (10.0.^logcolony_sheltervolume);
    max_shelter_volume_colony_litres_per_cm2 = (10.0.^maxlogcolony_sheltervolume);

    # convert from litres per cm2 to m3 per ha
    cm2_m3 = (10^-3) * 10^4 *10^4
    shelter_volume_colony_m3_per_ha = shelter_volume_colony_litres_per_cm2 * cm2_m3;
    max_shelter_volume_colony_m3_per_ha = max_shelter_volume_colony_litres_per_cm2 * cm2_m3;

    # calculate shelter volume of groups and size classes and multiply with covers
    sv = zeros(ntsteps, nspecies, nsites, nint, nreps);
    for sp = 1:nspecies
        sv[:,sp,:,:,:] = (shelter_volume_colony_m3_per_ha[sp] / max_shelter_volume_colony_m3_per_ha[sp]) .* X[:,sp,:,:,:];
    end

    #  sum over groups and size classes to estimate total shelter volume per ha
    return dropdims(sum(sv, dims=2), dims=2);
end
function shelter_volume(rs::ResultSet)::AbstractArray{<:Real}
    return shelter_volume(rs.raw, rs.inputs)
end


"""
    reef_condition_index(TC, E, SV, juveniles)
    reef_condition_index(rs)

Translates coral metrics in ADRIA to a Reef Condition Metrics

# Inputs
- TC        : Total relative coral cover across all groups
- E         : Evenness across four coral groups
- SV        : Shelter volume based coral sizes and abundances
- juveniles : Abundance of coral juveniles < 5 cm diameter

Input dimensions: timesteps, species, sites

# Outputs
Dimensions: timesteps, sites, interventions, repeats
"""
function reef_condition_index(TC::AbstractArray{<:Real}, E::AbstractArray, SV::AbstractArray{<:Real}, juveniles::AbstractArray{<:Real})::Array{<:Real}
    # Compare outputs against reef condition criteria provided by experts

    # These are median values for 7 experts. TODO: draw from distributions
    #  Condition        TC       E       SV      Juv
    # {'VeryGood'}      0.45     0.45    0.45    0.35
    # {'Good'    }      0.35     0.35    0.35    0.25
    # {'Fair'    }      0.25     0.25    0.30    0.25
    # {'Poor'    }      0.15     0.25    0.30    0.25
    # {'VeryPoor'}      0.05     0.15    0.18    0.15

    # Note that the scores for evenness and juveniles are slightly different
    lin_grid = Gridded(Linear())
    TC_func = interpolate(([0, 0.05, 0.15, 0.25, 0.35, 0.45, 1.0],), [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0], lin_grid)
    # E_func = interpolate(([0, 0.15, 0.25, 0.35, 0.45, 1.0],), [0, 0.1, 0.5, 0.7, 0.9, 1.0], lin_grid)
    SV_func = interpolate(([0, 0.18, 0.30, 0.35, 0.45, 1.0],), [0, 0.1, 0.3, 0.5, 0.9, 1.0], lin_grid)
    juv_func = interpolate(([0, 0.15, 0.25, 0.35, 1.0],), [0, 0.1, 0.5, 0.9, 1.0], lin_grid)

    TC_i = TC_func.(TC);
    # E_i = E_func.(E);
    SV_i = SV_func.(SV);
    juv_i = juv_func.(juveniles);

    # Original
    # Y = (TC_i + E_i + SV_i + juv_i) ./ 4;

    # Weighted, giving evenness 10#  weight
    # Y = (TC_i*0.3) + (E_i*0.1) + (SV_i*0.3) + (juv_i*0.3);

    # Removing evenness completely
    # Y = mean([TC_i, SV_i, juv_i])
    # Y = (TC_i .+ SV_i .+ juv_i) ./ 3;

    return mean([TC_i, SV_i, juv_i])
end
function reef_condition_index(rs::ResultSet)::Array{<:Real}
    cover = coral_cover(rs)
    TC = cover.total_cover
    juv = cover.juveniles
    E = coral_evenness(rs)
    SV = shelter_volume(rs)
    return reef_condition_index(TC, E, SV, juv)
end


function summarize_total_cover(rs::ResultSet)::NamedTuple
    cover = coral_cover(rs).total_cover

    x = Dict(Symbol(f) => dropdims(f(cover, dims=(3,2)), dims=(3,2)) for f in [mean, median, std, minimum, maximum])
    return (;zip(keys(x), values(x))...)
end

end