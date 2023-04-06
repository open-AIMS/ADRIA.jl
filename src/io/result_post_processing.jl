@doc """
Functions for post-processing ADRIA results.
"""

using ADRIA
using ADRIA: ResultSet
using Statistics
using CSV
using DataFrames

"""
    RCP_to_SSP(rcp)

Convert RCP scenarios to SSP scenarios.

# Arguments: 
- `rcp::String`: RCP scenario identifier
"""
function RCP_to_SSP(rcp::String)::String
    if rcp == "26"
        ssp = "SSP1"
    elseif rcp == "45"
        ssp = "SSP2"
    elseif rcp == "60"
        ssp = "SSP3"
    else
        throw("Unknown RCP value")
    end

    return ssp
end

"""
    create_BI_format_file(rs, file_loc)

Creates a tabular csv files from ADRIA ResultSet, suitable for processing in PowerBI etc.
Saves as csv at `file_loc`.

# Arguments:
- `rs::ResultSet`: ADRIA result set
- `file_loc::String`: directory to save the csv file in
"""
function create_BI_format_file(rs::ResultSet, file_loc::String)
    model = "ADRIA"

    # extract RCP and convert to SSP
    rcp = rs.RCP
    ssp = RCP_to_SSP(rcp)

    # calculate relative cover and shelter vol (other metrics to come)
    rel_cover = ADRIA.metrics.relative_cover(rs)
    #juveniles = dropdims(mean(cover[:juveniles],dims=:reps),dims=:reps)
    #rci = dropdims(mean(ADRIA.metrics.reef_condition_index(res),dims=:reps),dims=:reps)
    sheltervol = ADRIA.metrics.relative_shelter_volume(rs)

    tf, n_locations, n_scens = size(rel_cover)

    # extract key location data
    location_ids = rs.location_ids
    centroids = rs.location_centroids
    kvals = rs.location_max_coral_cover
    locationarea = rs.location_area

    # set up for extract at 5-yearly slices
    years = collect(2025:5:2025+tf)
    years_ints = collect(1:5:tf)
    n_years = length(years_ints)

    perms = n_years * n_locations * n_scens
    # storage dataframe
    data_sum_df = DataFrame(Model=repeat([""], outer=perms), SSP=repeat([""], outer=perms),
        SiteID=repeat([""], outer=perms), Latitude=repeat([0.0], outer=perms),
        Longitude=repeat([0.0], outer=perms), Year=repeat([0], outer=perms),
        DeployYear=repeat([0], outer=perms), SeedLevel=repeat([0.0], outer=perms),
        DHWenhancement=repeat([0.0], outer=perms), FogLevel=repeat([0], outer=perms),
        SiteArea=repeat([0.0], outer=perms), Kvalue=repeat([0.0], outer=perms),
        Guided=repeat([0], outer=perms), MeanCoralCoverProp=repeat([parse(Float32, "0")], outer=perms),
        DiffMeanCoralCoverProp=repeat([parse(Float32, "0")], outer=perms),
        ShelterVolume=repeat([parse(Float32, "0")], outer=perms),
        DiffShelterVolume=repeat([parse(Float32, "0")], outer=perms))
    #,Juveniles=Float64[],DiffJuveniles=Float64[],
    #RCI=Float64[],DiffRCI=Float64[])

    # find conterfactual index
    cond = ((rs.inputs.guided .== 0) .& (rs.inputs.seed_TA .== 0) .& (rs.inputs.seed_CA .== 0) .& (rs.inputs.fogging .== 0) .& (rs.inputs.SRM .== 0) .& (rs.inputs.a_adapt .== 0) .& (rs.inputs.n_adapt .== 0))
    counter_ind = rownumber.(eachrow(rs.inputs[cond, :]))
    count = 1
    for t in collect(1:n_years)
        for si in collect(1:n_locations)
            for sce in collect(1:n_scens)

                # guided or unguided
                rs.inputs.guided[sce] > 0 ? guided = 1 : guided = 0
                # seeding level including both species
                seed = rs.inputs.seed_TA[sce] + rs.inputs.seed_CA[sce]
                # fogging or no fogging
                rs.inputs.fogging[sce] > 0 ? fog = 1 : fog = 0

                # add scenario to structure
                data_sum_df[count, :] = (model, ssp, location_ids[si], centroids[si][2, 1], centroids[si][1, 1], years[t], Int(rs.inputs.seed_year_start[sce] + 2024),
                    seed, rs.inputs.a_adapt[sce], fog, locationarea[si], kvals[si], guided,
                    rel_cover[years_ints[t], si, sce] * 100, (rel_cover[years_ints[t], si, sce].-rel_cover[years_ints[t], si, counter_ind])[1] * 100,
                    sheltervol[years_ints[t], si, sce] * 100, (sheltervol[years_ints[t], si, sce].-sheltervol[years_ints[t], si, counter_ind])[1] * 100)
                #juveniles[years_ints[t],si,sce],juveniles[years_ints[t],si,sce]-juveniles[years_ints[t],si,0],
                #rci[years_ints[t],si,sce],rci[years_ints[t],si,sce]-rci[years_ints[t],si,0]))
                count = count + 1
            end
        end
    end


    cluster = rs.name
    file_loc = replace(file_loc, "\\" => "/")
    CSV.write(file_loc * "/" * "data_summary_" * "ADRIA" * rs.ADRIA_VERSION * "_" * rs.invoke_time * "_" * cluster * "_" * ssp * ".csv", data_sum_df)
end

