@doc """
Functions for post-processing ADRIA results into excel/csv
format suitable for PowerBI.

Argument: res::ResultSet, a set of results run in ADRIA
Output: a csv file with columns [Model, SSP, SiteID, Latitude, Longitude,
                                Year, DeployYear, SeedLevel, DHWenhancement, 
                                FogLevel, SiteArea, Kvalue, Guided,
                                MeanCoralCoverProp, DiffMeanCoralCoverProp,
                                ShelterVolume, DiffShelterVolume]
"""

using ADRIA
import ADRIA: ResultSet
using Statistics
using CSV
using DataFrames

function create_excel_format_file(res::ResultSet)

    model = "ADRIA"
    # extract RCP and convert to SSP
    rcp = res.RCP
    if rcp == "26"
        ssp = "SSP1"
    elseif rcp == "45"
        ssp = "SSP2"
    elseif rcp =="60"
        ssp = "SSP3"
    end

    # calculate relative cover and shelter vol (other metrics to come)
    rel_cover = ADRIA.metrics.relative_cover(res)
    #juveniles = dropdims(mean(cover[:juveniles],dims=:reps),dims=:reps)
    #rci = dropdims(mean(ADRIA.metrics.reef_condition_index(res),dims=:reps),dims=:reps)
    sheltervol = ADRIA.metrics.relative_shelter_volume(res)

    n_sites = size(rel_cover)[2]
    n_scens = size(rel_cover)[3]

    # extract key site data
    site_ids = res.site_ids
    centroids = res.site_centroids
    kvals = res.site_max_coral_cover
    sitearea = res.site_area

    # set up for extract at 5-yearly slices
    tf = res.sim_constants["tf"]
    years = collect(2025:5:2025+tf)
    years_ints = collect(1:5:size(rel_cover)[1])
    n_years = length(years_ints)

    # storage dataframe
    data_sum_df = DataFrame(Model = String[], SSP = String[],SiteID = String[],Latitude = Float64[],Longitude=Float64[],
                        Year=Int[],DeployYear=Int[],SeedLevel=Float64[],DHWenhancement=Float64[],FogLevel=Int[],
                        SiteArea=Float64[],Kvalue=Float64[],Guided=Int[],
                        MeanCoralCoverProp=Float32[],DiffMeanCoralCoverProp=Float32[],
                        ShelterVolume=Float32[],DiffShelterVolume=Float32[])
                        #,Juveniles=Float64[],DiffJuveniles=Float64[],
                        #RCI=Float64[],DiffRCI=Float64[])

    # find conterfactual index
    cond =((res.inputs.guided.==0) .& (res.inputs.seed_TA.==0) .& (res.inputs.seed_CA.==0) .& (res.inputs.fogging.==0) .& (res.inputs.SRM.==0) .& (res.inputs.a_adapt.==0) .& (res.inputs.n_adapt.==0))
    conter_ind = rownumber.(eachrow(res.inputs[cond,:]))

    for t in collect(1:n_years)
        for si in collect(1:n_sites)
            for sce in collect(1:n_scens)           

                # guided or unguided
                res.inputs.guided[sce] > 0 ? guided = 1 : guided = 0
                # seeding level including both species
                seed = res.inputs.seed_TA[sce] + res.inputs.seed_CA[sce]
                # fogging or no fogging
                res.inputs.fogging[sce]>0 ? fog = 1 : fog = 0

                # add scenario to structure
                push!(data_sum_df, (model, ssp, site_ids[si], centroids[si][2,1],centroids[si][1,1], years[t],Int(res.inputs.seed_year_start[sce]+2024),
                                    seed, res.inputs.a_adapt[sce], fog, sitearea[si], kvals[si], guided,
                                    rel_cover[years_ints[t],si,sce]*100,(rel_cover[years_ints[t],si,sce] .- rel_cover[years_ints[t],si,conter_ind])[1]*100,
                                    sheltervol[years_ints[t],si,sce]*100,(sheltervol[years_ints[t],si,sce] .- sheltervol[years_ints[t],si,conter_ind])[1]*100))
                                        #juveniles[years_ints[t],si,sce],juveniles[years_ints[t],si,sce]-juveniles[years_ints[t],si,0],
                                        #rci[years_ints[t],si,sce],rci[years_ints[t],si,sce]-rci[years_ints[t],si,0]))
                            
            end
        end
    end
    

    cluster = res.name
    CSV.write(cluster*"_"*ssp*"_data_summary.csv", data_sum_df)
end

