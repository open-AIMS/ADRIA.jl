using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using ADRIA
using SHA

const DATA_ROOT = normpath(joinpath(@__DIR__, "..", "data"))
const V1_ROOT = joinpath(DATA_ROOT, "Lizard_Historical_v0.1")
const V2_ROOT = joinpath(DATA_ROOT, "Lizard_Historical_v0.2")

function sha256_file(path::String)::String
    return bytes2hex(open(SHA.sha256, path))
end

previous_mode = get(ENV, "ADRIA_COTS_CONNECTIVITY_MODE", nothing)
previous_dataset = get(ENV, "ADRIA_COTS_CONNECTIVITY_DATASET", nothing)
ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = "cots"
ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = "owen_global_2018_2023"
try
    v2 = ADRIA.load_domain(ADRIA.LizardDomain, V2_ROOT, "historical")
    @assert length(v2.loc_ids) == 2895
    @assert v2.loc_ids == ["FNG_V2_$(lpad(index, 4, '0'))" for index in 1:2895]
    @assert length(unique(v2.loc_ids)) == 2895
    @assert length(unique(string.(v2.loc_data.UNIQUE_ID))) == 113
    @assert size(v2.conn) == (2895, 2895)
    @assert string.(v2.conn.Source.val.data) == v2.loc_ids
    @assert string.(v2.conn.Sink.val.data) == v2.loc_ids
    @assert all(isfinite, v2.conn.data)
    @assert minimum(v2.conn.data) >= 0.0
    @assert maximum(abs.(v2.conn.data .- transpose(v2.conn.data))) <= 1e-12
    @assert isapprox(sum(v2.loc_data.area), 408_769_510.674482; rtol=1e-10)
    @assert all(v2.loc_data.area .> 0.0)
    @assert size(v2.init_coral_cover) == (35, 2895)
    @assert size(v2.dhw_scens) == (40, 2895, 2)
    @assert maximum(v2.dhw_scens.data) > 0.0
    @assert maximum(abs.(v2.cots_conn.data .- v2.conn.data)) > 0.0
    @assert !isnothing(v2.cots_forcing)
    @assert v2.cots_forcing.years == collect(2018:2022)
    @assert length(v2.cots_forcing.reef_connectivity) == 5
    @assert length(v2.cots_forcing.site_to_reef) == 2895
    @assert size(v2.cots_forcing.source_survival) == (113, 5)
    @assert size(v2.cots_forcing.external_supply) == (113, 5)
    @assert v2.cots_forcing.outbreak_years == collect(1991:2025)
    @assert size(v2.cots_forcing.external_outbreak_potential) == (113, 35, 5)
    @assert size(v2.cots_forcing.external_outbreak_pelagic) == (113, 35, 5)
    @assert all(isfinite, v2.cots_forcing.external_outbreak_potential)
    @assert all(isfinite, v2.cots_forcing.external_outbreak_pelagic)

    # Loading V2 must not require or mutate the retained V1 control package.
    v1_provenance = joinpath(V1_ROOT, "cots_connectivity", "provenance.toml")
    v1_hash_before = sha256_file(v1_provenance)
    ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = "legacy"
    v1 = ADRIA.load_domain(ADRIA.LizardDomain, V1_ROOT, "historical")
    @assert length(v1.loc_ids) == 2914
    @assert sha256_file(v1_provenance) == v1_hash_before

    println(
        "Validated Lizard_Historical_v0.2: sites=", length(v2.loc_ids),
        ", reefs=", length(unique(string.(v2.loc_data.UNIQUE_ID))),
        ", area_m2=", sum(v2.loc_data.area),
        ", max_dhw=", maximum(v2.dhw_scens.data),
        ", v1_sites=", length(v1.loc_ids)
    )
finally
    if isnothing(previous_mode)
        delete!(ENV, "ADRIA_COTS_CONNECTIVITY_MODE")
    else
        ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = previous_mode
    end
    if isnothing(previous_dataset)
        delete!(ENV, "ADRIA_COTS_CONNECTIVITY_DATASET")
    else
        ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = previous_dataset
    end
end
