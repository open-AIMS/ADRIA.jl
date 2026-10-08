using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ADRIA

domain_version = get(ENV, "LIZARD_DOMAIN_VERSION", "Lizard_Historical_v0.1")
domain_path = normpath(joinpath(@__DIR__, "..", "data", domain_version))
domain = ADRIA.load_domain(ADRIA.LizardDomain, domain_path, "historical")

site_count = length(domain.loc_ids)
@assert size(domain.conn) == (site_count, site_count)
@assert size(domain.cots_conn) == size(domain.conn)
@assert domain.cots_conn.Source == domain.conn.Source
@assert domain.cots_conn.Sink == domain.conn.Sink
@assert maximum(abs.(domain.cots_conn.data .- domain.conn.data)) > 0.0
@assert !isnothing(domain.cots_forcing)
@assert length(domain.cots_forcing.years) == 6
@assert length(domain.cots_forcing.reef_connectivity) == 6
@assert length(domain.cots_forcing.site_to_reef) == site_count
@assert size(domain.cots_forcing.source_survival) == (113, 6)
@assert size(domain.cots_forcing.external_supply) == (113, 6)
@assert length(domain.cots_forcing.reef_habitable_area_ha) == 113
@assert all(domain.cots_forcing.reef_habitable_area_ha .> 0.0)
@assert domain.cots_forcing.outbreak_years == collect(1991:2025)
@assert size(domain.cots_forcing.external_outbreak_potential) == (113, 35, 6)
@assert size(domain.cots_forcing.external_outbreak_pelagic) == (113, 35, 6)
@assert all(0.0 .<= domain.cots_forcing.source_survival .<= 1.0)
@assert all(domain.cots_forcing.external_supply .>= 0.0)
@assert all(domain.cots_forcing.external_outbreak_potential .>= 0.0)
@assert all(domain.cots_forcing.external_outbreak_pelagic .>= 0.0)
@assert all(
    domain.cots_forcing.external_outbreak_pelagic .<=
    domain.cots_forcing.external_outbreak_potential
)

println(
    "Lizard COTS connectivity load validated; max=",
    maximum(domain.cots_conn.data),
    ", nonzero=",
    count(!iszero, domain.cots_conn.data),
    ", annual years=",
    join(domain.cots_forcing.years, ',')
)
