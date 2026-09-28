using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ADRIA

domain_path = normpath(joinpath(@__DIR__, "..", "data", "Lizard_Historical_v0.1"))
domain = ADRIA.load_domain(ADRIA.LizardDomain, domain_path, "historical")

@assert size(domain.conn) == (2914, 2914)
@assert size(domain.cots_conn) == size(domain.conn)
@assert domain.cots_conn.Source == domain.conn.Source
@assert domain.cots_conn.Sink == domain.conn.Sink
@assert maximum(abs.(domain.cots_conn.data .- domain.conn.data)) > 0.0

println(
    "Lizard COTS connectivity load validated; max=",
    maximum(domain.cots_conn.data),
    ", nonzero=",
    count(!iszero, domain.cots_conn.data)
)
