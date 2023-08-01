"""
    run_scenarios(scens::DataFrame, dom::ReefModDomain, rcp::String)

Run ADRIA coral ecosystem model with ReefMod data.
"""
function run_scenarios(::Type{ReefModDomain}, scens::DataFrame, dom::ReefModDomain, rcp::String)

    rs = run_scenarios(scens, dom, rcp)

    return rs
end
