"""
    run_scenarios(scens::DataFrame, dom::ReefModDomain, rcp::String)

Run ADRIA coral ecosystem model with ReefMod data.
"""
function run_scenarios(::Type{ReefModDomain}, scens::DataFrame, dom::ReefModDomain, rcp::String)
    # Convert ReefMod initial coral cover data into ADRIA equivalent

    # Values are relative to `k` area
    # Columns are 
    # - Enhanced Tabular Acropora
    # - Unenhanced Tabular Acropora, 
    # - Enhanced Corymbose Acropora, 
    # - Unenhanced Corymbose Acropora 
    # - Smalll Massives
    # - Large Massives
    # Values interpolated by K. Anthony.
    cover_prop = Array{Float64,2}([
        0.0 0.000257478 0.0 0.000261732 0.000160465 0.000161798
        0.0 0.001090911 0.0 0.001162885 0.001287053 0.00116821
        0.0 0.001947037 0.0 0.003690998 0.00727025 0.005169573
        0.0 0.011468529 0.0 0.019751902 0.029379472 0.024928862
        0.0 0.008814995 0.0 0.019410142 0.045405627 0.020003117
        0.0 0.031648876 0.0 0.029897563 0.05573267 0.060801976
    ])

    covers = repeat(cover_prop[:], 1, n_locations(dom))
    covers = covers .* (site_k_area(dom) ./ site_area(dom))'
    covers = NamedDimsArray(covers, species=1:size(covers, 1), sites=dom.site_ids)

    orig_cover = dom.init_coral_cover[:, :]  # take a copy of the original data set

    dom.init_coral_cover = covers
    rs = run_scenarios(scens, dom, rcp)
    dom.init_coral_cover = orig_cover

    return rs
end