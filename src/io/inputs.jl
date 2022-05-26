"""
    load_domain(path::String, rcp::Union{String, Integer})

Load domain specification from data package.

# Arguments
- path : location of data package
- rcp : RCP scenario to run
"""
function load_domain(path::String, rcp::Union{String, Integer})
    domain_name = basename(path)
    conn_path = joinpath(path, "connectivity/")
    site_data = joinpath(path, "site_data")

    dhw = joinpath(path, "DHWs", "dhwRCP$(rcp).mat")

    site_path = joinpath(site_data, "$(domain_name).gpkg")
    init_coral_cov = joinpath(site_data, "coral_cover.mat")
    wave = joinpath(path, "waves/wave_data.mat")

    dom = Domain(
        domain_name,
        rcp,
        site_path,
        "siteref",
        "reef_siteid",
        init_coral_cov,
        conn_path,
        dhw,
        wave
    );

    return dom
end
