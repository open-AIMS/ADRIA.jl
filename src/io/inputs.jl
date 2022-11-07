using JSON


"""
    _check_compat(dpkg_details::Dict)

Checks for version compatibility.

# Arguments
- dpkg_details : Datapackage spec
"""
function _check_compat(dpkg_details::Dict)
    if haskey(dpkg_details, "version") || haskey(dpkg_details, "dpkg_version")
        dpkg_version = dpkg_details["version"]
        if dpkg_version ∉ COMPAT_DPKG
            error("Incompatible Domain data package.")
        end
    else
        error("Incompatible Domain data package.")
    end
end

"""
    _load_dpkg(dpkg_path::String)

Load and parse datapackage.

# Arguments
- dpkg_path : path to datapackage
"""
function _load_dpkg(dpkg_path::String)
    dpkg_md = nothing
    open(joinpath(dpkg_path, "datapackage.json"), "r") do fp
        dpkg_md = JSON.parse(read(fp, String))
    end
    _check_compat(dpkg_md)

    return dpkg_md
end

"""
    load_domain(path::String, rcp::Int64)
    load_domain(path::String, rcp::String)
    load_domain(path::String)

Load domain specification from data package.

# Arguments
- path : location of data package
- rcp : RCP scenario to run. If none provided, no data path is set.
"""
function load_domain(path::String, rcp::String)::Domain
    domain_name::String = basename(path)
    if length(domain_name) == 0
        domain_name = basename(dirname(path))
    end

    dpkg_details = _load_dpkg(path)
    dpkg_version = dpkg_details["version"]

    # Handle compatibility
    this_version = parse(VersionNumber, dpkg_version)
    if this_version >= v"0.2.1"
        # Extract the time frame represented in this data package
        timeframe = dpkg_details["simulation_metadata"]["timeframe"]
    else
        # Default to 2025-2099
        timeframe = (2025, 2099)
    end

    if length(timeframe) == 2
        @assert timeframe[1] < timeframe[2] "Start date/year specified in data package must be < end date/year"
        # If only two elements, assume a range is specified.
        # Collate the time steps as a full list if necessary
        timeframe = collect(timeframe[1]:timeframe[2])
    end

    conn_path::String = joinpath(path, "connectivity/")
    site_data::String = joinpath(path, "site_data")

    site_path::String = joinpath(site_data, "$(domain_name).gpkg")
    init_coral_cov::String = joinpath(site_data, "coral_cover.nc")

    if !isempty(rcp)
        dhw::String = joinpath(path, "DHWs", "dhwRCP$(rcp).nc")
        wave::String = joinpath(path, "waves", "wave_RCP$(rcp).nc")
    else
        dhw = ""
        wave = ""
    end

    return Domain(
        domain_name,
        path,
        rcp,
        timeframe,
        site_path,
        "site_id",
        "reef_siteid",
        init_coral_cov,
        conn_path,
        dhw,
        wave
    )
end
function load_domain(path::String, rcp::Int)::Domain
    return load_domain(path, "$rcp")
end
function load_domain(path::String)::Domain
    return load_domain(path, "")
end


"""
    load_scenarios(domain::Domain, filepath::String)::DataFrame

Load and pre-process scenario values.
Parameters intended to be of Integer type or casted as such.
"""
function load_scenarios(domain::Domain, filepath::String)::DataFrame
    df = CSV.read(filepath, DataFrame, comment="#")

    if "RCP" in names(df) || :RCP in names(df)
        df = df[!, Not("RCP")]
    end
    process_inputs!(domain, df)

    return df
end

function process_inputs!(d::Domain, df::DataFrame)
    bnds = d.model[:bounds]
    p_types = d.model[:ptype]
    @inbounds for (i, dt) in enumerate(p_types)
        if dt == "integer"
            df[!, i] .= map_to_discrete.(df[!, i], bnds[i][2])
        end
    end
end
