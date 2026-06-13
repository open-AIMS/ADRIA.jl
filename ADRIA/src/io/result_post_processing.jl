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
