module viz

"""
Dummy interfaces for extensions to hook into.

Package extension precompilation fails if these are not defined
as extensions are expected to provide additional functionality
for existing functions/methods.
"""

# GUI
function explore() end

# Scenario plotting methods
function scenario() end
function scenario!() end

# Sensitivity analyses
function pawn() end
function pawn!() end
function tsa() end
function tsa!() end
function rsa() end
function rsa!() end

function outcome_map() end
function outcome_map!() end

# Clustering
function ts_cluster() end
function ts_cluster!() end

end  # module
