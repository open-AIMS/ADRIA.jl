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
function clustered_scenarios() end
function clustered_scenarios!() end

# Rule extraction
function rules_scatter() end
function rules_scatter!() end

# Spatial
function map() end
function map!() end

end  # module
