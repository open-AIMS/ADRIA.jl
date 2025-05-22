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
function scenarios() end
function scenarios!() end
function scenarios_legend!() end

# Plot a single scenario
function scenario_stratified() end

# Sensitivity analyses
function pawn() end
function pawn!() end
function tsa() end
function tsa!() end
function rsa() end
function rsa!() end
function convergence() end
function convergence!() end
function outcome_map() end
function outcome_map!() end

# Clustering
function clustered_scenarios() end
function clustered_scenarios!() end

# Rule extraction
function rules_scatter() end
function rules_scatter!() end

# Location selection plots
function ranks_to_frequencies() end
function ranks_to_frequencies!() end
function selection_criteria_map() end
function selection_criteria_map!() end

# Spatial
function map() end
function map!() end
function connectivity() end
function connectivity!() end

# Environment
function dhw_scenario() end
function dhw_scenarios() end
function cyclone_scenario() end

# Coral Dynamics
function taxonomy() end
function taxonomy!() end

# Economics
function data_envelopment_analysis() end
function data_envelopment_analysis!() end

end  # module
