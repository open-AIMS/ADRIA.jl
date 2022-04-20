module ADRIA

using MATLAB  # MATLAB interface
using MAT     # Julia package to read `.mat` files

using StaticArrays
using LinearAlgebra, Statistics
using DifferentialEquations

using Setfield, ModelParameters

using DataFrames, GeoDataFrames, Graphs

using CSV


include("utils/text_display.jl");  # need better name for this file

include("ecosystem/corals/growth.jl");
include("ecosystem/corals/fecundity_scope.jl");
include("ecosystem/corals/bleaching_mortality.jl");
include("ecosystem/corals/CoralGrowth.jl");
include("ecosystem/Ecosystem.jl");

# Generate base coral struct from default spec.
# Have to call this before including specification methods
create_coral_struct()

include("ecosystem/corals/spec.jl");
include("ecosystem/const_params.jl");

include("Domain.jl")
include("scenario.jl");

# Connectivity
include("sites/connectivity.jl")


export fecundity_scope!, bleaching_mortality!
export growthODE
export run_scenario, coral_spec
export create_coral_struct, Intervention, Criteria, Corals, SimConstants
export Domain

end
