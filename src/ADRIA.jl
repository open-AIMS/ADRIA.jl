module ADRIA

using Random, TOML, Hwloc, StaticArrays, SparseArrays, LinearAlgebra, Statistics
using Distributed

using MATLAB, MAT  # MATLAB interface and package to read in `.mat` files

using NamedArrays, SparseArrayKit, DifferentialEquations
using Setfield, ModelParameters, DataStructures
using DataFrames, GeoDataFrames, Graphs

using CSV


include("utils/text_display.jl");  # need better name for this file

include("ecosystem/corals/growth.jl");
include("ecosystem/corals/CoralGrowth.jl");
include("ecosystem/Ecosystem.jl");

# Generate base coral struct from default spec.
# Have to call this before including specification methods
create_coral_struct()

include("ecosystem/corals/spec.jl")
include("ecosystem/const_params.jl")

include("Domain.jl")
include("io/inputs.jl")
include("io/result_set.jl")
include("scenario.jl")

include("sites/connectivity.jl")
include("sites/dMCDA.jl")

include("metrics/metrics.jl")

include("main_app.jl")


export fecundity_scope!, bleaching_mortality!
export growthODE
export run_scenario, coral_spec
export create_coral_struct, Intervention, Criteria, Corals, SimConstants
export Domain, metrics, select

end
