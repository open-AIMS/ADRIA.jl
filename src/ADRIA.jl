module ADRIA

using MATLAB

using StaticArrays
using LinearAlgebra
using DifferentialEquations

using Setfield, ModelParameters

using DataFrames, GeoDataFrames

using CSV


include("utils/text_display.jl");  # need better name for this file

include("ecosystem/corals/growth.jl");
include("ecosystem/corals/fecundity_scope.jl");
include("ecosystem/corals/bleaching_mortality.jl");
include("ecosystem/corals/spec.jl");
include("ecosystem/Ecosystem.jl");
include("ecosystem/const_params.jl");

include("scenario.jl");


export fecundity_scope!, bleaching_mortality
export growthODE
export run_scenario, coral_spec 
export CoralParams, Intervention, Criteria, SimConstants

end
