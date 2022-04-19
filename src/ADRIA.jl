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
include("ecosystem/Ecosystem.jl");

# Generate base coral struct from default spec.
# Have to call this before including specification methods
create_coral_struct()

include("ecosystem/corals/spec.jl");
include("ecosystem/const_params.jl");

include("scenario.jl");


struct Domain
    TP_data     # site connectivity data
    site_ranks  # site rank
    strongpred  # strongest predecessor
    site_data   # table of site data (depth, carrying capacity, etc)
    init_coral_cover  # initial coral cover dataset
    connectivity_site_ids  # Site IDs as specified by the connectivity dataset (indicates order of `TP_data`)
    removed_sites  # indices of sites that were removed. Used to align site_data, DHW, connectivity, etc.
    dhw_scens  # DHW scenarios
    wave_scens # wave scenarios
end


export fecundity_scope!, bleaching_mortality
export growthODE
export run_scenario, coral_spec, to_spec
export create_coral_struct, Intervention, Criteria, Corals, SimConstants
export Domain

end
