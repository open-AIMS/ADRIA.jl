module ADRIA

using Random, TOML, Dates, CpuId
using StaticArrays, SparseArrays, LinearAlgebra, Statistics, Distributed
using NamedArrays, SparseArrayKit, DifferentialEquations

using MAT  # Package to read in `.mat` files

using Setfield, ModelParameters, DataStructures
using DataFrames, GeoDataFrames, Graphs, CSV
import ArchGDAL as AG

using PkgVersion

using ProgressMeter


include("utils/text_display.jl")  # need better name for this file
include("utils/setup.jl")

include("ecosystem/corals/growth.jl")
include("ecosystem/corals/CoralGrowth.jl")
include("ecosystem/Ecosystem.jl")

# Generate base coral struct from default spec.
# Have to call this before including specification methods
create_coral_struct()

include("ecosystem/corals/spec.jl")
include("ecosystem/const_params.jl")

include("Domain.jl")
include("io/inputs.jl")

include("sites/connectivity.jl")
include("sites/dMCDA.jl")

include("io/ResultSet.jl")
include("metrics/metrics.jl")
include("io/result_io.jl")
include("io/result_post_processing.jl")

include("scenario.jl")

# include("main_app.jl")


export fecundity_scope!, bleaching_mortality!
export growthODE
export run_scenario, coral_spec
export create_coral_struct, Intervention, Criteria, Corals, SimConstants
export site_area
export Domain, metrics, select

# metric helper methods
export dims, ndims


if ccall(:jl_generating_output, Cint, ()) == 1   # if we're precompiling the package
    precompile(load_results, (String, ))
    precompile(load_domain, (String, Int64))
    precompile(Domain, (String, Int64, String, String, String, String, String, String, String))
    precompile(EnvLayer, (String, String, String, String, String, String, String))
end


# Precompile as the final step of the module definition:
# if ccall(:jl_generating_output, Cint, ()) == 1   # if we're precompiling the package
#     precompile(load_domain, (String, Int64))
#     precompile(Domain, (String, Int64, String, String, String, String, String, String, String))
#     precompile(EnvLayer, (String, String, String, String, String, String, String))

#     precompile(Domain, (String, Int, EnvLayer, DataFrame, Vector{Float64}, Vector{Int64}, DataFrame, String, String, NamedMatrix, CoralGrowth,
#         Vector{String}, Vector{String}, NamedArray, NamedArray))

#     # let
#     #     here = @__DIR__
#     #     ex_dir = joinpath(here, "../examples")
#     #     @debug "Pre-running examples to reduce future spin-up time"

#     #     f() = begin 
#     #         @showprogress 1 for _ in 1:10
#     #         end
#     #     end
#     #     b = redirect_stdout(f, devnull);

#     #     ex_domain = ADRIA.load_domain(joinpath(ex_dir, "Example_domain"), 45)
#     #     p_df = ADRIA.load_scenarios(ex_domain, joinpath(ex_dir, "example_scenarios.csv"))

#     #     ENV["ADRIA_THRESHOLD"] = 1e-6
#     #     ex_domain.sim_constants.tf = 3
#     #     ds = (raw=nothing, site_ranks=nothing, seed_log=nothing, fog_log=nothing, shade_log=nothing)
#     #     run_scenario((1, p_df[1, :]), ex_domain, 1, ds)
#     #     run_scenario((1, p_df[end, :]), ex_domain, 1, ds)
#     #     delete!(ENV, "ADRIA_THRESHOLD")
#     # end
# end

end
