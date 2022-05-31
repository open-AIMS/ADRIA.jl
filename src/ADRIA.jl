module ADRIA

using Random, TOML, Dates, CpuId
using StaticArrays, SparseArrays, LinearAlgebra, Statistics, Distributed
using NamedArrays, SparseArrayKit, DifferentialEquations

using MAT  # Package to read in `.mat` files

using Setfield, ModelParameters, DataStructures
using DataFrames, GeoDataFrames, Graphs, CSV
using Plots

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


# Precompile as the final step of the module definition:
if ccall(:jl_generating_output, Cint, ()) == 1   # if we're precompiling the package
    let

        here = @__DIR__
        ex_dir = joinpath(here, "../examples")
        @debug "Pre-running examples to reduce future spin-up time"

        f() = begin 
            @showprogress 1 for _ in 1:10
            end
        end
        b = redirect_stdout(f, devnull);

        ex_domain = ADRIA.load_domain(joinpath(ex_dir, "Example_domain"), 45)
        p_df = CSV.read(joinpath(ex_dir, "example_scenarios.csv"), DataFrame, comment="#")

        ENV["ADRIA_THRESHOLD"] = 1e-6
        ex_domain.sim_constants.tf = 3
        ds = (raw=nothing, site_ranks=nothing, seed_log=nothing, fog_log=nothing, shade_log=nothing)
        run_scenario((1, p_df[1, :]), ex_domain, 1, ds)
        run_scenario((1, p_df[end, :]), ex_domain, 1, ds)
        delete!(ENV, "ADRIA_THRESHOLD")
    end
end

end
