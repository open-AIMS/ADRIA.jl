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

using SnoopPrecompile, RelocatableFolders

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

include("interventions/seeding.jl")

include("io/ResultSet.jl")
include("io/result_io.jl")
include("io/result_post_processing.jl")
include("metrics/metrics.jl")

include("scenario.jl")

# include("main_app.jl")


export fecundity_scope!, bleaching_mortality!
export growthODE
export run_scenario, coral_spec
export create_coral_struct, Intervention, Criteria, Corals, SimConstants
export site_area
export Domain, metrics, select, timesteps

# metric helper methods
export dims, ndims

# List out compatible domain datapackages
const COMPAT_DPKG = ["v0.2", "v0.2.1"]

@precompile_all_calls begin
    ex_dir = @path joinpath(@__DIR__, "../examples")

    f() = begin 
        @showprogress 1 for _ in 1:10
        end
    end
    b = redirect_stdout(f, devnull);

    ex_domain = ADRIA.load_domain(joinpath(ex_dir, "Example_domain"), 45)
    p_df = ADRIA.load_scenarios(ex_domain, joinpath(ex_dir, "example_scenarios.csv"))

    ENV["ADRIA_THRESHOLD"] = 1e-6
    run_scenario(p_df[1, :], ex_domain)
    run_scenario(p_df[end, :], ex_domain)
    delete!(ENV, "ADRIA_THRESHOLD")

    precompile(load_results, (String, ))
    precompile(EnvLayer, (String, String, String, String, String, String, String))
end

end
