using ADRIA
using ADRIA: site_selection
using DataFrames

@info "Loading data package"
here = @__DIR__
ex_domain = ADRIA.load_domain(joinpath(here, "Example_domain"))

@info "Create criteria structure using scenario DataFrame."
criteria_df = ADRIA.sample(ex_domain, 5) # get scenario dataframe

area_to_seed = 1.5 * 10^-6 # area of seeded corals in km^2
nreps = 30 # number of dhw and wave replicates you want to use
ts = 5 # time step to perform site selection at
depth = DataFrame(depth_min=5, depth_offset=10)

cover = 0.05 .* ones(sum(criteria_df.guided .> 0), 36, size(ex_domain.site_data, 1)) # cover has size scenarios * species * sites
ranks = site_selection_scens(ex_domain, criteria_df[criteria_df.guided.>0, :], depth, cover, area_to_seed, ts, n_reps)
