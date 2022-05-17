using PkgVersion
using Zarr

import Dates: now
import Setfield: @set!
import DataFrames: DataFrame
import ADRIA: Domain, EnvLayer


const RESULTS = "results"
const LOG_GRP = "logs"
const INPUTS = "inputs"


struct ResultSet
    name
    invoke_time
    env_layer_md

    inputs
    sim_constants

    raw
    ranks
    seed_log
    fog_log
    shade_log
end


"""
    setup_result_store!(domain::Domain, param_df::DataFrame)

Sets up an on-disk result store.

## Structure of Store

```
├───logs
│   ├───fog
│   ├───rankings
│   ├───seed
│   └───shade
└───results
└───inputs
```

`inputs` store includes domain specification metadata including what connectivity/DHW/wave data was used.

# Inputs
- `domain` : ADRIA scenario domain
- `param_df` : ADRIA scenario specification

# Notes
- `domain` is updated with scenario invoke time.
- -9999.0 is used as an arbitrary fill value.

"""
function setup_result_store!(domain::Domain, param_df::DataFrame, reps::Int)
    @set! domain.scenario_invoke_time = replace(string(now()), "T"=>"_", ":"=>"_", "."=>"_")
    log_location = joinpath(ENV["OUTPUT_DIR"], "$(domain.name)__$(domain.scenario_invoke_time)")

    z_store = DirectoryStore(log_location)

    # Store copy of inputs
    input_loc = joinpath(z_store.folder, INPUTS)
    input_dims = (nrow(param_df), ncol(param_df))
    attrs = Dict(
        :name => domain.name,
        :columns => names(param_df),
        :invoke_time => domain.scenario_invoke_time,
        :ADRIA_VERSION => PkgVersion.Version(@__MODULE__),
        
        :site_data_file => domain.env_layer_md.site_data_fn,
        :site_id_col => domain.env_layer_md.site_id_col,
        :unique_site_id_col => domain.env_layer_md.unique_site_id_col,
        :init_coral_cover_file => domain.env_layer_md.init_coral_cov_fn,
        :connectivity_file => domain.env_layer_md.connectivity_fn,
        :DHW_file => domain.env_layer_md.DHW_fn,
        :wave_file => domain.env_layer_md.wave_fn,
        :sim_constants => Dict(fn=>getfield(domain.sim_constants, fn) for fn ∈ fieldnames(typeof(domain.sim_constants)))
    )

    inputs = zcreate(Float64, input_dims..., fill_value=-9999.0, path=input_loc, chunks=input_dims, attrs=attrs)
    inputs[:, :] = Matrix(param_df)

    # Raw result store
    result_loc = joinpath(z_store.folder, RESULTS)
    tf, n_sites = domain.sim_constants.tf, domain.coral_growth.n_sites
    result_dims = (tf, domain.coral_growth.n_species, n_sites, reps, nrow(param_df))
    
    attrs = Dict(
        :structure => ("timesteps", "species", "sites", "reps", "scenarios"),
        :unique_site_ids => domain.unique_site_ids,
    )
    raw = zcreate(Float32, result_dims..., fill_value=-9999.0, path=result_loc, chunks=(result_dims[1:end-1]..., 1), attrs=attrs)

    # Set up logs for site ranks, seed/fog log
    zgroup(z_store, LOG_GRP)
    log_fn = joinpath(z_store.folder, LOG_GRP)

    # Store ranked sites
    n_scens = nrow(param_df)
    rank_dims = (n_sites, 2, reps, n_scens)  # sites, site id and rank, reps, no. scenarios
    fog_dims = (tf, n_sites, reps, n_scens)  # timeframe, sites, reps, no. scenarios

    # tf, no. intervention sites, site id and rank, no. scenarios
    seed_dims = (tf, 2, n_sites, reps, n_scens)

    attrs = Dict(
        :structure=> ("sites", "seed/fog/shade", "reps", "scenarios"),
        :unique_site_ids=>domain.unique_site_ids,
    )
    ranks = zcreate(Float32, rank_dims..., name="rankings", fill_value=-9999.0, path=log_fn, chunks=(rank_dims[1:3]..., 1), attrs=attrs)

    attrs = Dict(
        :structure=> ("timesteps", "intervened sites", "coral type", "scenarios"),
        :unique_site_ids=>domain.unique_site_ids,
    )
    seed_log = zcreate(Float32, seed_dims..., name="seed", fill_value=-9999.0, path=log_fn, chunks=(seed_dims[1:4]..., 1))
    fog_log = zcreate(Float32, fog_dims..., name="fog", fill_value=-9999.0, path=log_fn, chunks=(fog_dims[1:3]..., 1))
    shade_log = zcreate(Float32, fog_dims..., name="shade", fill_value=-9999.0, path=log_fn, chunks=(fog_dims[1:3]..., 1))

    return raw, ranks, seed_log, fog_log, shade_log
end


function load_results(result_loc)
    result_set = zopen(joinpath(result_loc, RESULTS))
    log_set = zopen(joinpath(result_loc, LOG_GRP))
    input_set = zopen(joinpath(result_loc, INPUTS))

    input_cols = input_set.attrs["columns"]
    inputs_used = DataFrame(input_set[:, :], input_cols)

    env_layer_md = EnvLayer(
        input_set.attrs["site_data_file"],
        input_set.attrs["site_id_col"],
        input_set.attrs["unique_site_id_col"],
        input_set.attrs["init_coral_cover_file"],
        input_set.attrs["connectivity_file"],
        input_set.attrs["DHW_file"],
        input_set.attrs["wave_file"],
    )

    return ResultSet(input_set.attrs["name"], 
                     input_set.attrs["invoke_time"],
                     env_layer_md,
                     inputs_used,
                     input_set.attrs["sim_constants"],
                     result_set,
                     log_set["rankings"],
                     log_set["seed"],
                     log_set["fog"],
                     log_set["shade"])

end
