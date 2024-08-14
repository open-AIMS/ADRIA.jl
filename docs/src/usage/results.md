# Loading Results

## Loading ReefModEngine Results

Results from ReefModEngine.jl can be loaded with the `load_results` function.

```julia
rs = ADRIA.load_results(RMEResultSet, "<path to data dir>")
```

Expected data directory structure:

```bash
data_dir
│
├───con_bin
│       CONNECT_ACRO_2010_11.bin
│       CONNECT_ACRO_2011_12.bin
│       CONNECT_ACRO_2012_13.bin
│       CONNECT_ACRO_2014_15.bin
│       CONNECT_ACRO_2015_16.bin
│       CONNECT_ACRO_2016_17.bin
│
├───id
│       id_list_2023_03_30.csv
│
├───region
│       reefmod_gbr.gpkg
│
└───results
        results.nc
        scenarios.csv
```
In order to reduce the duplication of geospatial and connectivity data, the data directory
and results directory can be supplied separately to avoid having copies for each result set
analysed.

```julia
rs = ADRIA.load_domain(RMEResultSet, "<path to data dir>", "<path to results dir>")
```

## Loading CScape Results

Results from CScape can be loaded with the `load_results` function.

```julia
# Assumes NetCDFs are contained in result subdirectory
rs = ADRIA.load_results(RMEResultSet, "<path to data dir>")

# Retrieves NetCDFs from seperate directory
rs = ADRIA.load_results(RMEResultSet, "<path to data dir>", "<path to result directory>")

# Manually passes a list of files to load as results
rs = ADRIA.load_results(RMEResultSet, "<path to data dir>", ["netcdf_fn1", "netcdf_fn2", ...])
```

The expected directory structure is
```bash
data_dir
│   ScenarioID.csv
│
├───connectivity
│       connectivity.csv
│
├───site_data
│       geospatial_data.gpkg
│
├───initial_cover
│       initial_cover.csv
│
└───results (optional)
        NetCDF_Scn_140001.nc
        NetCDF_Scn_140002.nc
        NetCDF_Scn_140003.nc
        ...
```
