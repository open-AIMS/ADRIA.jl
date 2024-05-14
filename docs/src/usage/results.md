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
In order to reduce the duplication of geospatial and conectivity data, the data directory
and results directory can be supplied seperately to avoid having copies for each resuilt set
analysed.

```julia
rs = ADRIA.load_domain(RMEResultSet, "<path to data dir>", "<path to results dir>")
```
