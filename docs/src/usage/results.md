# Loading Results

## Loading ReefModEngine Results

Results from ReefModEngine.jl can be loaded with the `load_results` function. There are two
options depending on the location of the results and accompanying data files. 

```julia
rs = ADRIA.load_results(RMEResultSet, "<path to data dir>", "<path to results dir>")
```

Expected data directory structure.

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
└───region
       reefmod_gbr.gpkg
```

The results data directory should contain a `results.nc` NetCDF file and `scenarios.csv`. If
the path to the results directory is not supplied, the results are expected to be
contained in a results subdirectory.

```julia
rs = ADRIA.load_domain(RMEResultSet, "<path to data dir>")
```

Expected data directory structure for defaulted result directory.

```bash
data_dir
├───con_bin
├───id
├───region
└───results
```
