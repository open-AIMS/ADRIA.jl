Example data package for discussion.

Data package structure:

Example_domain                   # Name of study area
|   datapackage.json             # Machine-readable metadata file
|   README.md                    # This file
|
+---connectivity                 # Directory for connectivity data, organized by year
|   \---2000
|           test_conn_data.csv
|
+---DHWs
|       dhwRCP45.mat             # Stochastic DHW data for each RCP
|
+---site_data                    # Directory of site data
|       coral_cover.mat          # Initial coral cover data (possibly stochastic)
|       Example_domain.gpkg      # Spatial data defining study area/reef/sites (the "domain")
|
\---waves                        # Directory for wave stress data
        wave_RCP45.mat           # Stochastic wave stress data for each RCP