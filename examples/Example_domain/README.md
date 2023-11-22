Overview of data included in this data package.

This is a release candidate and not intended for use in production.

**site_data**

Spatial data relating to the Moore cluster (from A. Cresswell).
Zone type data attached separately by T. Iwanaga, using zone data extracted by M. Puotinen

Consists of a geopackage file.

- Moore_2023-09-14.gpkg

Pertinent columns are:
- reef_siteid, (Note: unique site id, not to be confused with other "site_ids")
- area (in m^2)
- k, maximum carrying capacity for the associated site (in percentage values, gets divided by 100 in ADRIA)
- depth_med, median depth across site (meters)
- depth_mean, mean depth across site (meters)
- zone_type, management zone colors based on:
  - https://www2.gbrmpa.gov.au/access/zoning/interpreting-zones
  - https://elibrary.gbrmpa.gov.au/jspui/retrieve/f8cdc701-ef8f-41fc-b96f-12d7d55522d5/Map-5-Cairns.pdf

**connectivity**

Updated connectivity matrices from RECOM runs, provided by C. Ani.

**cyclone**

Provide coral mortality projections due to cyclones. Added in version 0.3.2-rc by P. R. Almeida.

**DHWs**

Degree Heating Weeks datasets v3.0 by C. Ani and V. Lago.

**Wave data**

Wave stress data provided by B. Robson.


**initial coral cover**

Downscaled initial coral cover values (0 to 1) that are relative to k area.
Initial coral cover is downscaled from ReefModEngine data (v1.0.18 / rme_ml_2023_03_30b).

RME data loaded by ADRIA distributes cover over size classes weighted by the PDF densities of ReefMod's size distribution (log-normal with mean log(700cm^2) and std log(4cm^2)).
This initial coral cover data is then distributed over sites according to the indicated available k area for each site, and exported in netCDF format.
