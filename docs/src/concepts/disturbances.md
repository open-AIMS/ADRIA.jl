# Disturbances

## Cyclones

Each cyclone mortality scenario is defined as a series of cyclone mortality rates for each timestep, location and functional group. They are the result of applying a set of cyclone stochastic generated category projections (Bozec et al., 2025) for each location and timestep, converted to windspeed (Bureau of Meteorology, 2025), to a regression that provides a coral mortality rate as a function of wind speed. The cyclone categories go from 0 (no cyclone) to 5 (maximum wind speed cyclone). The three regression models, for massive corals, branching corals deeper than 5 meters and branching corals shallower than 5 meters, were adjusted for a dataset extract from Fabricius et al. (2008). When the model is run, a cyclone mortality scenario is used, meaning that at each timestep, a mortality rate is applied to each location and functional group.

## Heat Stress

[To do]

## References

1. Bozec, Y. M., Adam, A. A., Nava, B. A., Cresswell, A. K., Haller-Bull, V., Iwanaga, T., ... & Mumby, P. J. (2025). A rapidly closing window for coral persistence under global warming. bioRxiv, 2025-01.

2. Bureau of Meteorology. (2025). Tropical cyclone categories. Australian Government. http://www.bom.gov.au/cyclone/tropical-cyclone-knowledge-centre/understanding/categories/

3. Fabricius, K. E., De'Ath, G., Puotinen, M. L., Done, T., Cooper, T. F., & Burgess, S. C. (2008). Disturbance gradients on inshore and offshore coral reefs caused by a severe tropical cyclone. Limnology and Oceanography, 53(2), 690-704.