# Disturbances

## Heat Stress

Heat stress in ADRIA is represented using **Degree Heating Weeks** (DHW), a standard
oceanographic metric that accumulates positive thermal anomalies above a bleaching
threshold over time. DHW data is supplied as part of the domain package as a
three-dimensional array (timesteps x locations x scenarios). The `dhw_scenario`
parameter in the scenario table selects which projection is applied in a given run.

## Bleaching mortality

Bleaching mortality is computed using a population-level probabilistic model. Each
combination of functional group, size class, and location has a **critical DHW tolerance
distribution** -- a truncated normal distribution describing the spread of thermal
tolerances within that sub-population. The lower bound of the distribution is fixed at
`HEAT_LB = 1.0` DHW-week; the mean cannot increase more than `HEAT_UB = 24.0` DHW-weeks
above its initial value.

At each timestep, locations where DHW exceeds `HEAT_LB` are considered thermally active.
For each active location the proportion of the sub-population that bleaches is given by
the cumulative distribution function (CDF) of the tolerance distribution evaluated at the
current DHW value. Bleaching mortality is then depth-adjusted:

```
mortality = bleaching_fraction * depth_coefficient(depth)
```

The depth coefficient follows Baird et al. (2018):

```
depth_coeff = exp(-0.07551 * (depth - 2.0))
```

clamped to [0, 1], so that deeper reefs experience proportionally less bleaching mortality.

## Functional group and size class differences

Each functional group and size class carries its own mean critical DHW threshold and
standard deviation, initialised from observational data (Hughes et al. 2018,
Bairos-Novak et al. 2021). Bleaching sensitivity therefore differs between functional
groups, with tabular Acropora being more sensitive than, for example, massive corals.

## Thermal adaptation

Population mean tolerance shifts each timestep through two mechanisms:

- **Survival selection**: After a bleaching event, individuals with higher tolerances are
  preferentially retained, shifting the surviving population mean upward. This is
  formalised via the Breeder's equation using a fixed heritability parameter.
- **Settler tolerance**: Newly seeded corals can carry a user-specified thermal tolerance
  offset (`a_adapt`), representing a given level of enhanced thermal tolerance, through
  assisted gene flow, adaptation or other enhancement process.

The mean tolerance is hard-capped at `initial_mean + HEAT_UB` to represent a biological
ceiling on adaptation.

## Fogging and Solar Radiation Management (SRM)

Interventions that reduce light and heat reaching corals (fogging and SRM) are modelled
as a direct multiplicative reduction of the DHW experienced at selected locations:

```
effective_DHW = DHW * (1 - fogging_effectiveness)
```

This reduction is applied before bleaching mortality is computed, so treated locations
experience reduced bleaching proportional to the intervention intensity.

## Cyclones

Each cyclone mortality scenario is defined as a series of cyclone mortality rates for each timestep, location and functional group. They are the result of applying a set of cyclone stochastic generated category projections (Bozec et al., 2025) for each location and timestep, converted to windspeed (Bureau of Meteorology, 2025), to a regression that provides a coral mortality rate as a function of wind speed. The cyclone categories go from 0 (no cyclone) to 5 (maximum wind speed cyclone). The three regression models, for massive corals, branching corals deeper than 5 meters and branching corals shallower than 5 meters, were adjusted for a dataset extract from Fabricius et al. (2008). When the model is run, a cyclone mortality scenario is used, meaning that at each timestep, a mortality rate is applied to each location and functional group.

## References

1. Baird, A., Madin, J., Alvarez-Noriega, M., Fontoura, L., Kerry, J., Kuo, C., Precoda, K., Torres-Pulliza, D., Woods, R., Zawada, K., & Hughes, T. (2018). A decline in bleaching suggests that depth can provide a refuge from global warming in most coral taxa. Marine Ecology Progress Series, 603, 257-264. https://doi.org/10.3354/meps12732

2. Bairos-Novak, K. R., Hoogenboom, M. O., van Oppen, M. J., & Connolly, S. R. (2021). Coral adaptation to climate change: Meta-analysis reveals high heritability across multiple traits. Global Change Biology, 27, 5694-5710. https://doi.org/10.1111/gcb.15829

3. Bozec, Y.-M., Hock, K., Mason, R. A. B., Baird, M. E., Castro-Sanguino, C., Condie, S. A., Puotinen, M., Thompson, A., & Mumby, P. J. (2022). Cumulative impacts across Australia's Great Barrier Reef: A mechanistic evaluation. Ecological Monographs, 92(1), e01494. https://doi.org/10.1002/ecm.1494

4. Bozec, Y. M., Adam, A. A., Nava, B. A., Cresswell, A. K., Haller-Bull, V., Iwanaga, T., ... & Mumby, P. J. (2025). A rapidly closing window for coral persistence under global warming. bioRxiv, 2025-01.

5. Bureau of Meteorology. (2025). Tropical cyclone categories. Australian Government. http://www.bom.gov.au/cyclone/tropical-cyclone-knowledge-centre/understanding/categories/

6. Fabricius, K. E., De'Ath, G., Puotinen, M. L., Done, T., Cooper, T. F., & Burgess, S. C. (2008). Disturbance gradients on inshore and offshore coral reefs caused by a severe tropical cyclone. Limnology and Oceanography, 53(2), 690-704.

7. Hughes, T. P., Kerry, J. T., Baird, A. H., Connolly, S. R., Dietzel, A., Eakin, C. M., Heron, S. F., Hoey, A. S., Hoogenboom, M. O., Liu, G., McWilliam, M. J., Pears, R. J., Pratchett, M. S., Skirving, W. J., Stella, J. S., & Torda, G. (2018). Global warming transforms coral reef assemblages. Nature, 556, 492-496. https://doi.org/10.1038/s41586-018-0041-2