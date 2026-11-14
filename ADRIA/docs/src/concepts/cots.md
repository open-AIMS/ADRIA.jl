# Crown-of-thorns starfish (COTS)

ADRIA can run a crown-of-thorns starfish (COTS, *Acanthaster* spp.) population at every
location, which grazes coral each year. The population model lives in the
[COTSMod.jl](https://github.com/open-AIMS/COTSMod.jl) package. ADRIA passes it the coral
cover, takes back the coral it ate, and logs the COTS state.

The submodel is **off by default**. With it off, ADRIA runs the coral-only model, and
results are identical to versions of ADRIA without COTS.

## Switching it on

COTS is controlled by the `cots_on` factor (0 = off, 1 = on). Like the coral components, all
COTS factors are fixed by default in `sample()`. Set them in the scenario table, or fix
them on the domain before sampling:

```julia
dom = ADRIA.load_domain("path/to/domain", "45")
ADRIA.fix_factor!(dom; cots_on=1)
scens = ADRIA.sample(dom, 16)
rs = ADRIA.run_scenarios(dom, scens, "45")
```

## What happens each year

COTS act after the year's bleaching and cyclone mortality:

1. Each location's COTS population advances one year. It recruits, matures, feeds and
   dies, with feeding and starvation set by the cover of fast-growing coral (the three
   *Acropora* groups) and slow-growing coral (the two massive groups).
2. The coral eaten is removed from that location's cover.
3. COTS larvae disperse between locations. By default they use the domain's connectivity:
   a data package may supply a dedicated COTS connectivity matrix in `cots_connectivity/`,
   otherwise the coral connectivity is used. With the optional reef-scale inputs (see
   [Optional data-package inputs](#optional-data-package-inputs)) larvae can instead
   disperse through annual reef-scale matrices (`cots_larval_conn_mode`), settle only
   through that dispersal (`cots_separated_recruitment`), and be supplemented by
   outbreaks outside the domain (`cots_external_outbreak_production`).

Populations start from COTSMod's default (a random 25% of locations hold 0.1 adults per
hectare) unless the data package supplies starting densities. The default draw uses the
scenario's own random number generator, so it is reproducible.

## Factors

The `CotsParams` component holds every COTS setting as a model factor, so each one is
recorded in the scenario table and the result set.

| Group | Factors |
| --- | --- |
| switch | `cots_on` |
| recruitment | `cots_bh_a`, `cots_bh_b`, `cots_ricker_a`, `cots_ricker_b`, `cots_allee_threshold`, `cots_fecundity_gate` |
| mortality and starvation | `cots_mortality_recruit`, `cots_mortality_juvenile`, `cots_mortality_adult`, `cots_starvation_max_fraction`, `cots_max_food_cover`, `cots_starvation_exponent`, `cots_condition_timescale` |
| feeding | `cots_consumption_rate_fast`, `cots_consumption_rate_slow`, `cots_handling_time`, `cots_consumption_exponent_fast`, `cots_consumption_exponent_slow` |
| immigration | `cots_immigration_base`, `cots_immigration_threshold`, `cots_immigration_exponent` |
| optional mechanisms (0 = off) | `cots_adult_senescence`, `cots_allee_stretched`, `cots_juvenile_predation`, `cots_per_capita_consumption`, `cots_preferred_condition`, `cots_preferred_prey_starvation` |
| reef-scale dispersal and inputs | `cots_separated_recruitment`, `cots_larval_fecundity`, `cots_larval_conn_mode`, `cots_larval_conn_seed`, `cots_external_outbreak_production`, `cots_initial_state_scale` |
| mechanism parameters | `cots_senescence_age`, `cots_senescence_mortality`, `cots_allee_exponent`, `cots_juvenile_predation_capacity`, `cots_consumption_m2_per_adult_year`, `cots_condition_saturation`, `cots_condition_memory_years`, `cots_starvation_preferred_threshold`, `cots_starvation_memory_years` |

Every factor name starts with `cots_`. The factors that are COTSMod parameters are mapped
to COTSMod's own parameter names by `COTS_PARAM_MAP` (for example `cots_mortality_recruit`
is COTSMod's `m1`).

Densities and the Allee threshold are COTS per hectare. The factor bounds are broad
placeholders for sampling, not evidence-derived ranges.

### Calibrated values

Calibrated values can be read from the same calibration NetCDF as the coral parameters
(`calib_params_fn` in `load_domain`). Store them in a variable named `cots` with one
dimension, `cots_param`, whose labels are the factor names. Values are matched by label,
and factors missing from the file keep their defaults.

### Random seeds

When `ADRIA_RNG_SEED` is not set, ADRIA derives each scenario's seed from the sum of its
parameter values. The COTS factors are left out of that sum, so adding the COTS component
does not change the seed, or the results, of any coral-only scenario. The one exception is
`cots_larval_conn_seed`, which seeds its own generator for the annual-matrix draws.

## Outputs

When at least one scenario has `cots_on = 1`, the result set holds two extra logs:

- `cots_pop_log`: COTS density by age class (recruits, juveniles, adults), with dimensions
  timesteps x age class x locations x scenarios, in COTS per hectare;
- `cots_condition_log`: COTS body condition, with dimensions timesteps x locations x
  scenarios.

Scenarios with COTS off hold zeros in these logs. When no scenario runs COTS, the logs are
not created and both fields are `nothing`.

## Optional data-package inputs

A data package can carry extra COTS inputs. All are optional; without them COTS runs as
described above. Their layout and formats are in the `CotsInputs` docstring.

- `cots_connectivity/annual_forcing/`: reef-scale larval connectivity, one matrix per
  year, plus larval supply from outbreaks outside the domain. Needed by separated
  recruitment, the annual connectivity modes and external outbreak supply.
- `cots_initial_state.csv`: starting densities per location (scaled by
  `cots_initial_state_scale`).
- `cots_larval_survival.csv`: larval survival multipliers per location and year.
- `cots_settlement_modifier.csv`: settlement multipliers per location.

`cots_larval_conn_mode` selects the connectivity: 0 uses the mean matrix, 1 steps through
the annual matrices in turn, and 2 draws one each year using `cots_larval_conn_seed`.
Unsupported combinations (for example separated recruitment without annual connectivity)
are rejected before any scenario runs.

## Not yet included

COTS control (culling) is not yet part of ADRIA. It is being developed and calibrated in
[COTSMod-Calibration](https://github.com/open-AIMS/COTSMod-Calibration) for a fine-scale
Lizard Island domain.
