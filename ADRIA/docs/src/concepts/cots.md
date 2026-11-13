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
3. COTS larvae disperse between locations through the domain's connectivity. A domain
   may supply a dedicated COTS connectivity matrix (`cots_conn`); otherwise the coral
   connectivity is used.

Populations start from COTSMod's default: a random 25% of locations hold 0.1 adults per
hectare. The draw uses the scenario's own random number generator, so it is
reproducible.

## Factors

The `CotsParams` component holds every COTS setting as a model factor, so each one is
recorded in the scenario table and the result set.

| Group | Factors |
| --- | --- |
| switch | `cots_on` |
| recruitment | `a`, `b`, `a_ricker`, `b_ricker`, `allee_threshold`, `fecundity_gate` |
| mortality and starvation | `m1`, `m2`, `m3`, `p_tilde`, `C_max`, `eta_starve`, `tau_condition` |
| feeding | `a_F`, `a_S`, `h`, `eta_F`, `eta_S` |
| immigration | `IMM`, `imm_threshold`, `eta_imm` |
| optional mechanisms (0 = off) | `adult_senescence`, `allee_stretched`, `juvenile_predation`, `per_capita_consumption`, `preferred_condition`, `preferred_prey_starvation` |
| mechanism parameters | `senescence_age`, `senescence_mortality`, `allee_exponent`, `juvenile_predation_capacity`, `consumption_m2_per_adult_year`, `condition_saturation`, `condition_memory_years`, `starvation_preferred_threshold`, `starvation_memory_years` |

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
does not change the seed, or the results, of any coral-only scenario.

## Outputs

When at least one scenario has `cots_on = 1`, the result set holds two extra logs:

- `cots_pop_log`: COTS density by age class (recruits, juveniles, adults), with dimensions
  timesteps x age class x locations x scenarios, in COTS per hectare;
- `cots_condition_log`: COTS body condition, with dimensions timesteps x locations x
  scenarios.

Scenarios with COTS off hold zeros in these logs. When no scenario runs COTS, the logs are
not created and both fields are `nothing`.

## Not yet included

Reef-scale and annual larval connectivity, larvae arriving from outside the domain,
starting densities supplied by the domain, and COTS control (culling) are not yet part of
ADRIA. They are being developed and calibrated in
[COTSMod-Calibration](https://github.com/open-AIMS/COTSMod-Calibration) for a fine-scale
Lizard Island domain.
