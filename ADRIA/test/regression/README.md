# Baseline regression harness (opt-in)

Developer tool for refactor stages: shows whether a change altered simulation results.
It is **not** wired into `runtests.jl` and adds no dependencies (it uses Zarr/JSON from the
ADRIA environment). Nothing is written inside the repository.

## Use

```sh
# 1. Before the refactor: build and cache the reference from a rev (default HEAD).
#    Runs the driver in a throwaway git worktree, which is removed afterwards.
ADRIA/test/regression/run_regression.sh ref <rev>

# 2. After changes: run the driver on a checkout and compare against the reference.
ADRIA/test/regression/run_regression.sh compare [checkout] [rev]
ADRIA/test/regression/run_regression.sh compare . --tol 1e-4   # fallback, see below

# 3. Check the comparison itself.
ADRIA/test/regression/run_regression.sh selftest
```

Environment: `ADRIA_REGRESSION_REF` (cache, default `~/.cache/adria-regression`),
`ADRIA_REGRESSION_DRIVER` (driver script, default `driver.jl`), `ADRIA_REGRESSION_THREADS`
(default 4; use the same value for reference and candidate), `JULIA`, `TMPDIR`.

## What is compared

`driver.jl` runs a fixed-seed ensemble on `Test_domain` (RCP 45): 8 guided scenarios with
COTS off, the same 8 with COTS on (juvenile predation on for half), and 4 unguided
scenarios with COTS off and on (`mixed` case, 24 scenarios), plus the 8 guided COTS-off
scenarios alone (`cots_off` case, a store without COTS logs).

`compare.jl` compares decoded contents, not file bytes: for every Zarr array the path set,
dtype, shape, chunks, fill value and all values (NaN-aware exact equality); attrs against an
allow-list (ignored: `invoke_time`, `ADRIA_VERSION`; `*_file` attrs by basename);
`model_spec/model_spec.csv` as text, ignoring its `#` header comment (version, time). The geopackage is not compared.

Default mode requires everything to be identical. `--tol X` is a documented fallback for a
stage that deliberately re-orders floating-point work: it relaxes only outcome arrays and
non-intervention logs to a max abs difference of `X`. Ranks, intervention logs
(`logs/rankings`, `coral_aquaculture`, `larval_methods`, `shading_log`), the inputs table and `model_spec.csv` always stay exact. Do not use
it for result-changing stages (compare ensemble statistics instead).

## Per-stage driver variant

If a stage renames an accessor the driver uses, copy `driver.jl` to `driver_<stage>.jl`, keep the ensemble (seed, sampling calls) identical, and set
`ADRIA_REGRESSION_DRIVER` for the compare run. The reference is always built with the
driver that existed when it was built.
