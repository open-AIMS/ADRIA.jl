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
ADRIA/test/regression/run_regression.sh compare --tol 1e-4     # fallback, see below
#    (CHECKOUT and REV are optional and positional; --tol may come before or after them)

# 3. Check the comparison itself.
ADRIA/test/regression/run_regression.sh selftest
```

Environment: `ADRIA_REGRESSION_REF` (cache, default `~/.cache/adria-regression`),
`ADRIA_REGRESSION_DRIVER` (driver script, default `driver.jl`), `ADRIA_REGRESSION_THREADS`
(default 4; use the same value for reference and candidate), `JULIA`, `TMPDIR`.

`JULIA` is the Julia command and is split on spaces, so a juliaup channel can be selected
(`JULIA="julia +1.12"`) as well as a binary (`JULIA=/path/to/julia-1.12.7/bin/julia`, a
path without spaces). Use the same Julia for `ref` and `compare`. The default `julia` is
whatever juliaup (or `PATH`) resolves, which may not be 1.12; to match the 1.12.7 reference
pass `JULIA="julia +1.12"` explicitly.

## Julia version and run metadata

Exact equality holds only for the same Julia version: the last digits of results can
differ between versions (for example a fixed-seed cover total is 163.98244844000132 on
Julia 1.11.9 and 1.12.7, 163.9824484400013 on 1.13.4). The refactor gate is therefore
"same Julia, exact"; the reference used so far was built with Julia 1.12.7, but any
version works as long as the reference and the candidate use the same one.

Every run directory (reference and candidate) gets a `META` file: `julia_version` and
`threads` are written by `driver.jl` itself (`VERSION` and `Threads.nthreads()` of the
process that produced the results), and `run_regression.sh` appends `driver`,
`driver_sha256` (`sha256sum`, or `shasum -a 256`) and `source` (rev or checkout). The
runner fails if any value is empty, and `compare.jl` treats a missing or empty value like
a missing `META` (two empty values never match). `compare.jl` prints both and then:

- different Julia version: fails with exit code 2 in exact mode (rebuild the reference
  with the candidate's Julia); with `--tol X` it warns and compares under the tolerance,
  which is the way to compare across versions;
- different thread count or driver: warning (a per-stage driver variant is expected);
- no `META` (a reference built before it was recorded): fails with exit code 2 in exact
  mode; a warning under `--tol X`.

`ref` refuses to reuse a cached reference that has no `META`, or whose `META` does not
match the current Julia, thread count and driver (exit code 2); remove it to rebuild, or
point `ADRIA_REGRESSION_REF` elsewhere.

**Existing caches must be rebuilt after this change**: a reference cached by an earlier
version of the harness has no `META` (or an older one) and is refused by `ref`. Remove the
cached `<rev>` directory and run `ref <rev>` again.

## What is compared

`driver.jl` runs a fixed-seed ensemble on `Test_domain` (RCP 45): 8 guided scenarios with
COTS off, the same 8 with COTS on (juvenile predation on for half), and 4 unguided
scenarios with COTS off and on (`mixed` case, 24 scenarios), the 8 guided COTS-off
scenarios alone (`cots_off` case, a store without COTS logs), and 4 guided COTS-off plus 4
COTS-on scenarios on a domain with the `ShrinkMock(0.1)` test submodel attached
(`attached` case: the submodel hook path, the per-cell decrease-only check and the
propagation of the reduced cover). The submodel's own logs are not written to the store;
the case checks its effect on the outcomes.

The `attached` case needs `_attach_submodel!` and `ShrinkMock` in the checkout, so a
reference with it can only be built from revs that have the submodel hook plumbing.

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
