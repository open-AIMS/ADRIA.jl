# Baseline regression harness (opt-in)

Developer tool for refactor stages: shows whether a change altered simulation results.
The harness itself is **not** part of the test suite and adds no dependencies (the runner
uses the standard library only; the comparison uses Zarr/JSON from the ADRIA environment).
Nothing is written inside the repository. The runner is a plain Julia script, so the same
commands work on Linux, macOS and Windows, locally and in CI. Only the logic of the
comparison and the runner (exit codes, `META` rules, argument handling) is tested by the
suite, on synthetic data (`ADRIA/test/regression_harness.jl`).

## Use

```sh
# 1. Before the refactor: build and cache the reference from a rev (default HEAD).
#    Runs the driver in a throwaway git worktree, which is removed afterwards.
julia ADRIA/test/regression/run_regression.jl ref <rev>

# 2. After changes: run the driver on a checkout and compare against the reference.
julia ADRIA/test/regression/run_regression.jl compare [checkout] [rev]
julia ADRIA/test/regression/run_regression.jl compare --tol 1e-4     # fallback, see below
#    (CHECKOUT and REV are optional and positional; --tol may come before or after them)

# 3. Check the comparison itself.
julia ADRIA/test/regression/run_regression.jl selftest
```

The runner needs no project environment of its own (it uses the standard library only), so
`julia ADRIA/test/regression/run_regression.jl ...` works from any directory. It runs under
the Julia that launched it, and starts the driver and the comparison with that same Julia:
to use another version, launch the runner with it (`julia +1.12 ADRIA/test/regression/run_regression.jl ...`
with juliaup, or the path of a Julia binary). Use the same Julia for `ref` and `compare`.

`compare` and `selftest` use the rev recorded by the last `ref` (the cache's `CURRENT`)
unless `compare` is given a REV, which is resolved in this repository like the one of
`ref`; a missing `CURRENT`, reference folder or `META` is reported as an error (exit 2).

The runner exits with 0 (identical, or within `--tol`), 1 (differences), 2 (usage error, or
a refusal: an environment mismatch, missing run metadata or no reference, see below) or 3
(a crash: the driver, `compare.jl` or the runner itself failed). `selftest` checks the exact
exit code of each case, so a crash cannot pass as an expected difference; it covers value,
array, attr, `model_spec.csv`, case-set and `attached`-case perturbations, the always-exact
arrays under `--tol`, and the `META` rules. It runs every case through
`run_regression.jl compare-dirs REF_DIR CAND_DIR [--tol X] [--driver-variant]`, which
compares two existing run folders without running the driver.

Environment: `ADRIA_REGRESSION_REF` (cache, default `~/.cache/adria-regression`),
`ADRIA_REGRESSION_DRIVER` (driver script, default `driver.jl`), `ADRIA_REGRESSION_THREADS`
(default 4; use the same value for reference and candidate), `TMPDIR` (`TEMP` on Windows;
scratch folders, keep it short on Windows because of the path length limit).

Files: `run_regression.jl` (the CLI), `compare.jl` (the comparison, run in the checkout's
ADRIA environment), `meta.jl` (the `META` helpers both use), `driver.jl` (the pinned
workload) and `selftest.jl`.

## Julia version and run metadata

Exact equality holds only for the same Julia version: the last digits of results can
differ between versions (for example a fixed-seed cover total is 163.98244844000132 on
Julia 1.11.9 and 1.12.7, 163.9824484400013 on 1.13.4). The refactor gate is therefore
"same Julia, exact"; the reference used so far was built with Julia 1.12.7, but any
version works as long as the reference and the candidate use the same one.

Every run directory (reference and candidate) gets a `META` file: `julia_version` and
`threads` are written by `driver.jl` itself (`VERSION` and `Threads.nthreads()` of the
process that produced the results), and `run_regression.jl` appends `driver`,
`driver_sha256`, `source` (rev or checkout), `os` (`Sys.KERNEL`) and `arch` (`Sys.ARCH`).
`driver_sha256` is the sha256 of the driver with CRLF line ends normalised to LF, so
`core.autocrlf` on Windows does not change it (the recorded value for `driver.jl` is
`d48b2075d6e1...`). The runner fails if any value is empty, and `compare.jl` treats a
missing or empty value like a missing `META` (two empty values never match).
`compare.jl` prints both and then:

- different Julia version, thread count or driver: fails with exit code 2 in exact mode
  (rebuild the reference with the candidate's Julia and `ADRIA_REGRESSION_THREADS`, or use
  the driver it was built with); with `--tol X` it warns and compares under the
  tolerance, which is the way to compare across Julia versions;
- `--driver-variant` accepts a different driver in exact mode (a per-stage driver variant,
  see below); thread count and Julia version must still match;
- different OS or CPU architecture: always fails with exit code 2, also under `--tol`: a
  reference built on one OS must not be accepted on another (build one per OS);
- no `META`, or a `META` with a missing or empty value: fails with exit code 2 in exact
  mode; a warning under `--tol X`.

Keys and values are compared with surrounding whitespace trimmed, by `compare.jl` and by
the runner alike.

A cached reference is keyed on its rev (the folder `<cache>/<rev>`) and validated by its
`META`: `ref` refuses to reuse one that has no `META`, a missing or empty value, or whose
`META` does not match the current Julia version, thread count, driver checksum, OS and
architecture (exit code 2); remove it to rebuild, or point `ADRIA_REGRESSION_REF`
elsewhere.

**Existing caches must be rebuilt after this change**: a reference cached by an earlier
version of the harness has no `META`, or a `META` without `os` and `arch`, and is refused.
Remove the cached `<rev>` directory and run `ref <rev>` again.

## Windows

The runner passes commands to child processes as argument lists (never through a shell),
uses no `sha256sum`/`mktemp`/`/`-paths, waits for the driver to exit before removing the
worktree (falling back to deleting the folder and `git worktree prune`), keeps scratch
paths short (`<temp>\adr-XXXXXX\w`) and enables `core.longpaths` for the worktree. Run it
from a shell that has `git` and `julia` on `PATH`.

## What is compared

`driver.jl` runs a fixed-seed ensemble on `Test_domain` (RCP 45): 8 guided scenarios with
COTS off, the same 8 with COTS on (juvenile predation on for half), and 4 unguided
scenarios with COTS off and on (`mixed` case, 24 scenarios), the 8 guided COTS-off
scenarios alone (`cots_off` case, a store without COTS logs), and 4 guided COTS-off plus 4
COTS-on scenarios on a domain with the `ShrinkMock(0.1)` test submodel attached
(`attached` case: the submodel hook path, the per-cell decrease-only check and the
propagation of the reduced cover). The submodel's own logs are not written to the store;
the case checks its effect on the outcomes.

`ShrinkMock` is defined in `driver.jl` itself (a pinned copy of the test mock in
`ADRIAinterfaces/test/MockSubmodels.jl`), so the reference and the candidate run identical
mock code even if a stage changes the checkout's mocks. The `attached` case needs
`_attach_submodel!` and the ADRIAinterfaces hook API in the checkout, so a reference with
it can only be built from revs that have the submodel hook plumbing (331795c88 or later).

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
`ADRIA_REGRESSION_DRIVER` for the compare run, passing `--driver-variant`
(`compare --driver-variant`), since a driver change is otherwise refused in exact mode. The reference is always built with the
driver that existed when it was built.
