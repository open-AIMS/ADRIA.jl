# Testing the COTS integration

How to check that a change to the crown-of-thorns starfish (COTS) integration hasn't broken
ADRIA or changed its default behaviour. The integration comprises a low-DHW bleaching
reduction, the COTS submodel and optional COTS inputs in the data package.

## Which check to run

Run the quick test for the part you changed while you work. Before pushing, run the full
suite, the formatter and the comparison with `main`. The times below assume ADRIA has
already compiled; the first `using ADRIA` in a new checkout takes 5 to 10 minutes.

| You changed | Run | Time |
| --- | --- | --- |
| COTS factors or the adapter (`ecosystem/cots/cots_factors.jl`, `ecosystem/cots/cots_dynamics.jl`) | `test/ecosystem/cots.jl` | about 30 seconds |
| Data-package loading (`ecosystem/cots/cots_inputs.jl`, `Domain.jl`) | `test/io/cots_inputs.jl` | about 1 minute |
| Calibration-file loading or the low-DHW bleaching reduction (`calib_params.jl`, `thermal_tolerance.jl`) | `test/calib_params_loading.jl` | about 4 minutes |
| The COTS step in `run_model` (`model_phases.jl`) or the COTS result logs | `test/ecosystem/cots_scenarios.jl` | about 3 minutes |
| Anything that could change COTS results (on or off) | `test/ecosystem/cots_regression.jl` | about 40 seconds |
| Anything that could change results with COTS off | the comparison with `main` | about 10 minutes |
| Anything, before pushing | the full test suite and the formatter | 20 to 30 minutes |

## Quick checks

Each COTS test file runs on its own. Run these from the `ADRIA/` folder:

```bash
julia --project=. -e 'include("test/ecosystem/cots.jl")'
julia --project=. -e 'include("test/io/cots_inputs.jl")'
julia --project=. -e 'using Test, ADRIA; include("test/calib_params_loading.jl")'
julia --project=. -e 'include("test/ecosystem/cots_scenarios.jl")'
julia --project=. -e 'include("test/ecosystem/cots_regression.jl")'
```

`cots_regression.jl` is a fixed-seed regression test. It runs paired COTS-off and COTS-on
scenarios on Test_domain (derived and fixed `ADRIA_RNG_SEED` seeds) and on a copy of it
with every COTS input in the data package (all annual-forcing modes, derived and fixed
larval seeds). It checks the sums and extrema of relative cover and the COTS logs against
values recorded before the COTS code was restructured (`rtol=1e-12`), and that COTS lowers
coral cover. Its results do not depend on the number of threads. To record new reference
values after a deliberate change in results, run it with `COTS_REGRESSION_PRINT=true`.

Each prints a `Test Summary` table. Any `Fail` or `Error` column means something broke.
When this guide was last updated (November 2026) the files ran 34, 37, 60, 52 and 27
checks, in the order listed above.

The scenario tests write their result stores to temporary folders (see
`with_temp_output_dir` in `test/output_sandbox.jl`), which are removed when Julia exits, so
nothing is left in the current folder. If an `Outputs/` folder appears after a test run, a
test is missing that wrapper.

## Full test suite and formatting

The full suite is what CI's `test` check runs. From the `ADRIA/` folder:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

It runs five tiers in order and stops at the first tier that fails. Results when this guide
was last updated:

| Tier | Passed |
| --- | --- |
| 1: Pre-Smoke | 185 |
| 2: I/O Smokescreen | 254 |
| 3: Spec & Domain Loads | 3670, plus 2 "Broken" |
| 4: Simulation Gate | 10 |
| 5: Simulation-Dependent | 961 |

The 2 "Broken" results in Tier 3 are deliberate skips (`@test_skip` in `decisions/mcda.jl`
and `sampling.jl`), not failures.

CI's `code-style` check runs JuliaFormatter 2.5.1 with the repository's
`.JuliaFormatter.toml`. To run it locally, from the repository root:

```bash
julia -e 'using Pkg; Pkg.activate(temp=true); Pkg.add(name="JuliaFormatter", version="2.5.1"); using JuliaFormatter; format(".")'
git status
```

Any file it changes would fail the check, so commit the formatted version. On Windows, git
may also warn about line endings; those warnings are harmless.

## Is COTS-off still identical to `main`?

With every COTS factor at its default, a branch must give exactly the same results as
`main`. Two short scripts check this. The first runs 32 sampled Test_domain scenarios with a
fixed sampling seed; the second compares two runs value by value.

1. Put a checkout of `main` next to the branch. From the repository root:
   `git worktree add ../ADRIA-main origin/main`. If `main` has moved on, merge it into the
   branch first, so both sides start from the same `main`.
2. Save the two scripts below in a folder outside the repository, for example
   `../cots_check/`.
3. From the repository root, run `main`, then the branch, then the comparison:

```bash
julia --project=../ADRIA-main/ADRIA ../cots_check/cots_legacy_run.jl ../cots_check/out_main
julia --project=ADRIA ../cots_check/cots_legacy_run.jl ../cots_check/out_branch
julia --project=ADRIA ../cots_check/cots_legacy_compare.jl ../cots_check/out_main ../cots_check/out_branch
```

A pass ends with `RESULT: identical`, after seven `same` lines, one per outcome. The first
line counts the columns the branch adds to the scenario table; none of them may vary.

| What the comparison reports | Likely cause |
| --- | --- |
| Shared scenario columns differ | `sample()` changed, for example because a new factor isn't constant by default |
| A new column varies | A new factor is being sampled; collapse its sampling bounds to its value |
| Tables match, but outcomes differ | The COTS-off path of `run_model` changed, or the per-scenario random seed did |

`cots_legacy_run.jl`:

```julia
# Run 32 sampled Test_domain scenarios and save the scenario table and every outcome.
#   julia --project=<ADRIA folder> cots_legacy_run.jl <output folder>
using ADRIA, CSV, DataFrames, Random

out = mkpath(ARGS[1])
ADRIA.setup()  # sets the output folder to ./Outputs, so override it afterwards
ENV["ADRIA_OUTPUT_DIR"] = mktempdir()  # keep result stores out of the repository
dom = ADRIA.load_domain(joinpath(pkgdir(ADRIA), "test", "data", "Test_domain"), "45")
Random.seed!(20261022)  # sample() draws from the global RNG
scens = ADRIA.sample(dom, 32)
CSV.write(joinpath(out, "scenarios.csv"), scens)

rs = ADRIA.run_scenarios(dom, scens, "45")
rows = DataFrame(; outcome=String[], index=Int[], value=Float64[])
for (name, cube) in pairs(rs.outcomes)
    v = vec(Float64.(Array(cube)))
    append!(rows, DataFrame(; outcome=string(name), index=eachindex(v), value=v))
end
CSV.write(joinpath(out, "outcomes.csv"), rows)
println("ADRIA at ", pkgdir(ADRIA), ": ", nrow(scens), " scenarios written to ", out)
```

`cots_legacy_compare.jl`:

```julia
# Compare two cots_legacy_run.jl outputs (main first). Exits 1 if anything differs.
#   julia --project=<ADRIA folder> cots_legacy_compare.jl <main output> <branch output>
using CSV, DataFrames

main_dir, branch_dir = ARGS
sa = CSV.read(joinpath(main_dir, "scenarios.csv"), DataFrame)
sb = CSV.read(joinpath(branch_dir, "scenarios.csv"), DataFrame)
shared = intersect(names(sa), names(sb))
differ = [c for c in shared if !isequal(sa[!, c], sb[!, c])]
varying = [c for c in setdiff(names(sb), names(sa)) if length(unique(sb[!, c])) > 1]
println("Scenario table: $(length(shared)) shared columns, $(length(differ)) differ; ",
    "$(length(names(sb)) - length(shared)) new columns, $(length(varying)) of them vary")
isempty(differ) || println("  differ: ", join(first(differ, 10), ", "), length(differ) > 10 ? ", ..." : "")
isempty(varying) || println("  vary: ", join(first(varying, 10), ", "), length(varying) > 10 ? ", ..." : "")

oa = CSV.read(joinpath(main_dir, "outcomes.csv"), DataFrame)
ob = CSV.read(joinpath(branch_dir, "outcomes.csv"), DataFrame)
ok = isempty(differ) && isempty(varying)
for name in sort(unique(vcat(oa.outcome, ob.outcome)))
    x = sort(oa[oa.outcome .== name, :], :index).value
    y = sort(ob[ob.outcome .== name, :], :index).value
    same = isequal(x, y)  # exact, with NaN equal to NaN
    global ok &= same
    println(same ? "same       " : "DIFFERENT  ", name)
end
println(ok ? "RESULT: identical" : "RESULT: differences found")
exit(ok ? 0 : 1)
```

## Adding a new integration point

These rules keep the comparison with `main` identical and keep new settings calibratable:

- **Add each setting as a factor.** Add a `ParamSpec` to `COTS_PARAM_SPECS` (from which `COTS_FACTOR_SPECS` is built) and a field to
  `CotsParams` in `ecosystem/cots/cots_factors.jl`. `_cots_factor` then makes it constant by
  default, keeps its full bounds for `set_factor_bounds!`, excludes it from the derived
  per-scenario seed, and lets the `cots` variable of a calibration NetCDF set it.
- **Give it a default that changes nothing:** off (0), a multiplier of 1, or a rate of 0.
- **Read settings from `param_set`,** never from environment variables.
- **Keep COTS code inside `if cots_on` in `run_model`.** Don't draw from the scenario's
  `rng` outside that block, or coral-only scenarios will change.
- **Make new data optional.** Add it as a field of `CotsInputs` and load it in
  `load_cots_inputs` (`ecosystem/cots/cots_inputs.jl`). A missing file gives `nothing`, not an error.
- **Give any new domain type a `CotsParams` component.** `run_model` reads `cots_on` in
  every scenario. `ADRIADomain`, `RMEDomain` and `ReefModDomain` get it from
  `load_calib_params`.
- **A new factor outside COTS,** like the low-DHW bleaching reduction (`LowDHWBleaching`),
  must collapse its sampling bounds to its value and default to 0. Otherwise `sample()` and
  the derived seed change.
- **Add a test** next to the existing ones: `test/ecosystem/cots.jl` for the adapter and
  factors, `test/io/cots_inputs.jl` for data loading, and `test/ecosystem/cots_scenarios.jl`
  for a change to the COTS step.
