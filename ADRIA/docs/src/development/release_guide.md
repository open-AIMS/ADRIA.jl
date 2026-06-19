# Release Guide

All releases are done on `main`.

Note: version numbers should follow [Semantic Versioning](https://semver.org/).

## Mono-Repo Structure

This repository contains multiple packages:
- `ADRIA/` — main package
- `ADRIAviz/` — visualization package
- `ADRIAanalysis/` — analysis package

Each package has its own `Project.toml` and version. Releases are coordinated per-package.

## Public "Final" Releases

1. Increase the version number (manually) following SemVer in the appropriate `Project.toml` file:
   - For ADRIA: `ADRIA/Project.toml`
   - For ADRIAviz: `ADRIAviz/Project.toml`
   - For ADRIAanalysis: `ADRIAanalysis/Project.toml`
2. Update the environment and run tests locally to ensure all pass:
   ```bash
   cd ADRIA && julia --project=. -e 'using Pkg; Pkg.test()'
   cd ADRIAviz && julia --project=. -e 'using Pkg; Pkg.test()'
   cd ADRIAanalysis && julia --project=. -e 'using Pkg; Pkg.test()'
   ```
   See [Testing](@ref) for more details.
3. Submit a PR with the updated version number(s).
   Wait for approval.
4. Once approved, go to the [releases page](https://github.com/open-AIMS/ADRIA.jl/releases)
   and click the "Draft a new release" button at the top right.
5. Under "Choose a tag" on the left, enter the new version number with a package prefix (e.g., "ADRIA-v0.99.0", "ADRIAviz-v1.0.0", or "ADRIAanalysis-v2.0.0") and then
   select "Create new tag: <TAG-NAME> on publish".
6. At the top-right of the textbox, select the last full release then click the
   "Generate release notes" button (at the top-right of the textbox).
7. Under "What's Changed" add a short description of the major changes.
   Explicitly note any major breaking changes (i.e., anything that results obtained with previous versions of ADRIA incompatible)
   Copy the release notes (for step 8).
   Click "Publish release".
8. Register the updated package by opening a new issue in the [General registry](https://github.com/JuliaRegistries/General) with the title "Register [package] [version number]"
   e.g., `Register ADRIA v1.0`, `Register ADRIAviz v2.1.0`, or `Register ADRIAanalysis v1.5.0`
9. State in the comment: `@JuliaRegistrator register subdir=ADRIA`
   Paste in the generated text from step 6 (an example is shown below)
10. Submit the issue. The JuliaRegistrator bot should take care of the rest.

```
@JuliaRegistrator register subdir=ADRIA

Release notes:

Paste the generated release notes here.
```

For other packages, use:
- `@JuliaRegistrator register subdir=ADRIAviz`
- `@JuliaRegistrator register subdir=ADRIAanalysis`


See Julia Registrator usage notes [here](https://github.com/JuliaComputing/Registrator.jl?installation_id=32448289&setup_action=install#details-for-triggering-juliaregistrator-for-step-2-above) for more details.


!!! note "Issues can block release"
    The JuliaRegistrator bot submits a corresponding Pull Request with the Julia package registry.
    Registration may be blocked for a number of reasons. Keep an eye on the auto-submitted
    Pull Request and resolve any issues reported there. Otherwise the package version will never be
    released.


## Development Release

Development releases provide users with the most recent "working" version and may still have some known bugs.
It provides users a chance to try new features and/or provide feedback before a public release.

Deploying a Development Release follows the same steps as "Public" releases, except:

- Add "-dev.x" to the version number with package prefix.
   e.g., `ADRIA-v1.2.3-dev.1`, `ADRIAviz-v1.2.3-dev.1`, `ADRIAanalysis-v1.2.3-dev.1`, etc.
- Untick "Set as the latest release" and tick the "Set as a pre-release" option.
- Ignore Step 8 through 10; DO NOT trigger the `JuliaRegistrator` bot.


## Release Candidates

Release candidates are releases that are not yet "final" but are close to it. Release candidates provide a "last chance" opportunity
for users to report bugs prior to a "final" release.

Deploying a Release Candidate follows the same steps as "Public" releases, except:

- Add "-rc.x" to the version number with package prefix.
   e.g., `ADRIA-v1.2.3-rc.1`, `ADRIAviz-v1.2.3-rc.2`, `ADRIAanalysis-v1.0.0-rc.1`, etc.
- Untick "Set as the latest release" and tick the "Set as a pre-release" option.
- Ignore Step 8 through 10; DO NOT trigger the `JuliaRegistrator` bot.
