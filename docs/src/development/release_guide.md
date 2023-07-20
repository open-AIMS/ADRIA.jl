# Release Guide

All releases are done on `main`.

Note: version numbers should follow [Semantic Versioning](https://semver.org/).

## Public "Final" Releases

1. Increase the version number (manually) following SemVer in the Project.toml file
2. Update environment and run tests locally to ensure all pass.
   [Testing](@ref)
3. Submit a PR with the updated version number.
   Wait for approval.
4. Once approved, go to the [releases page](https://github.com/open-AIMS/ADRIA.jl/releases)  
   and click the "Draft a new release" button at the top right
5. Under "Choose a tag" on the left, enter the new version number (e.g., "v0.99.0") and then  
   select "Create new tag: <TAG-NAME> on publish"
6. At the top-right of the textbox, select the last full release then click the 
   "Generate release notes" button (at the top-right of the textbox).
7. Under "What's Changed" add a short description of the major changes.  
   Explicitly note any major breaking changes (i.e., anything that results obtained with previous versions of ADRIA incompatible)
   Copy the release notes (for step 8).
   Click "Publish release".
8. Register the updated package by opening a new issue with the title "Register [version number]"  
   e.g., `Register v1.0`
9. State in the comment: `@JuliaRegistrator register`
   Paste in the generated text from step 6 (an example is shown below)
10. Submit the issue. The JuliaRegistrator bot should take care of the rest.

```
@JuliaRegistrator register

Release notes:

Paste the generated release notes here.
```


See Julia Registrator usage notes [here](https://github.com/JuliaComputing/Registrator.jl?installation_id=32448289&setup_action=install#details-for-triggering-juliaregistrator-for-step-2-above) for more details.


## Development Release

Development releases provide users with the most recent "working" version of ADRIA and may still have some known bugs.
It provides users a chance to try new features and/or provide feedback before a public release.

Deploying a Development Release follows the same steps as "Public" releases, except:

- Add "-dev.x" to the version number.  
   e.g., v1.2.3-dev.1; v1.2.3-dev.2 for the second development release, etc.
- Untick "Set as the latest release" and tick the "Set as a pre-release" option.
- Ignore Step 8 through 10; DO NOT trigger the `JuliaRegistrator` bot.

## Release Candidates

Release candidates are releases that are not yet "final" but are close to it. Release candidates provide a "last chance" opportunity
for users to report bugs prior to a "final" release.

Deploying a Release Candidate follows the same steps as "Public" releases, except:

- Add "-rc.x" to the version number.  
   e.g., v1.2.3-rc.1; v1.2.3-rc.2 for the second release candidate, etc.
- Untick "Set as the latest release" and tick the "Set as a pre-release" option.
- Ignore Step 8 through 10; DO NOT trigger the `JuliaRegistrator` bot.

