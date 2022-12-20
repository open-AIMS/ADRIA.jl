# Release Guide

## Public "Final" Releases

1. Run tests locally and ensure all pass.

2. Ensure all version numbers have been updated (check Project.toml file)

3. Submit PR from development branch into `main` and request code review/approval

4. Once PR is merged into main, go to the [releases page](https://github.com/open-AIMS/ADRIA.jl/releases) and draft a new release

5. Under "choose a tag" create a new tag "on publish"
   Note version numbers should follow Semantic Versioning: https://semver.org/

6. Click the "Generate release notes" button (top-right of textbox). 
   Under "Whats new" add a short description of the major changes.
   Explicitly note any major breaking changes (i.e., anything that results obtained with previous versions of ADRIA incompatible)

7. Release title should match the version number.

Click "Publish release" (green button at the bottom)


## Development Release

Development releases provide users with the most recent "working" version of ADRIA and may still have some known bugs.
It provides users a chance to try new features and/or provide feedback before a public release.

Deploying a Development Release follows the same steps as "Public" releases, except:

1. Add "-dev.x" to the version number.
   e.g., v1.2.3-dev.1; v1.2.3-dev.2 for the second development release, etc.

2. Untick the "Set as the latest release" option and tick the "Set as a pre-release" option.


## Release Candidates

Release candidates are releases that are not yet "final" but are close to it. Release candidates provide a "last chance" opportunity
for users to report bugs prior to a "final" release.

Deploying a Release Candidate follows the same steps as "Public" releases, except:

1. Add "-rc.x" to the version number.
   e.g., v1.2.3-rc.1; v1.2.3-rc.2 for the second release candidate, etc.

2. Untick the "Set as the latest release" option and tick the "Set as a pre-release" option.

