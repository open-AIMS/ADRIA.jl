Script to create a Julia system image (a "sysimage") that holds precompiled versions of package dependencies.
Helps reduce start up time.

## Usage

Instructions here assume you are in this `build` directory.

The very first time (and only the first time!):

```bash
$ julia --project=.
julia> ] add PackageCompiler
julia> exit()
```

Then:

```bash
$ julia --project=. make.jl
```

To create a development sysimage which includes Revise, Infiltrator and other helpful development packages:

```bash
$ julia --project=. make.jl dev
```

The process will create an ADRIA_sysimage.dll file (or ADRIA_sysimage_dev.dll if the `dev` flag was used).
Copy this file to the ADRIA project folder (or your "sandbox" directory if using the dev sysimage).

When starting Julia, use the `-J` flag and point to the sysimage location, for example:

```bash
# From project directory
$ julia --project=. -J ADRIA_sysimage.dll
```

This will reduce startup time from minutes to seconds.