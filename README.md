# ADRIAjl

Julia implementation of ADRIA (in progress).

In the short-term, functionality may leverage the existing MATLAB version and so a working copy of MATLAB and ADRIA is recommended.
Interaction is handled via the MATLAB.jl interface (See intro here: https://github.com/JuliaInterop/MATLAB.jl).

The `example.ipynb` in the `examples` directory showcases an example use of MATLAB.jl


## Development setup

```bash
# Start julia specifying the current directory as the project
$ julia --project=.

# Instantiate project. Sets up project packages. Only need to do this once.
julia> ]instantiate
```