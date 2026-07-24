const _BACKEND_UUIDS = Dict{String,Base.UUID}(
    "WGLMakie" => Base.UUID("276b4fcb-3e11-5398-bf8b-a0c2d153d008"),
    "GLMakie" => Base.UUID("e9467ef8-e4e7-5192-8a1a-b1aee30e663a"),
    "CairoMakie" => Base.UUID("13f3f980-e62b-5c42-98c6-ff1f3baf88f0")
)

const _PLOTLY_PKG_ID = Base.PkgId(
    Base.UUID("a03496cd-edff-5a9b-9e67-9cda94a718b5"), "PlotlyBase"
)

const _KALEIDO_PKG_ID = Base.PkgId(
    Base.UUID("f2990250-8cf9-495f-b13a-cce12b45703c"), "PlotlyKaleido"
)

"""
    activate(backend::Symbol)

Convenience method. The only supported symbol is `:plotly`.

Loads `PlotlyBase` (required) and `PlotlyKaleido` (optional — enables
`ADRIA.viz.savefig`). If `PlotlyKaleido` is not installed, a warning is printed
and only the interactive backend is activated.

```julia
using ADRIA, ADRIAviz
ADRIAviz.activate("plotly")   # loads PlotlyBase + PlotlyKaleido (if installed)
```
"""
function activate(backend::Symbol)
    backend in (:plotly,) || throw(
        ArgumentError(
            "Symbol shorthand only supports :plotly. " *
            "For Makie backends use `activate(\"WGLMakie\")` etc."
        )
    )
    isnothing(Base.find_package("PlotlyBase")) && throw(
        ArgumentError(
            "PlotlyBase is not installed. Install it first:\n\n    pkg> add PlotlyBase\n"
        )
    )
    Base.require(_PLOTLY_PKG_ID)

    if isnothing(Base.find_package("PlotlyKaleido"))
        @warn "PlotlyKaleido is not installed — `ADRIA.viz.savefig` will not be available.\n" *
            "To enable static export: pkg> add PlotlyKaleido"
    else
        kaleido = Base.require(_KALEIDO_PKG_ID)
        Base.invokelatest(kaleido.start)
    end

    return nothing
end

"""
    activate(backend::AbstractString="WGLMakie")

Load `backend` together with GeoMakie and GraphMakie, triggering the
`ADRIAvizMakieExt` extension so that all `ADRIA.viz.*` methods become available.

Supported backends: `"WGLMakie"`, `"GLMakie"`, `"CairoMakie"`.

The packages must already be installed in the active environment.
If they are not, install them first:

```julia
pkg> add WGLMakie GeoMakie GraphMakie
```

# Example

```julia
using ADRIA, ADRIAviz
ADRIAviz.activate()               # defaults to WGLMakie
ADRIAviz.activate("CairoMakie")   # non-interactive / CI
```
"""
function activate(backend::AbstractString="WGLMakie")
    # Convenience alias — delegate to Symbol overload
    if lowercase(backend) == "plotly"
        return activate(:plotly)
    end

    haskey(_BACKEND_UUIDS, backend) || throw(
        ArgumentError(
            "Unknown backend \"$(backend)\". Supported: $(join(sort(collect(keys(_BACKEND_UUIDS))), ", "))"
        )
    )

    needed = [backend, "GeoMakie", "GraphMakie"]
    missing_pkgs = filter(p -> isnothing(Base.find_package(p)), needed)
    if !isempty(missing_pkgs)
        throw(
            ArgumentError(
                "The following packages are not installed: $(join(missing_pkgs, ", "))\n" *
                "Install them first:\n\n" *
                "    pkg> add $(join(needed, " "))\n"
            )
        )
    end

    _req(name, uuid) = Base.require(Base.PkgId(Base.UUID(uuid), name))

    _req(backend, string(_BACKEND_UUIDS[backend]))
    _req("GeoMakie", "db073c08-6b98-4ee5-b6a4-5efafb3259c6")
    _req("GraphMakie", "1ecd5474-83a3-4783-bb4f-06765db800d2")

    return nothing
end
