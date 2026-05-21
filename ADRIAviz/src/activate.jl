const _BACKEND_UUIDS = Dict{String,Base.UUID}(
    "WGLMakie" => Base.UUID("276b4fcb-3e11-5398-bf8b-a0c2d153d008"),
    "GLMakie" => Base.UUID("e9467ef8-e4e7-5192-8a1a-b1aee30e663a"),
    "CairoMakie" => Base.UUID("13f3f980-e62b-5c42-98c6-ff1f3baf88f0")
)

"""
    activate(backend::AbstractString="WGLMakie")

Load `backend` together with GeoMakie and GraphMakie, triggering the
`ADRIAvizMakieExt` extension so that all `ADRIA.viz.*` methods become available.

Supported backends: `"WGLMakie"`, `"GLMakie"`, `"CairoMakie"`.

The packages must already be installed in the active environment.
If they are not, install them first:

```julia-repl
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
