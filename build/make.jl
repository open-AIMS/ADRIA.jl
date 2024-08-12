using Pkg, PackageCompiler

# Removed ADRIA from sysimage project so line below is not needed (but could be useful later).
# project_deps = filter!(i -> i != "ADRIA", project_deps)  # Remove ADRIA (we don't want to make this static for dev purposes)
sysimage_fn = "ADRIA_sysimage.dll"
if "dev" in ARGS
    @info "Adding dev packages"
    dev_pkgs = ["Revise", "Infiltrator", "Statistics", "Plots", "GR", "GR_jll"]
    Pkg.add(dev_pkgs)

    sysimage_fn = "ADRIA_sysimage_dev.dll"
end

project_deps = collect(keys(Pkg.project().dependencies))
create_sysimage(
    project_deps;
    sysimage_path=sysimage_fn,
    precompile_execution_file="precompile_script.jl"
)

if "dev" in ARGS
    @info "Removing dev packages from project spec"
    Pkg.rm(dev_pkgs)
end

# Copy the ADRIA_sysimage.dll file created with the above into your development folder (e.g., `sandbox`)

# A julia session with the sysimage can be started with:
# julia --project=. -J ADRIA_sysimage.dll

exit()
