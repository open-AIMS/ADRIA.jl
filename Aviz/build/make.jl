"""Build script that produces an app or sysimage for the ADRIA front-end"""

using Pkg, PackageCompiler

here = @__DIR__
cd(here)

if "app" in ARGS
    @info "Compiling application..."

    cd("../..")

    create_app(
        "Aviz", "adria_aviz"; include_lazy_artifacts=true, force=true, incremental=true
    )
    exit()
end

if "sysimage" in ARGS
    @info "Creating sysimage..."

    sysimage_fn = "app_image.dll"
    if "dev" in ARGS
        @info "Adding dev packages"
        dev_pkgs = ["Revise", "Infiltrator", "Plots"]
        Pkg.add(dev_pkgs)

        sysimage_fn = "app_image_dev.dll"
    end

    project_deps = collect(keys(Pkg.project().dependencies))
    try
        create_sysimage(
            project_deps;
            sysimage_path=sysimage_fn,
            precompile_execution_file="precompile_script.jl"
        )
    catch
        @info "Sysimage build failed..."
    end

    if "dev" in ARGS
        @info "Removing dev packages from project spec"
        Pkg.rm(dev_pkgs)
    end
end

# Copy the app_image.dll or app_image_dev.dll file created into your project root or development folder (e.g., `sandbox`)

# A julia session with the sysimage can be started with:
# julia --project=. -J app_image.dll

exit()
