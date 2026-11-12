# Regression guard for precompilation time.
#
# Background: a Flatten.jl/ModelParameters specialisation blow-up in the domain-model
# assembly made `@compile_workload` take ~3 minutes (cold JIT ~365 s for the heavy block),
# turning a fresh `Pkg.precompile()` into a 20-30 minute wait. See
# `src/ecosystem/corals/flatten_overrides.jl`.
#
# The workload hides this cost when it works, so we measure the cold JIT cost directly in
# a subprocess with ADRIA's workload disabled via a LocalPreferences.toml in a temp copy
# of the active environment.
# Healthy: ~15-40 s. Broken: 150-400 s. Budget is deliberately generous (machine variance)
# but far below the broken value. Override with ADRIA_PRECOMPILE_BUDGET_S.

using Test

@testset "precompile budget" begin
    budget = parse(Float64, get(ENV, "ADRIA_PRECOMPILE_BUDGET_S", "120"))
    script = joinpath(@__DIR__, "cold_model_assembly.jl")

    mktempdir() do envdir
        # Temporary copy of the active environment that carries the preference. (A
        # stacked load-path env is NOT honoured when ADRIA is precompiled, so the
        # LocalPreferences.toml must sit next to the active project.)
        src_project = Base.active_project()
        src_dir = dirname(src_project)
        manifest_src = Base.project_file_manifest_path(src_project)
        absolutise(txt) = replace(
            txt, r"(?m)^(\s*(?:path\s*=\s*))\"([^\"]+)\"" => function (m)
                mm = match(r"^(\s*path\s*=\s*)\"([^\"]+)\"$", m)
                p = expanduser(mm.captures[2])
                p = isabspath(p) ? p : abspath(joinpath(src_dir, p))
                return string(mm.captures[1], "\"", p, "\"")
            end
        )
        write(joinpath(envdir, "Project.toml"), absolutise(read(src_project, String)))
        if manifest_src !== nothing
            write(
                joinpath(envdir, basename(manifest_src)),
                absolutise(read(manifest_src, String))
            )
        end
        write(
            joinpath(envdir, "LocalPreferences.toml"),
            "[ADRIA]\nprecompile_workload = false\n"
        )

        julia = joinpath(Sys.BINDIR, Base.julia_exename())
        project = envdir

        # Precompile without the workload first, so only JIT time is measured below.
        run(`$julia --project=$project -e "using ADRIA"`)

        out = IOBuffer()
        cmd = `$julia --project=$project $script $budget`
        proc = run(pipeline(ignorestatus(cmd); stdout=out, stderr=stderr))
        output = String(take!(out))
        m = match(r"COLD_MODEL_ASSEMBLY_SECONDS=([\d.]+)", output)
        secs = m === nothing ? NaN : parse(Float64, m.captures[1])

        @info "Cold model assembly: $(secs)s (budget $(budget)s)"
        if !(secs <= budget)
            @error """
            ################################################################
            ##  PRECOMPILATION TIME REGRESSION                            ##
            ################################################################
            Cold model assembly took $(secs)s (budget $(budget)s).
            This work runs inside ADRIA's @compile_workload, so a fresh
            `Pkg.precompile()` will now be many minutes slower for every user.
            Look at Flatten/ModelParameters overrides in
            src/ecosystem/corals/flatten_overrides.jl and any recently added
            Model components / Param types.
            ################################################################
            """
        end
        @test success(proc)
        @test secs <= budget
    end
end
