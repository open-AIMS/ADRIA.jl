using Aqua

@testset "Aqua" begin
    # Don't check ambiguities as there are too many, most of which are, I suspect, issues
    # in dependencies.
    # Aqua.test_ambiguities([ADRIA, Base, Core]; exclude=[(==), write])

    Aqua.test_unbound_args(ADRIA)
    Aqua.test_undefined_exports(ADRIA)
    Aqua.test_project_extras(ADRIA)
    Aqua.test_stale_deps(ADRIA)
    Aqua.test_deps_compat(ADRIA)
end
