# Keep the result stores written by tests out of the repository.
#
# `ADRIA.run_scenarios` (and `ADRIA.combine_results`) write a result store under
# `ENV["ADRIA_OUTPUT_DIR"]`. `ADRIA.setup()` sets that from `./config.toml`, or to
# `./Outputs` when there is none, both relative to the current directory. And
# `run_scenarios` itself calls `setup()`, which on its first call in a process overwrites any
# ADRIA_OUTPUT_DIR set before it. So `setup()` has to run before the variable is redirected.

using ADRIA

"""
    with_temp_output_dir(f)

Call `f()` with result stores written to a new temporary directory, restoring the previous
`ADRIA_OUTPUT_DIR` afterwards. The directory is kept until the Julia process exits, so
result sets returned by `f` can still be read.
"""
function with_temp_output_dir(f)
    ADRIA.setup()
    return withenv(f, "ADRIA_OUTPUT_DIR" => mktempdir())
end
