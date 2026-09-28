using Test
using TOML

include(joinpath(@__DIR__, "calibration_run.jl"))

const REPO_ROOT_FOR_TEST = normpath(joinpath(@__DIR__, "..", ".."))

@testset "Calibration run metadata" begin
    mktempdir(joinpath(REPO_ROOT_FOR_TEST, "sandbox", "calibration")) do output_dir
        config = CalibrationRunConfig(
            mode="FOCUSED", seed=42, max_steps=3, max_evals=2, max_time=0.0,
            method=:random_search,
            run_id="metadata_test", resume=false, output_dir=output_dir
        )
        metadata_path = joinpath(output_dir, "run_metadata.toml")
        input_path = joinpath(REPO_ROOT_FOR_TEST, "sandbox", "data", "reef_cots.csv")
        write_run_metadata(
            metadata_path, config, REPO_ROOT_FOR_TEST,
            ["a_F", "a_S"], [(0.1, 2.0), (0.01, 0.9)], [input_path]
        )
        metadata = TOML.parsefile(metadata_path)
        @test metadata["run"]["seed"] == 42
        @test metadata["search"]["parameter_names"] == ["a_F", "a_S"]
        @test metadata["inputs"][input_path] == file_sha256(input_path)

        complete_run_metadata(metadata_path; best_loss=1.25, evaluations=4)
        completed = TOML.parsefile(metadata_path)
        @test completed["result"]["best_loss"] == 1.25
        @test completed["result"]["logged_evaluations"] == 4
    end
end
