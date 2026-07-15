using Test
using DataFrames
using ADRIA
using ADRIAanalysis

# `feature_set(rs::ResultSet)` needs a real ResultSet (dhw_stats, ranks, seed_log,
# mc_log, etc. are all populated from an actual simulation run) -- there is no
# lightweight mock for this in ADRIAanalysis/test. Reuse the same gated-fixture
# pattern as `integration_cluster_rules.jl`: run a small scenario set against the
# on-disk `Test_domain` fixture, skipping (with a @warn, not a failure) if that
# domain data isn't available in this environment.
domain_path = get(
    ENV, "ADRIA_TEST_DOMAIN_PATH",
    joinpath(@__DIR__, "..", "..", "ADRIA", "test", "data", "Test_domain")
)

if !isdir(domain_path)
    @warn "No domain path (ADRIA_TEST_DOMAIN_PATH) -- feature_set metadata tests skipped"
else
    @testset "feature_set colmetadata tagging" begin
        dom = ADRIA.load_domain(domain_path, "45")
        scens = ADRIA.sample(dom, 8)
        rs = ADRIA.run_scenarios(dom, scens, "45")

        fs = ADRIAanalysis.feature_set(rs)
        @test fs isa DataFrame

        fs_cols = Symbol.(names(fs))

        # Explicitly-tagged derived columns from feature_set.jl. Some may be
        # filtered out by `_filter_constants` in a small synthetic/short run
        # (e.g. a single-functional-group domain, or a factor that happened
        # not to vary across the sampled scenarios) -- only assert on columns
        # that actually survive into the output.
        candidate_tagged_cols = [
            :dhw_mean, :dhw_stdev, :dhw_complexity,
            :n_loc_seed_mean, :n_loc_fog_mean, :n_loc_mc_mean,
            :depth_max,
            :seed_total_deployed_coral_M, :mc_total_deployed_coral_M
        ]
        present_tagged_cols = [c for c in candidate_tagged_cols if c in fs_cols]

        # At minimum, DHW stat columns and depth_max should always be present
        # (they don't depend on functional-group/deployment configuration).
        @test :dhw_mean in present_tagged_cols
        @test :dhw_stdev in present_tagged_cols
        @test :dhw_complexity in present_tagged_cols
        @test :depth_max in present_tagged_cols

        for col in present_tagged_cols
            @test colmetadata(fs, col, "ptype", "MISSING") != "MISSING"
            @test colmetadata(fs, col, "label", "MISSING") != "MISSING"
        end

        # At least one seed_volume_*/mc_volume_* per-functional-group column
        # should also be tagged, if any survived `_filter_constants`.
        volume_cols = [
            c for c in fs_cols if
            startswith(string(c), "seed_") || startswith(string(c), "mc_")
        ]
        volume_cols = [
            c for c in volume_cols if
            c ∉ (:seed_total_deployed_coral_M, :mc_total_deployed_coral_M)
        ]
        if !isempty(volume_cols)
            for col in volume_cols
                @test colmetadata(fs, col, "ptype", "MISSING") != "MISSING"
                @test colmetadata(fs, col, "label", "MISSING") != "MISSING"
            end
        else
            @warn "No seed_/mc_ volume columns survived _filter_constants; " *
                "skipping per-functional-group metadata assertions."
        end

        # Normalized "actual effort" columns: min-max in [0, 1], tagged, and
        # skipped for any source column that was constant across the sample.
        effort_cols = [c for c in fs_cols if endswith(string(c), "_effort")]
        if !isempty(effort_cols)
            for col in effort_cols
                @test colmetadata(fs, col, "ptype", "MISSING") != "MISSING"
                @test colmetadata(fs, col, "label", "MISSING") != "MISSING"
                vals = fs[!, col]
                @test all(0.0 .<= vals .<= 1.0)
                @test minimum(vals) == 0.0
                @test maximum(vals) == 1.0

                source_col = Symbol(replace(string(col), "_effort" => ""))
                @test source_col in fs_cols
            end
        else
            @warn "No *_effort columns survived (all source columns constant " *
                "in this small sample); skipping effort-metric assertions."
        end
    end
end

# ----------------------------------------------------------------------------
# Lightweight synthetic check: colmetadata survives the same kinds of
# column/row selection operations feature_set performs internally
# (Not(...) column drops, hcat!, row masking) -- doesn't require a ResultSet.
# ----------------------------------------------------------------------------
@testset "colmetadata survives selection operations used by feature_set" begin
    df = DataFrame(; a=rand(10), b=rand(10), c=rand(10))
    colmetadata!(df, :a, "ptype", "continuous"; style=:note)
    colmetadata!(df, :a, "label", "A label"; style=:note)

    # Symbol-based column drop (as used for :depth_offset, :dhw_scenario)
    df_dropped = df[:, Not(:b)]
    @test colmetadata(df_dropped, :a, "ptype", "MISSING") == "continuous"
    @test colmetadata(df_dropped, :a, "label", "MISSING") == "A label"

    # Row masking (as used per-stratum / per-cf_mask slicing)
    mask = vcat(trues(5), falses(5))
    df_masked = df[mask, :]
    @test colmetadata(df_masked, :a, "ptype", "MISSING") == "continuous"

    # hcat! with a second frame (as used to attach seed/mc volume columns)
    extra = DataFrame(; d=rand(10))
    colmetadata!(extra, :d, "ptype", "continuous"; style=:note)
    DataFrames.hcat!(df, extra)
    @test colmetadata(df, :a, "ptype", "MISSING") == "continuous"
    @test colmetadata(df, :d, "ptype", "MISSING") == "continuous"
end
