using ADRIA.DataFrames: colmetadata, colmetadata!, Not, nrow, insertcols!

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

@testset "ResultSet inputs colmetadata" begin
    rs = TEST_RS

    # Pick a second known column (besides :guided) that is both a real model_spec
    # factor and present in rs.inputs, so the test doesn't rely on a guessed name.
    other_col = first(
        f for f in rs.model_spec.fieldname if
              f != :guided && hasproperty(rs.inputs, f)
    )

    @testset "metadata present on known columns" begin
        for col in (:guided, other_col)
            ptype = colmetadata(rs.inputs, col, "ptype", "MISSING")
            label = colmetadata(rs.inputs, col, "label", "MISSING")

            @test ptype != "MISSING"
            @test label != "MISSING"

            spec_row = rs.model_spec[rs.model_spec.fieldname .== col, :]
            @test ptype == spec_row[1, :ptype]
            @test label == spec_row[1, :name]
        end
    end

    @testset "metadata survives feature_set-style DataFrame operations" begin
        scens = copy(rs.inputs)
        @test colmetadata(scens, :guided, "ptype", "MISSING") != "MISSING"
        @test colmetadata(scens, :guided, "label", "MISSING") != "MISSING"

        # Boolean row masking
        mask = trues(nrow(scens))
        scens2 = scens[mask, :]
        @test colmetadata(scens2, :guided, "ptype", "MISSING") != "MISSING"
        @test colmetadata(scens2, :guided, "label", "MISSING") != "MISSING"

        # Not() column selection (drop a real droppable column, keep :guided)
        scens3 = scens[:, Not(:dhw_scenario)]
        @test colmetadata(scens3, :guided, "ptype", "MISSING") != "MISSING"
        @test colmetadata(scens3, :guided, "label", "MISSING") != "MISSING"
    end

    @testset "metadata does NOT auto-appear on genuinely new columns" begin
        scens = copy(rs.inputs)

        # insertcols! with a brand new column
        insertcols!(scens, :dummy_new_col => 1.0)
        @test colmetadata(scens, :dummy_new_col, "ptype", "MISSING") == "MISSING"
        @test colmetadata(scens, :dummy_new_col, "label", "MISSING") == "MISSING"

        # Arithmetic derived column also does not inherit metadata from its source
        scens[!, :guided_plus_one] = scens.guided .+ 1
        @test colmetadata(scens, :guided_plus_one, "ptype", "MISSING") == "MISSING"
        @test colmetadata(scens, :guided_plus_one, "label", "MISSING") == "MISSING"

        # Original source column metadata is untouched
        @test colmetadata(scens, :guided, "ptype", "MISSING") != "MISSING"
    end
end
