using Test
using ADRIA.Distributions
using ADRIA.decision.JMcDM
using ADRIA.decision: subtypes


@testset "Validate included MCDA methods" begin
    """
    Identifies MCDA methods that pass a simple test to inform whether they should be included
    in ADRIA.

    - A [10 ⋅ 5] decision matrix is set up (10 alternates, 5 criteria).
    - Only criteria 2 and 3 are considered, where Criteria 3 is given maximum weight.
      The others are given weights of zero.
    - The desired ranks 10 to 1 based on the values for Criteria 3
    - The desired ranks should be flipped when the directionality of optimization is reversed
      (e.g., minimize to maximize)
    - Each method is tested to ensure the known ranks are produced in both directions.

    MCDA methods that are included for use with ADRIA are those that pass the above.

    This validation accepts methods that were shown to produce nonsensical results
    when considering only a single criterion, but in practice such cases should not
    eventuate (see notes below).

    Notes:
        The MooraMethod only works if there are combinations of min and max directions.

    Listed below are methods which failed a separate simple test where only heat stress (DHW)
    is considered and the method selected the hottest locations (marked with "**")
    or produced a nonsensical selection (e.g., ranks following location order such as
    1, 2, 3, 4, 5; marked with "*^").

        - ** JMcDM.Topsis.TopsisMethod,
        - ** JMcDM.ARAS.ArasMethod,
        - *^ JMcDM.COCOSO.CocosoMethod,
        - ** JMcDM.CODAS.CodasMethod,
        - ** JMcDM.EDAS.EdasMethod,
        - ** JMcDM.GREY.GreyMethod,
        - ** JMcDM.MABAC.MabacMethod,
        - ** JMcDM.MARCOS.MarcosMethod,
        - *^ JMcDM.MOORA.MooraMethod,
        - ** JMcDM.SAW.SawMethod,
        - ** JMcDM.WASPAS.WaspasMethod,
        - ** JMcDM.WPM.WPMMethod

    Current identified methods (as of 2024-03-17) are:

    - CoCoSo
    - Mairca
    - Moora
    - PIV
    - VIKOR

    The PSI method passes the simple test based on heat stress mentioned above, but does not
    pass the assessment applied in this test.

    Moora and CoCoSo also fail the heat stress test but is included based on the assessment
    here as the conditions in which only a single criteria is examined should not occur in
    practice.
    """

    mcda_methods = subtypes(MCDMMethod)

    # Set up 10 alternates, 5 criteria
    dm = rand(10, 5)

    # Only two valid criteria (criteria 2 and 3)
    max_dir_result = collect(1:10)
    min_dir_result = reverse(max_dir_result)
    dm[:, 2] .= max_dir_result
    dm[:, 3] .= min_dir_result

    # Set weights (only really consider 3rd criteria, with some influence from 2nd)
    # Known ranks ordered by 10 to 1 (when minimizing)
    # i.e., Location 10 should be rank 1, and Location 1 should be rank 10
    weights = [0.0, 0.25, 1.0, 0.0, 0.0]
    weights = normalize(weights)

    # Define direction of interest for each criteria
    min_directions = [minimum, maximum, minimum, minimum, maximum]
    flip_directions = [maximum, minimum, maximum, maximum, minimum]

    # Results of:
    # - Scores are not all NaN
    # - Expected ranks when minimizing Criteria 3
    # - Expected ranks when maximizing Criteria 3
    test_results = zeros(length(mcda_methods), 3)

    for (i, method) in enumerate(mcda_methods)
        local res
        try
            res = mcdm(MCDMSetting(dm, weights, min_directions), method())
        catch err
            # @info "$(method) failed or not supported by JMcDM yet."
            continue
        end

        try
            res.scores
        catch err
            # @info "$(method) does not support a scoring approach"
            continue
        end

        # Test that results are not all NaN
        test_results[i, 1] = !all(isnan.(res.scores))
        test_results[i, 2] = all(
            sortperm(res.scores; rev=true) .== sortperm(min_dir_result; rev=true)
        )

        # Repeat, but switch desired direction
        try
            res = mcdm(MCDMSetting(dm, weights, flip_directions), method())
        catch
            # @info "$(method) requires at least one criteria to be minimized"
            continue
        end
        test_results[i, 3] = all(sortperm(res.scores) .== sortperm(max_dir_result))
    end

    valid_method_idx = map(x -> all(x .> 0.0), eachrow(test_results))
    valid_methods = mcda_methods[valid_method_idx]

    msg = """
    The number/order of valid MCDA methods has changed!
    Either JMcDM has updated to change the order, methods have been updated, or
    new methods have been added.

    If methods have been updated/added, these should be re-evaluated.
    """

    @test ADRIA.mcda_methods() == valid_methods || msg
end
